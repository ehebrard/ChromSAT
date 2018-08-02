#include <iostream>

#include "brancher.hpp"
#include "dimacs.hpp"
#include "fillin.hpp"
#include "graph.hpp"
#include "mycielski.hpp"
#include "options.hpp"
#include "prop.hpp"
#include "rewriter.hpp"
#include "sparse_dynamic_graph.hpp"
#include "statistics.hpp"
#include "utils.hpp"
#include "vcsolver.hpp"

#include <minicsp/core/cons.hpp>
#include <minicsp/core/solver.hpp>
#include <minicsp/core/utils.hpp>

enum class vertex_status : uint8_t {
    in_graph,
    low_degree_removed,
    indset_removed,
};

template <class adjacency_struct> void histogram(gc::graph<adjacency_struct>& g)
{
    std::vector<int> degrees;
    for (auto v : g.nodes) {
        degrees.push_back(g.matrix[v].size());
    }
    std::sort(begin(degrees), end(degrees));
    int psize = degrees.size() / 10;
    int pstart{0};
    while (static_cast<size_t>(pstart) < degrees.size()) {
        int pend
            = std::min(static_cast<size_t>(pstart + psize), degrees.size() - 1);
        while (static_cast<size_t>(pend) < degrees.size() - 1
            && degrees[pend + 1] == degrees[pend])
            ++pend;
        std::cout << (pend - pstart + 1) << " vertices: degree "
                  << degrees[pstart] << "-" << degrees[pend] << "\n";
        pstart = pend + 1;
    }

    double sum_degrees(0);
    for (auto d : degrees) {
        sum_degrees += d;
    }
    std::cout << "Density : "
              << sum_degrees / ((double)(g.size()) * (double)(g.size() - 1))
            * 100
              << "%" << std::endl;
}

template< class adjacency_struct >
struct graph_reduction {
    const gc::graph<adjacency_struct>& g;
    const gc::statistics& statistics;
    std::vector<int> removed_vertices;
    std::vector<gc::indset_constraint> constraints;
    std::vector<vertex_status> status;
    gc::bitset nodeset;
    gc::bitset util_set;

    explicit graph_reduction(
        const gc::graph<adjacency_struct>& g, const gc::statistics& statistics)
        : g(g)
        , statistics(statistics)
        , status(g.capacity(), vertex_status::in_graph)
        , nodeset(0, g.capacity(), gc::bitset::empt)
        , util_set(0, g.capacity(), gc::bitset::empt)
    {
    }

    int extend_solution(std::vector<int>& col)
    {
        int maxc{0};
        nodeset.copy(g.nodeset);
        for (auto v : g.nodes)
            maxc = std::max(maxc, col[v]);
        for (auto i = removed_vertices.rbegin(), iend = removed_vertices.rend();
             i != iend; ++i) {
            auto v = *i;
            assert(status[v] != vertex_status::in_graph);
            util_set.clear();
            for (auto u : g.matrix[v]) {
                if (!nodeset.fast_contain(u))
                    continue;
                util_set.fast_add(col[u]);
            }
            for (int q = 0; q != g.capacity(); ++q) {
                if (util_set.fast_contain(q))
                    continue;
                assert(q < statistics.best_ub);
                maxc = std::max(maxc, q);
                col[v] = q;
                break;
            }
            nodeset.fast_add(v);
        }
        return maxc + 1;
    }
};

template< class adjacency_struct >
struct gc_model {
    int lb{0}, ub{-1};
    const gc::options& options;
    gc::statistics& statistics;

    std::vector<int>
        vertex_map; // maps vertices of the original graph to vertices of g

    gc::graph<adjacency_struct>& original;
    graph_reduction<adjacency_struct> reduction;
    gc::dense_graph g;

    boost::optional<std::vector<std::pair<int, int>>> fillin;

    minicsp::Solver s;
    gc::varmap vars;
    gc::cons_base* cons;
    std::vector<minicsp::cspvar> xvars;
    gc::rewriter rewriter;

    std::unique_ptr<gc::Brancher> brancher;

    // std::vector<edge> var_map;

    gc::varmap create_chord_vars()
    {
        gc::minfill_buffer<gc::dense_graph> mb{g};
        mb.minfill();
        std::cout << "[modeling] "<< mb.fillin.size()
                  << " edges in minfill, width = " << mb.width << "\n";
        fillin = std::move(mb.fillin);

        using std::begin;
        using std::end;
        gc::varmap vars(begin(g.nodes), end(g.nodes));
        vars.vars.resize(g.size());
        for (auto e : *fillin) {
            int i = e.first, j = e.second;
            Var v = s.newVar();
            vars.vars[i][j] = v;
            vars.vars[j][i] = v;
            if (options.trace) {
                using namespace std::string_literals;
                using std::to_string;
                auto n = "e"s + to_string(i) + "-"s + to_string(j);
                s.setVarName(vars.vars[i][j], n);
            }
        }

        std::cout << "[modeling] created " << s.nVars() << " chord variables\n\n";

        return vars;
    }

    gc::varmap create_all_vars()
    {
        using std::begin;
        using std::end;
        gc::varmap vars(begin(g.nodes), end(g.nodes));
        vars.vars.resize(g.size());
        for (auto i : g.nodes) {
             for (auto j : g.nodes) {
                if (j < i)
                    continue;
                if (g.matrix[i].fast_contain(j))
                    continue;
                vars.vars[i][j] = s.newVar();
                vars.vars[j][i] = vars.vars[i][j];
                if (options.trace) {
                    using namespace std::string_literals;
                    using std::to_string;
                    auto n = "e"s + to_string(i) + "-"s + to_string(j);
                    s.setVarName(vars.vars[i][j], n);
                }
            }
        }

        std::cout << "[modeling] created " << s.nVars()
                  << " classic variables\n\n";

        return vars;
    }

    gc::varmap create_vars()
    {
        if (options.fillin)
            return create_chord_vars();
        else
            return create_all_vars();
    }


    graph_reduction<adjacency_struct> degeneracy_peeling(
        gc::graph<adjacency_struct>& g, std::pair<int, int> bounds,
        bool myciel = false)
    {
        graph_reduction<adjacency_struct> gr(g, statistics);
        if (options.preprocessing == gc::options::NO_PREPROCESSING)
            return gr;

        std::cout << "[preprocessing] start peeling\n";

        lb = bounds.first;
        ub = bounds.second;

        gc::clique_finder<adjacency_struct> cf{g, std::min(100, g.size())};
        gc::mycielskan_subgraph_finder<adjacency_struct> mf(g, cf, false);
        gc::degeneracy_finder<gc::graph<adjacency_struct>> df{g};

        adjacency_struct toremove;
        toremove.initialise(0, g.capacity(), gc::bitset::empt);

        int plb{lb-1};

        do {

            if (lb == plb) {
                break;
            }
            plb = lb;
            if (g.size() == 0) {
                break;
            }

            // histogram(g);

            cf.clear();
            df.clear();
            // mf.clear();

            std::cout << "[preprocessing] compute degeneracy ordering\n";
            int degeneracy = 0;
            df.degeneracy_ordering();

            // std::cout << df.order.size() << " / " << g.size() << std::endl;

            // for(auto v : df.order) {
            //      std::cout << " " << df.degrees[v];
            // }
            // std::cout << std::endl;

            std::vector<int>
                reverse; // TODO change that by using iterators in clique finder
            for (auto rit = df.order.rbegin(); rit != df.order.rend(); ++rit) {
                if (df.degrees[*rit] > degeneracy) {
                    degeneracy = df.degrees[*rit];
                }
                reverse.push_back(*rit);
            }

            if (ub > degeneracy + 1) {
                ub = degeneracy + 1;
                statistics.notify_ub(ub);
                statistics.display(std::cout);
            }

            std::cout << "[preprocessing] compute lower bound\n";

            lb = cf.find_cliques(reverse);
            if (false and g.size() < 1000)
                lb = mf.improve_cliques_larger_than(lb);
            if (lb < plb)
                lb = plb;
            statistics.notify_lb(lb);
            statistics.display(std::cout);

            // std::cout << ", lb = " << lb << std::endl;

            std::cout << "[preprocessing] remove low degree nodes\n";

            bool removal = false;
            for (auto v : df.order) {
                if (df.degrees[v] >= lb)
                    break;
                toremove.add(v);
                gr.removed_vertices.push_back(v);
                removal = true;
            }

            if (removal) {
                toremove.canonize();
                g.remove(toremove);
                for (auto u : toremove) {
                    gr.status[u] = vertex_status::low_degree_removed;
                }

                toremove.clear();
                statistics.notify_removals(g.size());
                statistics.display(std::cout);
            }

        } while (true);

        return gr;
    }


    graph_reduction<adjacency_struct> core_reduction(
        gc::graph<adjacency_struct>& g, std::pair<int, int> bounds,
        bool myciel = false)
    {
        graph_reduction<adjacency_struct> gr(g, statistics);
        if (options.preprocessing == gc::options::NO_PREPROCESSING)
            return gr;

        // std::cout << "CORE REDUCTION: " << g.size() << "(" << (int*)(&g) << ")"
        //           << std::endl;

        lb = bounds.first;
        ub = bounds.second;
        int hlb{0};
        gc::clique_finder<adjacency_struct> cf(g);
        gc::mycielskan_subgraph_finder<adjacency_struct> mf(g, cf, false);
        gc::degeneracy_finder<gc::graph<adjacency_struct>> df{g};

        gc::bitset forbidden(0, g.capacity(), gc::bitset::empt);
        gc::bitset util_set(0, g.capacity(), gc::bitset::empt);
        adjacency_struct removedv(0, g.capacity(), gc::bitset::empt);
        std::vector<int> toremove;
        bool removed{false};
        int niteration{0};
        do {
            ++niteration;
            removed = false;
            auto sol{gc::brelaz_color(g)};
            for (auto u : g.nodes)
                for (auto v : g.matrix[u])
                    assert(sol[u] != sol[v]);
            int hub{*max_element(begin(sol), end(sol)) + 1};

            // df.degeneracy_ordering();
            //             int hub{*max_element(
            //                 begin(df.degrees), end(df.degrees))};
            if (ub < 0 || (hub < ub && hub >= lb)) {
                ub = hub;
                statistics.notify_ub(ub);
            }

            hlb = cf.find_cliques(g.nodes);
            if (myciel)
                hlb = mf.improve_cliques_larger_than(lb);

            if (hlb > lb) {
                lb = hlb;
                statistics.notify_lb(lb);
            }
            statistics.display(std::cout);

            forbidden.clear();
            toremove.clear();
            for (auto u : g.nodes) {
                if (forbidden.fast_contain(u))
                    continue;
                util_set.copy(g.matrix[u]);
                util_set.intersect_with(g.nodeset);
                if (util_set.size() >= static_cast<size_t>(lb))
                    continue;
                removed = true;
                removedv.fast_add(u);
                // ++statistics.num_vertex_removals;
                toremove.push_back(u);
                gr.removed_vertices.push_back(u);
                gr.status[u] = vertex_status::low_degree_removed;
                forbidden.union_with(g.matrix[u]);
            }
            for (auto u : toremove) {
                g.nodes.remove(u);
                g.nodeset.remove(u);
            }
        } while (removed);
        if (removedv.size() > 0) {
            for (auto v : g.nodes) {
                g.matrix[v].setminus_with(removedv);
                g.origmatrix[v].setminus_with(removedv);
            }
            statistics.notify_removals(g.size());
            statistics.display(std::cout);
        }

        return gr;
    }

    // template< class adjacency_struct >
    void find_is_constraints(gc::graph<adjacency_struct>& g, graph_reduction<adjacency_struct>& gr)
    {
        gc::degeneracy_vc_solver<gc::graph<adjacency_struct>> vc(g);
        auto bs = vc.find_is();
        std::cout << "IS size = " << bs.size() << "\n";
        for (auto v : bs) {
            gr.removed_vertices.push_back(v);
            gr.status[v] = vertex_status::indset_removed;
            gr.constraints.emplace_back(gc::indset_constraint{g.matrix[v], v});
            g.nodes.remove(v);
            g.nodeset.remove(v);
        }
        for (auto v : g.nodes)
            g.matrix[v].intersect_with(g.nodeset);

        statistics.notify_removals(g.size());
        statistics.display(std::cout);
    }

    // template< class adjacency_struct >
    graph_reduction<adjacency_struct> preprocess(
        gc::graph<adjacency_struct>& g, std::pair<int, int> bounds, bool myciel = false)
    {

        // std::cout << "PREPROCESS: " << g.size() << "(" << (int*)(&g) << ")"
        //           << std::endl;

        auto gr{core_reduction(g, bounds, myciel)};
        if (g.size() > 0 and options.indset_constraints) {
            std::cout << "[preprocessing] search for indset constraints\n";
            find_is_constraints(g, gr);
            for (auto& c : gr.constraints) {
                std::cout << "constraint " << c.vs << "\n";
            }
        }
        std::cout << "[preprocessing] finished at " << minicsp::cpuTime() << "\n\n";
        return gr;
    }

    graph_reduction<adjacency_struct> qpreprocess(
        gc::graph<adjacency_struct>& g, std::pair<int, int> bounds,
        bool myciel = false)
    {
        auto gr{degeneracy_peeling(g, bounds, myciel)};

        if (g.size() > 0) {
            std::cout << "[preprocessing] launch dsatur\n";

            gc::dyngraph dg(g);
            // dg.sort(false);
            gc::coloring col;

            int maxiter = 1, iter = 0;
            do {
                col.brelaz_color(dg, 0);

                auto ncol{
                    *std::max_element(begin(col.color), end(col.color)) + 1};

                if (ub > ncol) {
                    ub = ncol;
                    statistics.notify_ub(ub);
                    statistics.display(std::cout);
                }

                if (++iter >= maxiter) {
                    // std::cout << "stop because of limit\n";
                    break;
                }

                dg.undo();
                col.clear();
            } while (lb < ub);

            if (g.size() > 0 and options.indset_constraints)
                find_is_constraints(g, gr);
        }
        std::cout << "[preprocessing] finished at " << minicsp::cpuTime() << "\n";
        return gr;
    }

    // template< class graph_struct >
    gc_model(gc::graph<adjacency_struct>& ig, const gc::options& options,
        gc::statistics& statistics, std::pair<int, int> bounds)
        : options(options)
        , statistics(statistics)
        , vertex_map(ig.capacity())
        , original(ig)
        , reduction(qpreprocess(ig, bounds))
        // , reduction(preprocess(ig, bounds))
        , g(ig, vertex_map)
        , vars(create_vars())
        , cons(gc::post_gc_constraint(
              s, g, fillin, vars, reduction.constraints, options, statistics))
        , rewriter(s, g, cons, vars, xvars)
    {
        setup_signal_handlers(&s);
        s.trace = options.trace;
        s.polarity_mode = options.polarity;

        if (options.learning == gc::options::NO_LEARNING)
            s.learning = false;

        auto plb = bounds.first;
        auto pub = bounds.second;

        lb = std::max(lb, plb);
        if (ub < 0 || pub < ub)
            ub = pub;

        if (cons) {

            cons->bestlb = std::max(lb, cons->bestlb);
            cons->ub = std::min(ub, cons->ub);

            if (options.xvars) {
                xvars = s.newCSPVarArray(g.capacity(), 0, cons->ub - 2);
                for (size_t i = 0; i != xvars.size(); ++i) {
                    if (!g.nodes.contain(i))
                        continue;
                    for (size_t j = i + 1; j != xvars.size(); ++j) {
                        if (!g.nodes.contain(j))
                            continue;
                        if (g.matrix[i].fast_contain(j))
                            minicsp::post_neq(s, xvars[i], xvars[j], 0);
                        else
                            minicsp::post_eq_re(s, xvars[i], xvars[j], 0,
                                minicsp::Lit(vars[i][j]));
                    }
                }

                // rewrite clauses to not use x literals
                s.use_clause_callback([this](vec<minicsp::Lit>& clause,
                    int btlvl) { return rewriter.rewrite(clause, btlvl); });
            }

            auto sum = [](int x, int y) { return x + y; };
            auto prod = [](int x, int y) { return x * y; };

            switch (options.branching) {
            case gc::options::VSIDS:
                if (options.branching_low_degree) {
                    brancher = std::make_unique<gc::VSIDSBrancher>(
                        s, g, vars, xvars, *cons, options);
                    brancher->use();
                } else
                    s.varbranch = minicsp::VAR_VSIDS;
                break;
            case gc::options::VSIDS_GUIDED:
                if (options.branching_low_degree) {
                    brancher = std::make_unique<gc::VSIDSBrancher>(
                        s, g, vars, xvars, *cons, options);
                    brancher->use();
                } else
                    s.varbranch = minicsp::VAR_VSIDS;
                s.phase_saving = false;
                s.solution_phase_saving = true;
                break;
            case gc::options::VSIDS_PHASED:
                brancher = std::make_unique<gc::VSIDSPhaseBrancher>(
                    s, g, vars, xvars, *cons, options, -1, -1);
                brancher->use();
                break;
            case gc::options::VSIDS_CLIQUE:
                brancher = std::make_unique<gc::VSIDSCliqueBrancher>(
                    s, g, vars, xvars, *cons, options);
                break;
            case gc::options::VSIDS_COLORS_POSITIVE:
                if (!options.xvars) {
                    std::cout << "VSIDS_COLORS_POSITIVE needs --xvars\n";
                    exit(1);
                }
                brancher = std::make_unique<gc::VSIDSColorBrancher>(
                    s, g, vars, xvars, *cons, options);
                brancher->use();
                break;
            case gc::options::BRELAZ:
                brancher = std::make_unique<gc::BrelazBrancher>(
                    s, g, vars, xvars, *cons, options);
                brancher->use();
                break;
            case gc::options::PARTITION_SUM:
                brancher = gc::make_partition_brancher<-1, -1>(
                    s, g, vars, xvars, *cons, options, sum);
                brancher->use();
                break;
            case gc::options::PARTITION_PRODUCT:
                brancher = gc::make_partition_brancher<-1, -1>(
                    s, g, vars, xvars, *cons, options, prod);
                brancher->use();
                break;
            case gc::options::DEGREE_SUM:
                brancher = gc::make_degree_brancher<-1, -1>(
                    s, g, vars, xvars, *cons, options, sum);
                brancher->use();
                break;
            case gc::options::DEGREE_PRODUCT:
                brancher = gc::make_degree_brancher<-1, -1>(
                    s, g, vars, xvars, *cons, options, prod);
                brancher->use();
                break;
            case gc::options::DEGREE_UNION:
                brancher = std::make_unique<gc::DegreeUnionBrancher<-1, -1>>(
                    s, g, vars, xvars, *cons, options);
                brancher->use();
                break;
            case gc::options::PARTITION_SUM_DYN:
                brancher = gc::make_partition_brancher<2, 3>(
                    s, g, vars, xvars, *cons, options, sum);
                brancher->use();
                break;
            case gc::options::PARTITION_PRODUCT_DYN:
                brancher = gc::make_partition_brancher<2, 3>(
                    s, g, vars, xvars, *cons, options, prod);
                brancher->use();
                break;
            case gc::options::DEGREE_SUM_DYN:
                brancher = gc::make_degree_brancher<2, 3>(
                    s, g, vars, xvars, *cons, options, sum);
                brancher->use();
                break;
            case gc::options::DEGREE_PRODUCT_DYN:
                brancher = gc::make_degree_brancher<2, 3>(
                    s, g, vars, xvars, *cons, options, prod);
                brancher->use();
                break;
            case gc::options::DEGREE_UNION_DYN:
                brancher = std::make_unique<gc::DegreeUnionBrancher<2, 3>>(
                    s, g, vars, xvars, *cons, options);
                brancher->use();
                break;
            }
        }
    }

    std::vector<int> get_solution()
    {
        std::vector<int> col(original.capacity());
        int next{0};
        for (auto u : g.nodes) {
            for (auto v : g.partition[u])
                col[original.nodes[v]] = next;
            ++next;
        }
        return col;
    }

    std::pair<int, int> solve()
    {
        using minicsp::l_False;
        using minicsp::l_True;
        using minicsp::l_Undef;

        auto llb = lb;
        auto lub = ub;

        minicsp::lbool sat{l_True};
        while (sat != l_False && lb < ub) {
            sat = s.solveBudget();
            if (sat == l_True) {
                auto col_r = get_solution();
                int solub = g.nodes.size();
                assert(solub < cons->ub);
                cons->sync_graph();
                int actualub = reduction.extend_solution(col_r);
                statistics.notify_ub(actualub);
                statistics.display(std::cout);

                cons->ub = solub;
                cons->actualub = actualub;
                if (options.xvars) {
                    for (auto v : xvars)
                        v.setmax(s, cons->ub - 2, minicsp::NO_REASON);
                }
                assert(lub >= cons->ub);
                lub = cons->ub;
            } else if (sat == l_Undef) {
                std::cout << "*** INTERRUPTED ***\n";
                break;
            } else {
                cons->bestlb = cons->ub;
                statistics.display(std::cout);
                assert(llb <= cons->bestlb);
                llb = cons->bestlb;
            }
        }
        return std::make_pair(llb,lub);
    }

    void print_stats()
    {
        assert(!cons
            or (cons->bestlb == statistics.best_lb
                   and cons->ub == statistics.best_ub));
        // assert(lb == statistics.best_lb and ub == statistics.best_ub);

        if (statistics.best_lb >= statistics.best_ub)
            std::cout << "OPTIMUM " << statistics.best_ub << "\n";
        else
            std::cout << "Best bounds [" << statistics.best_lb << ", "
                      << statistics.best_ub << "]\n";
        minicsp::printStats(s);
        statistics.display(std::cout);
        std::cout << std::endl;
    }
};

template<class adjacency_struct>
std::pair<int, int> initial_bounds(
    const gc::graph<adjacency_struct>& g, gc::statistics& stat, bool myciel = false)
{
    // gc::degeneracy_finder df{g};
    // df.degeneracy_ordering();
    // for( auto u : df.order ) {
    //  std::cout << u << "(" << g.matrix[u].size() << ") ";
    // }
    // std::cout << std::endl;

    gc::clique_finder<adjacency_struct> cf{g};
    gc::mycielskan_subgraph_finder<adjacency_struct> mf(g, cf, false);
    int lb{cf.find_cliques(g.nodes)};

    if (myciel)
        lb = mf.improve_cliques_larger_than(lb - 1);

    auto sol{gc::brelaz_color(g)};
    for (auto u : g.nodes)
        for (auto v : g.matrix[u])
            assert(sol[u] != sol[v]);

    int ub{*max_element(begin(sol), end(sol)) + 1};
    stat.notify_lb(lb);
    stat.notify_ub(ub);
    stat.display(std::cout);
    return std::make_pair(lb, ub);
}

int main(int argc, char* argv[])
{
    auto options = gc::parse(argc, argv);
    options.describe(std::cout);

    std::cout << "[reading] ";

    gc::graph<gc::vertices_vec> g;
    int num_edges = 0;
    dimacs::read_graph(options.instance_file.c_str(),
        [&](int nv, int) { g = gc::graph<gc::vertices_vec>{nv}; },
        [&](int u, int v) {
            if (u != v) {
                g.add_edge(u - 1, v - 1);
                ++num_edges;
            }
        },
        [&](int, gc::weight) {});
    g.canonize();

    g.describe(std::cout, num_edges);
    // histogram(g);

    gc::statistics statistics(g.capacity());
    if (options.preprocessing)
        statistics.update_ub = false;

    statistics.describe(std::cout);

    // std::cout << "MAIN (READ): " << g.size() << "(" << (int*)(&g) << ")"
    //           << std::endl;
    switch (options.strategy) {
    case gc::options::BNB: {
        std::pair<int, int> bounds{0, g.capacity()};
        if (options.preprocessing == gc::options::NO_PREPROCESSING) {
            std::cout << "compute init bounds\n";
            bounds = initial_bounds(
                g, statistics, options.boundalg != gc::options::CLIQUES);
        }

        gc_model<gc::vertices_vec> model(g, options, statistics, bounds);
        model.solve();
        model.print_stats();
        break;
    }
    case gc::options::BOTTOMUP: {
        statistics.update_ub = false;
        auto bounds = initial_bounds(
            g, statistics, options.boundalg != gc::options::CLIQUES);
        auto lb = bounds.first;
        auto ub = bounds.second;
        for (int i = lb; i < ub; ++i) {
            gc::graph<gc::bitset> gcopy{g};
            gc_model<gc::bitset> model(
                gcopy, options, statistics, std::make_pair(i, i + 1));
            auto ibounds = model.solve();
            auto ilb = ibounds.first;
            auto iub = ibounds.second;
            if (iub == i) {
                statistics.notify_ub(iub);
                statistics.display(std::cout);
                break;
            } else if (ilb != iub) {
                statistics.notify_lb(ilb);
                statistics.display(std::cout);
                std::cout << "INTERRUPTED\n";
            } else {
                statistics.notify_lb(ilb);
                statistics.display(std::cout);
            }
            statistics.unbinds();
        }
    } break;
    case gc::options::TOPDOWN: {
        statistics.update_lb = false;
        auto bounds = initial_bounds(
            g, statistics, options.boundalg != gc::options::CLIQUES);
        auto lb = bounds.first;
        auto ub = bounds.second;
        for (int i = ub - 1; i >= lb; --i) {
            gc::graph<gc::bitset> gcopy{g};
            gc_model<gc::bitset> model(
                gcopy, options, statistics, std::make_pair(i, i + 1));
            auto ibounds = model.solve();
            auto ilb = ibounds.first;
            auto iub = ibounds.second;
            if (ilb == i + 1) {
                statistics.notify_lb(ilb);
                statistics.display(std::cout);
                break;
            } else if (ilb != iub) {
                statistics.notify_ub(iub);
                statistics.display(std::cout);
                std::cout << "INTERRUPTED\n";
            } else {
                statistics.notify_ub(iub);
                statistics.display(std::cout);
            }
            statistics.unbinds();
        }
    } break;
    case gc::options::BOUNDS: {
        std::pair<int, int> bounds{0, g.capacity()};
        bounds = initial_bounds(
            g, statistics, options.boundalg != gc::options::CLIQUES);
    } break;
    }
}
