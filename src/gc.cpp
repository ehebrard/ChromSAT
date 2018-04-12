#include <iostream>

#include "brancher.hpp"
#include "dimacs.hpp"
#include "graph.hpp"
#include "options.hpp"
#include "prop.hpp"
#include "rewriter.hpp"
#include "statistics.hpp"
#include "utils.hpp"
#include <minicsp/core/cons.hpp>
#include <minicsp/core/solver.hpp>
#include <minicsp/core/utils.hpp>

struct graph_reduction {
    const gc::graph& g;
    std::vector<int> removed_vertices;
    gc::bitset nodeset;
    gc::bitset util_set;

    explicit graph_reduction(const gc::graph& g)
        : g(g)
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
            util_set.clear();
            for (auto u : g.matrix[v]) {
                if (!nodeset.fast_contain(u))
                    continue;
                util_set.fast_add(col[u]);
            }
            for (int q = 0; q != g.capacity(); ++q) {
                if (util_set.fast_contain(q))
                    continue;
                maxc = std::max(maxc, q);
                col[v] = q;
                break;
            }
            nodeset.fast_add(v);
        }
        return maxc + 1;
    }
};

struct gc_model {
    int lb{0}, ub{-1};
    const gc::options& options;
    gc::statistics& statistics;

    graph_reduction reduction;
    gc::graph& g;
    minicsp::Solver s;
    std::vector<std::vector<minicsp::Var>> vars;
    gc::cons_base* cons;
    std::vector<minicsp::cspvar> xvars;
    gc::rewriter rewriter;

    std::unique_ptr<gc::Brancher> brancher;

    // std::vector<edge> var_map;

    std::vector<std::vector<minicsp::Var>> create_vars()
    {
        std::vector<std::vector<minicsp::Var>> vars;
        vars.resize(g.capacity());
        for (size_t i = 0; i != vars.size(); ++i) {
            if (!g.nodes.contain(i))
                continue;
            vars[i].resize(g.capacity());
            vars[i][i] = minicsp::var_Undef;
            for (size_t j = 0; j != i; ++j) {
                if (!g.nodes.contain(j))
                    continue;
                vars[i][j] = vars[j][i];
            }
            for (size_t j = i + 1; j != vars.size(); ++j) {
                if (!g.nodes.contain(j))
                    continue;
                if (g.matrix[i].fast_contain(j))
                    vars[i][j] = minicsp::var_Undef;
                else {
                    vars[i][j] = s.newVar();
                    if (options.trace) {
                        using namespace std::string_literals;
                        using std::to_string;
                        auto n = "e"s + to_string(i) + "-"s + to_string(j);
                        s.setVarName(vars[i][j], n);
                    }
                }
            }
        }
        return vars;
    }

    graph_reduction preprocess(gc::graph& g, std::pair<int, int> bounds)
    {
        graph_reduction gr(g);
        if (options.preprocessing == gc::options::NO_PREPROCESSING)
            return gr;

        lb = bounds.first;
        ub = bounds.second;
        int hlb{0};
        gc::clique_finder cf{g};
        gc::bitset forbidden(0, g.capacity(), gc::bitset::empt);
        gc::bitset util_set(0, g.capacity(), gc::bitset::empt);
        gc::bitset removedv(0, g.capacity(), gc::bitset::empt);
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
            if (ub < 0 || (hub < ub && hub >= lb)) {
                ub = hub;
                statistics.notify_ub(ub);
            }

            hlb = cf.find_cliques(g.nodes);
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
                ++statistics.num_vertex_removals;
                toremove.push_back(u);
                gr.removed_vertices.push_back(u);
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
        }
        return gr;
    }

    gc_model(gc::graph& g, const gc::options& options,
        gc::statistics& statistics, std::pair<int, int> bounds)
        : options(options)
        , statistics(statistics)
        , reduction(preprocess(g, bounds))
        , g(g)
        , vars(create_vars())
        , cons(gc::post_gc_constraint(s, g, vars, options, statistics))
        , rewriter(s, g, *cons, vars, xvars)
    {
        setup_signal_handlers(&s);
        s.trace = options.trace;
        s.polarity_mode = options.polarity;

        if (options.learning == gc::options::NO_LEARNING)
            s.learning = false;

        auto [plb, pub] = bounds;
        lb = std::max(lb, plb);
        if (ub < 0 || pub < ub)
            ub = pub;

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
                        minicsp::post_eq_re(
                            s, xvars[i], xvars[j], 0, minicsp::Lit(vars[i][j]));
                }
            }

            // rewrite clauses to not use x literals
            s.use_clause_callback([this](vec<minicsp::Lit>& clause, int btlvl) {
                return rewriter.rewrite(clause, btlvl);
            });
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
        case gc::options::BRELAZ:
            if (!options.xvars) {
                std::cout << "Cannot use Brelaz ordering without xvars\n";
                exit(1);
            }
            brancher = std::make_unique<gc::BrelazBrancher>(
                s, g, vars, xvars, *cons, options);
            brancher->use();
            break;
						//         case gc::options::BRELAZ_GUIDED:
						//             if (!options.xvars) {
						//                 std::cout << "Cannot use Brelaz ordering without xvars\n";
						//                 exit(1);
						//             }
						//             brancher = std::make_unique<gc::BrelazBrancher>(
						//                 s, g, vars, xvars, *cons, options);
						//             brancher->use();
						// s.phase_saving = false;
						// s.solution_phase_saving = true;
						//             break;
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

    std::vector<int> get_solution()
    {
        std::vector<int> col(g.capacity());
        int next{0};
        for (auto u : g.nodes) {
            for (auto v : g.partition[u])
                col[v] = next;
            ++next;
        }
        return col;
    }

    std::pair<int, int> solve()
    {
        using minicsp::l_False;
        using minicsp::l_True;
        using minicsp::l_Undef;

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
                if (actualub != solub)
                    std::cout << " UB in reduced graph = " << solub << std::endl;
                cons->ub = solub;
                if (options.xvars) {
                    for (auto v : xvars)
                        v.setmax(s, cons->ub - 2, minicsp::NO_REASON);
                }
            } else if (sat == l_Undef) {
                std::cout << "*** INTERRUPTED ***\n";
                break;
            } else {
                cons->bestlb = cons->ub;
                statistics.display(std::cout);
            }
        }
        return std::make_pair(cons->bestlb, cons->ub);
    }

    void print_stats()
    {
        if (cons->bestlb >= cons->ub)
            std::cout << "OPTIMUM " << cons->ub << "\n";
        else
            std::cout << "Best bounds [" << cons->bestlb << ", " << cons->ub
                      << "]\n";
        minicsp::printStats(s);
        statistics.display(std::cout);
        std::cout << std::endl;
    }
};

std::pair<int, int> initial_bounds(const gc::graph& g, gc::statistics& stat, bool myciel=false)
{
    gc::clique_finder cf{g};
    gc::mycielskan_subgraph_finder mf(g, cf, false);
    int lb{cf.find_cliques(g.nodes)};

		if(myciel)
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

void histogram(gc::graph& g)
{
    std::vector<int> degrees;
    gc::bitset N(0, g.capacity()-1, gc::bitset::empt);
    for (auto v : g.nodes) {
        N.copy(g.matrix[v]);
        N.intersect_with(g.nodeset);
        degrees.push_back(N.size());
    }
    std::sort(begin(degrees), end(degrees));
    int psize = degrees.size()/10;
    int pstart{0};
    while (static_cast<size_t>(pstart) < degrees.size()) {
        int pend
            = std::min(static_cast<size_t>(pstart + psize), degrees.size() - 1);
        while (static_cast<size_t>(pend) < degrees.size() - 1
            && degrees[pend + 1] == degrees[pend])
            ++pend;
        std::cout << (pend - pstart + 1) << " vertices: degree "
                  << degrees[pstart] << "-" << degrees[pend] << "\n";
        pstart = pend+1;
    }
}

int main(int argc, char* argv[])
{
    auto options = gc::parse(argc, argv);
    options.describe(std::cout);

    gc::graph g;
    dimacs::read_graph(options.instance_file.c_str(),
        [&](int nv, int) { g = gc::graph{nv}; },
        [&](int u, int v) {
            if (u != v)
                g.add_edge(u - 1, v - 1);
        },
        [&](int, gc::weight) {});
    g.describe(std::cout);
    histogram(g);

    gc::statistics statistics(g.capacity());
    if (options.preprocessing)
        statistics.update_ub = false;
		
		statistics.describe(std::cout);
		

    switch (options.strategy) {
    case gc::options::BNB: {
        std::pair<int, int> bounds{0, g.capacity()};
        if (options.preprocessing == gc::options::NO_PREPROCESSING)
            bounds = initial_bounds(g, statistics, options.boundalg!=gc::options::CLIQUES);
        gc_model model(g, options, statistics, bounds);
        model.solve();
        model.print_stats();
        break;
    }
    case gc::options::BOTTOMUP: {
        statistics.update_ub = false;
        auto [lb, ub] = initial_bounds(g, statistics, options.boundalg!=gc::options::CLIQUES);
        for (int i = lb; i < ub; ++i) {
            gc::graph gcopy{g};
            gc_model model(
                gcopy, options, statistics, std::make_pair(i, i + 1));
            auto [ilb, iub] = model.solve();
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
        auto [lb, ub] = initial_bounds(g, statistics, options.boundalg!=gc::options::CLIQUES);
        for (int i = ub - 1; i >= lb; --i) {
            gc::graph gcopy{g};
            gc_model model(
                gcopy, options, statistics, std::make_pair(i, i + 1));
            auto [ilb, iub] = model.solve();
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
				bounds = initial_bounds(g, statistics, options.boundalg!=gc::options::CLIQUES);
    } break;
    }

}
