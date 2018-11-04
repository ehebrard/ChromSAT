#include <iostream>

#include "brancher.hpp"
#include "dimacs.hpp"
#include "edgeformat.hpp"
#include "fillin.hpp"
#include "graph.hpp"
#include "mycielski.hpp"
#include "options.hpp"
#include "prop.hpp"
#include "reduction.hpp"
#include "rewriter.hpp"
#include "snap.hpp"
#include "sparse_dynamic_graph.hpp"
#include "statistics.hpp"
#include "utils.hpp"
#include "vcsolver.hpp"
#include "interval_list.hpp"
#include "dsatur.hpp"
#include "cliquesampler.hpp"

#include <minicsp/core/cons.cpp>
#include <minicsp/core/solver.hpp>
#include <minicsp/core/utils.hpp>

#include "./sota/Segundo/DSATUR/dsatur_algo.h"
#include "./sota/Segundo/DSATUR/graphe.h"

// #define DEBUG_DSIT

template <class graph_struct> void print(graph_struct& g)
{
    std::ofstream outfile("debug.col", std::ios_base::out);

    std::vector<int> vmap(g.capacity());
    graph_struct gc(g, vmap);

    outfile << "p edge " << gc.size() << " " << gc.count_edges() << std::endl;
    for (auto u : gc.nodes) {
        for (auto v : gc.matrix[u]) {
            if (u < v)
                outfile << "e " << (u + 1) << " " << (v + 1) << std::endl;
        }
    }
    outfile.close();
}

int mineq(const int N, const int k)
{
    int a = N / k;
    return a * (N - k * (a + 1) / 2);
}

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
struct gc_model {
    const gc::options& options;
    gc::statistics& statistics;

    bool degeneracy_sol;
    bool dsatur_sol;
    bool search_sol;

    int lb, ub;

    std::vector<int> debug_sol;

    // stores the coloring
    std::vector<int> solution;

    // maps vertices of the original graph to vertices of g
    std::vector<int> vertex_map;

    gc::graph<adjacency_struct>& original;
    gc::degeneracy_finder<gc::graph<adjacency_struct>> df;
    gc::clique_sampler cs;
    gc::dsatur col;
    adjacency_struct toremove;
    gc::graph_reduction<adjacency_struct> reduction;
    gc::dense_graph final;

    // gc::dyngraph dg;

    boost::optional<gc::fillin_info> fillin;

    minicsp::Solver solver;
    gc::varmap vars;
    gc::cons_base* cons{NULL};
    std::vector<minicsp::cspvar> xvars;
    std::unique_ptr<gc::rewriter> rewriter;

    std::unique_ptr<gc::Brancher> brancher;

    minicsp::cons_pb* mineq_constraint{NULL};

    gc::varmap create_chord_vars(minicsp::Solver& s, gc::dense_graph& g)
    {
        gc::minfill_buffer<gc::dense_graph> mb{g};
        mb.minfill();
        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[modeling] " << mb.fillin.size()
                      << " edges in minfill, width = " << mb.width << "\n";

        fillin = gc::fillin_info{std::move(mb.fillin), std::move(mb.order)};

        using std::begin;
        using std::end;
        gc::varmap vars(begin(g.nodes), end(g.nodes));
        vars.vars.resize(g.size());
        for (auto e : fillin->edges) {
            int i = e.first, j = e.second;
            Var v = s.newVar();
            vars.vars[i][j] = v;
            vars.vars[j][i] = v;
            if (options.trace) {
							
                using namespace std::string_literals;
                using std::to_string;
                auto n = "e"s + to_string(i) + "-"s + to_string(j);

                std::cout << std::setw(7) << n << " (" << i << "." << j << ")"
                          << std::endl;

                s.setVarName(vars.vars[i][j], n);
            }
        }

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[modeling] created " << s.nVars()
                      << " chord variables at " << minicsp::cpuTime() << "\n\n";

        return vars;
    }

    gc::varmap create_all_vars(minicsp::Solver& s, gc::dense_graph& g)
    {
        using std::begin;
        using std::end;

        std::vector<minicsp::Var> all_vars;

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
                if (options.equalities)
                    all_vars.push_back(vars.vars[i][j]);
                if (options.trace) {
                    using namespace std::string_literals;
                    using std::to_string;
                    auto n = "e"s + to_string(i) + "-"s + to_string(j);
                    s.setVarName(vars.vars[i][j], n);
                }
            }
        }

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[modeling] created " << s.nVars()
                      << " classic variables\n\n";

        if (options.equalities) {

            // std::cout << "UB = " << ub << std::endl;
            // int k = ub-1;
            // int a = final.capacity() / k;
            // int mineq = a * (final.capacity() - k * (a + 1) / 2);

            std::vector<int> w;
            w.resize(all_vars.size(), 1);
            auto m{mineq(final.capacity(), ub - 1)};
            mineq_constraint
                = new cons_pb(s, all_vars, w, m);

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[modeling] create implied constraint #eq >= " << m
                          << " / " << all_vars.size() << "\n\n";
        }

        return vars;
    }

    template <class map_struct>
    gc::varmap create_vars(
        minicsp::Solver& s, gc::dense_graph& g, map_struct& vmap)
    {

        if (g.size() > 0 and options.ddsaturiter > 0 and lb < ub) {

					std::cout << "HERE\n";

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[modeling] launch dense dsatur ("
                          << options.ddsaturiter << " times) at "
                          << minicsp::cpuTime() << "\n";

            // we have a dense graph now, so let's compute dsatur the other way
            // just in case
            for (int i = 0; i < options.ddsaturiter; ++i) {
                auto sol{gc::brelaz_color(g, (options.ddsaturiter > 1))};
                int ncol{*max_element(begin(sol), end(sol)) + 1};

                if (ub > ncol) {
                    assert(g.size() <= original.size());
                    for (int i = 0; i < g.size(); ++i) {
                        solution[vmap[i]] = sol[i];
                    }

                    auto actualncol = reduction.extend_solution(solution, ub, true);

                    if (ub > actualncol) {
												dsatur_sol = true;
                        ub = actualncol;
                        statistics.notify_ub(ub);
                        statistics.display(std::cout);
                    }
                }
            }

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[preprocessing] finished at "
                          << minicsp::cpuTime() << "\n";
        }

        if (options.fillin)
            return create_chord_vars(s, g);
        else
            return create_all_vars(s, g);
    }

    /*

    - compute k-cores, degeneracy and greedy ub

    - find a clique (lb) by sampling

    - remove all but the lb-core

    - run sparse dsatur


    */

    void probe_lb(const int samplebase_)
    {
        int samplebase = std::min(samplebase_, original.size());
        // assert(begin(df.order) + samplebase <= end(df.order));

        cs.set_domain((end(df.order) - samplebase), end(df.order),
            original.capacity(), true);

        size_t width{1};
        while (width <= options.probewidth) {
            auto nlb{cs.find_clique(
                original, lb, end(df.order), end(df.order), samplebase, width)};
            if (nlb > lb) {
                lb = nlb;

                statistics.notify_lb(lb);
                statistics.display(std::cout);
            } else {
                width *= 2;
            }
        }
    }

    void reduce(
        gc::graph_reduction<adjacency_struct>& gr, const int threshold, int& k)
    {

        if (df.core_degree[k] >= threshold)
            return;

        auto prev{k};
        while (df.core_degree[k] < threshold) {
            ++k;
        }

        for (auto vp{df.core[prev]}; vp != df.core[k]; ++vp) {
            auto v{*vp};

            toremove.add(v);
            gr.removed_vertices.push_back(v);
        }

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[preprocessing] remove " << toremove.size()
                      << " low degree nodes (<" << threshold << ")\n";

        toremove.canonize();

        original.remove(toremove);
        for (auto u : toremove) {
            gr.status[u] = gc::vertex_status::low_degree_removed;
        }

        toremove.clear();
    }

    gc::dense_graph dsatur_reduced()
    {
        col.get_core(original, options.core, lb, ub);

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[search] extract dsatur-core " << col.core.size()
                      << " / " << (col.frontier - begin(col.order) + 1) << " / "
                      << original.size() << ", current ub = " << ub
                      << std::endl;

        size_t stamp{0};
        original.nodes.save(stamp);
        original.nodes.clear();

        gc::bitset utilset(original.nodeset);

        original.nodeset.clear();

        // for (auto vp{begin(col.core)}; vp != end(col.core); ++vp) {
        for (auto vp{begin(col.core)}; vp < end(col.core); ++vp) {
            auto v{*vp};
            original.add_node(v);
        }

        gc::graph<gc::bitset> g(original, vertex_map);

        original.nodes.restore(stamp);
        // for (auto v : original.nodes)
        //     original.nodeset.fast_add(v);
        original.nodeset.copy(utilset);

        return g;
    }

    gc::graph<gc::vertices_vec> dsatur_sparse_reduced()
    {
        col.get_core(original, options.core, lb, ub);

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[search] extract dsatur-core " << col.core.size()
                      << " / " << (col.frontier - begin(col.order) + 1) << " / "
                      << original.size() << ", current ub = " << ub
                      << std::endl;

        size_t stamp{0};
        original.nodes.save(stamp);
        original.nodes.clear();

        gc::bitset utilset(original.nodeset);

        original.nodeset.clear();

        // for (auto vp{begin(col.core)}; vp != end(col.core); ++vp) {
        for (auto vp{begin(col.core)}; vp < end(col.core); ++vp) {
            auto v{*vp};
            original.add_node(v);
        }

        gc::graph<gc::vertices_vec> g(original, vertex_map);

        original.nodes.restore(stamp);
        // for (auto v : original.nodes)
        //     original.nodeset.fast_add(v);
        original.nodeset.copy(utilset);

        return g;
    }

    void upper_bound()
    {
        if (options.sdsaturiter > 0 and !dsatur_sol and lb < ub) {

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[preprocessing] launch sparse dsatur ("
                          << options.sdsaturiter << " times) at "
                          << minicsp::cpuTime() << "\n";

            col.clear();
            int niter{options.sdsaturiter};
            do {
                if (options.verbosity >= gc::options::YACKING) {
                    std::cout << ".";
                    std::cout.flush();
                }
                auto ncol{col.brelaz_color(original, ub - 1,
                    (1 << (options.sdsaturiter + 1 - niter)), 12345 + niter)};

                //                 std::cout << " ==> " << ncol << std::endl;
                //                 col.get_core(original, ub - 1, true);
                //                 col.get_core(original, ub - 1, false);
                // col.select();

                if (ncol < ub) {
                    for (int i = 0; i < original.size(); ++i) {
                        auto v{original.nodes[i]};
                        solution[v] = col.color[v];
                    }
                    dsatur_sol = true;

                    auto actualncol = reduction.extend_solution(solution, true);

                    if (ub > actualncol) {
                        ub = actualncol;
                        statistics.notify_ub(ub);

                        if (options.verbosity >= gc::options::NORMAL) {
                            if (options.verbosity >= gc::options::YACKING)
                                std::cout << std::endl;
                            statistics.display(std::cout);
                        }

                        // auto n_rel{1 + col.frontier - begin(col.order)};

                        // std::cout << "relevant vertices: " << n_rel
                        //           << std::endl;
                    }
                }
				if (--niter > 0) // Do not clear last col
                	col.clear();
            } while (niter > 0);
        }
        if (options.verbosity >= gc::options::YACKING)
            std::cout << std::endl;
    }

    bool peeling(
        gc::graph_reduction<adjacency_struct>& gr, const int k_core_threshold)
    {
        int original_size{original.size()};
        // std::cout << "[preprocessing] start peeling (" << k_core_threshold
        //           << ")\n";

        auto threshold = std::max(lb, k_core_threshold);

        // bool ub_safe{(threshold <= lb)};
        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[preprocessing] compute degeneracy ordering\n";
        df.degeneracy_ordering();

        for (int i{1}; i < df.core.size(); ++i) {

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[info] " << df.core_degree[i - 1]
                          << "-core: " << (df.core[i] - df.core[i - 1])
                          << " vertices [" << *(df.core[i - 1]) << "--"
                          << *(df.core[i] - 1) << "]\n";
        }

        // std::cout << ub << " <= " << (df.degeneracy + 1) << std::endl;

        if (ub > df.degeneracy + 1) {

            int maxc{gr.greedy_solution(
                solution, rbegin(df.order), rend(df.order), df.degeneracy + 1)};
            degeneracy_sol = true;

            assert(maxc <= df.degeneracy + 1);
            // assert(maxc >= lb);

            ub = maxc;
            // if (ub_safe) {

            std::cout << "DEGENERACY UB = " << ub << std::endl;
            statistics.notify_ub(ub);

            if (options.verbosity >= gc::options::NORMAL)
                statistics.display(std::cout);
            // }
        }

        // adjacency_struct toremove;
        // toremove.initialise(0, g.capacity(), gc::bitset::empt);

        auto samplebase{(original.size() > options.samplebase
                ? (options.samplebase
                      + (original.size() - options.samplebase) / 100)
                : original.size())};

        int prev_size{original_size + 1};
        int k{0};

        while (original.size() and prev_size > original.size() and lb < ub) {

            prev_size = original.size();

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[preprocessing] compute lower bound\n";
            probe_lb(samplebase);

            threshold = std::max(lb, threshold);
            reduce(gr, threshold, k);

            if (prev_size > original.size() or prev_size == original_size) {
                if (options.preprocessing == gc::options::FULL)
                    neighborhood_dominance(gr);
            }
        }

        // for(auto vi{df.core[k]}; vi!=end(df.order); ++vi) {
        // 	std::cout << (end(df.order) - vi) << ": " << df.degrees[*vi] <<
        // std::endl;
        // }

        // minimal_core(g, gr, ub, 200);

        return (original_size > original.size());
    }

    // reduce the graph so that it has either no more than 'size_threshold'
    // vertices left, or a minimal 'k_core_threshold'-core
    void minimal_core(gc::graph_reduction<adjacency_struct>& gr,
        const int k_core_threshold, const int size_threshold)
    {
        std::cout << "[search] search for minimal core (n=" << size_threshold
                  << ", k=" << k_core_threshold << ")\n";

        auto vi{end(df.order)};
        bool is_k_core{false};
        for (; vi-- > begin(df.order);) {
            // std::cout << (end(df.order) - vi) << ": " << *vi << " " <<
            // df.degrees[*vi] << " "
            //           << is_k_core << std::endl;
            if (is_k_core) {
                if ((end(df.order) - vi) > size_threshold
                    // or df.degrees[*vi] < lb
                    ) {
                    // std::cout << "KEEP!" << std::endl;
                    ++vi;
                    break;
                }
            } else {
                is_k_core = (df.degrees[*vi] >= k_core_threshold);
            }
        }

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[search] computed minimal core (n="
                      << (end(df.order) - vi) << ", k=" << df.degrees[*vi]
                      << ")\n";

        // std::cout << "keep";
        // for (auto tk{vi}; tk != end(df.order); ++tk) {
        //     std::cout << " " << *tk;
        // }
        // std::cout << std::endl;

        assert(toremove.size() == 0);

        for (auto vp{begin(df.order)}; vp != vi; ++vp) {
            auto v{*vp};

            if (original.nodeset.fast_contain(v)) {
                toremove.add(v);
                gr.removed_vertices.push_back(v);
            }
        }

        toremove.canonize();

        original.remove(toremove);
        for (auto u : toremove) {
            gr.status[u] = gc::vertex_status::low_degree_removed;
        }

        toremove.clear();

        if (options.verbosity >= gc::options::YACKING) {
            std::cout << "[search] minimal core: ";
            original.describe(std::cout, -1);
            std::cout << std::endl << std::endl;
        }
    }


    bool neighborhood_dominance(gc::graph_reduction<adjacency_struct>& gr)
    {

        int size_before{original.size()};

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[preprocessing] neighborhood dominance\n";

        std::vector<int> nodes;
        for (auto u : original.nodes)
            nodes.push_back(u);

        gc::bitset removed(0, original.capacity() - 1, gc::bitset::empt);
        gc::bitset utilset(0, original.capacity() - 1, gc::bitset::empt);

        int limit = 100000;
        auto psize{original.size()};
        for (auto u : nodes) {
            for (auto v : nodes)
                // check if u dominates v
                if (u != v and !original.matrix[u].fast_contain(v)
                    and original.nodeset.fast_contain(v)
                    and original.nodeset.fast_contain(u)) {
                    if (--limit <= 0)
                        break;

                    // if (limit % 1000 == 0)
                    //     std::cout << limit << std::endl;

                    utilset.copy(original.matrix[v]);
                    utilset.setminus_with(original.matrix[u]);
                    if (!utilset.intersect(original.nodeset)) {
                        // N(v) <= N(U)s

                        assert(!removed.fast_contain(v));

                        ++limit;
                        removed.add(v);
                        //
                        // std::cout << "\nrm " << v << " " <<
                        // original.matrix[v]
                        // <<
                        // std::endl;
                        // std::cout << "bc " << u << " " << original.matrix[u]
                        // <<
                        // std::endl;

                        original.remove(v);
                        gr.removed_vertices.push_back(v);
                        gr.dominator.push_back(u);
                        gr.status[v] = gc::vertex_status::dominated_removed;
                    }
                }
            if (limit <= 0)
                break;
        }

        if (psize > original.size()) {
            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[preprocessing] remove "
                          << (psize - original.size()) << " dominated nodes\n";
            if (options.verbosity >= gc::options::NORMAL)
                statistics.display(std::cout);
        }

        return (size_before > original.size());
    }

    // template< class adjacency_struct >
    void find_is_constraints(gc::graph_reduction<adjacency_struct>& gr)
    {
        gc::degeneracy_vc_solver<gc::graph<adjacency_struct>> vc(original);

        auto bs = vc.find_is();

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[preprocessing] extract IS constraint size = "
                      << bs.size() << "\n";

        for (auto v : bs) {
            gr.removed_vertices.push_back(v);
            gr.status[v] = gc::vertex_status::indset_removed;
            gr.constraints.emplace_back(
                gc::indset_constraint{original.matrix[v], v});
            original.nodes.remove(v);
            original.nodeset.remove(v);
        }
        for (auto v : original.nodes)
            original.matrix[v].intersect_with(original.nodeset);

        statistics.notify_removals(original.size());

        if (options.verbosity >= gc::options::NORMAL)
            statistics.display(std::cout);
    }

    gc::graph_reduction<adjacency_struct> preprocess(const int k_core_threshold)
    {
        // auto gr{degeneracy_peeling(original)};

        gc::graph_reduction<adjacency_struct> gr(
            original, statistics, solution);
        if (options.preprocessing == gc::options::NO_PREPROCESSING)
            return gr;
        else
            peeling(gr, k_core_threshold);

        if (options.verbosity >= gc::options::YACKING) {
            std::cout << "[preprocessing] preprocessed graph: ";
            original.describe(std::cout, -1);
            std::cout << std::endl << std::endl;
        }

        if (original.size() > 0 and lb < ub) {
            upper_bound();

            if (options.indset_constraints
                and options.strategy != gc::options::BOUNDS) {
                find_is_constraints(gr);
            }
        }

        if (options.verbosity >= gc::options::YACKING) {
            std::cout << "[preprocessing] finished at " << minicsp::cpuTime()
                      << "\n\n[modeling] preprocessed graph: ";
            original.describe(std::cout, -1);
            std::cout << std::endl << std::endl;
        }

        return gr;
    }

    void post_eqvar_debug_sol(Solver& s, const std::vector<int>& coloring)
    {
        std::vector<int> vertex_revmap(final.capacity());
        for (size_t i = 0; i != vertex_map.size(); ++i)
            if (vertex_map[i] >= 0)
                vertex_revmap[vertex_map[i]] = i;
        std::vector<int> sol(s.nVars());
        for (int i = 0; i != final.capacity(); ++i)
            for (int j = i + 1; j < final.capacity(); ++j)
                if (vars[i][j] != var_Undef) {
                    int ci = coloring[vertex_revmap[i]];
                    int cj = coloring[vertex_revmap[j]];
                    if (ci == cj)
                        sol[vars[i][j]] = 1;
                    else
                        sol[vars[i][j]] = 0;
                }
        s.debug_solution = sol;
    }

    template <class map_struct>
    void init_search(minicsp::Solver& s, gc::dense_graph& g, map_struct& vmap)
    {
        if (!options.dsatur) {
            vars = gc::varmap(create_vars(s, g, vmap));						
            cons = gc::post_gc_constraint(s, g, fillin, vars,
                reduction.constraints, vertex_map, options, statistics);

            double vm_usage;
            double resident_set;
            gc::process_mem_usage(vm_usage, resident_set);

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[modeling] created coloring constraint #nodes = "
                          << g.size() << ", #vars = " << s.nVars()
                          << ", memory = " << (long)resident_set << " \n";
        }

        if (!debug_sol.empty())
            post_eqvar_debug_sol(solver, debug_sol);

        // g.tell_class();
        // cons->g.tell_class();
        // cons->cf.g.tell_class();
        //

        rewriter = std::make_unique<gc::rewriter>(s, g, cons, vars, xvars);

        setup_signal_handlers(&s);
        s.trace = options.trace;
        s.polarity_mode = options.polarity;
        s.verbosity = (options.verbosity >= gc::options::SOLVERINFO);

        if (options.learning == gc::options::NO_LEARNING)
            s.learning = false;

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
                    int btlvl) { return rewriter->rewrite(clause, btlvl); });
            }

            auto sum = [](int x, int y) { return x + y; };
            auto prod = [](int x, int y) { return x * y; };

            switch (options.branching) {
            case gc::options::VSIDS:
                if (options.branching_low_degree) {
                    brancher = std::make_unique<gc::VSIDSBrancher>(
                        s, g, cons->fg, vars, xvars, *cons, options);
                    brancher->use();
                } else
                    s.varbranch = minicsp::VAR_VSIDS;
                break;
            case gc::options::VSIDS_GUIDED:
                if (options.branching_low_degree) {
                    brancher = std::make_unique<gc::VSIDSBrancher>(
                        s, g, cons->fg, vars, xvars, *cons, options);
                    brancher->use();
                } else
                    s.varbranch = minicsp::VAR_VSIDS;
                s.phase_saving = false;
                s.solution_phase_saving = true;
                break;
            case gc::options::VSIDS_PHASED:
                brancher = std::make_unique<gc::VSIDSPhaseBrancher>(
                    s, g, cons->fg, vars, xvars, *cons, options, -1, -1);
                brancher->use();
                break;
            case gc::options::VSIDS_CLIQUE:
                brancher = std::make_unique<gc::VSIDSCliqueBrancher>(
                    s, g, cons->fg, vars, xvars, *cons, options);
                break;
            case gc::options::VSIDS_COLORS_POSITIVE:
                if (!options.xvars) {
                    std::cout << "VSIDS_COLORS_"
                                 "POSITIVE needs "
                                 "--xvars\n";
                    exit(1);
                }
                brancher = std::make_unique<gc::VSIDSColorBrancher>(
                    s, g, cons->fg, vars, xvars, *cons, options);
                brancher->use();
                break;
            case gc::options::BRELAZ:
                brancher = std::make_unique<gc::BrelazBrancher>(
                    s, g, cons->fg, vars, xvars, *cons, options);
                brancher->use();
                break;
            case gc::options::PARTITION_SUM:
                brancher = gc::make_partition_brancher<-1, -1>(
                    s, g, cons->fg, vars, xvars, *cons, options, sum);
                brancher->use();
                break;
            case gc::options::PARTITION_PRODUCT:
                brancher = gc::make_partition_brancher<-1, -1>(
                    s, g, cons->fg, vars, xvars, *cons, options, prod);
                brancher->use();
                break;
            case gc::options::DEGREE_SUM:
                brancher = gc::make_degree_brancher<-1, -1>(
                    s, g, cons->fg, vars, xvars, *cons, options, sum);
                brancher->use();
                break;
            case gc::options::DEGREE_PRODUCT:
                brancher = gc::make_degree_brancher<-1, -1>(
                    s, g, cons->fg, vars, xvars, *cons, options, prod);
                brancher->use();
                break;
            case gc::options::DEGREE_UNION:
                brancher = std::make_unique<gc::DegreeUnionBrancher<-1, -1>>(
                    s, g, cons->fg, vars, xvars, *cons, options);
                brancher->use();
                break;
            case gc::options::PARTITION_SUM_DYN:
                brancher = gc::make_partition_brancher<2, 3>(
                    s, g, cons->fg, vars, xvars, *cons, options, sum);
                brancher->use();
                break;
            case gc::options::PARTITION_PRODUCT_DYN:
                brancher = gc::make_partition_brancher<2, 3>(
                    s, g, cons->fg, vars, xvars, *cons, options, prod);
                brancher->use();
                break;
            case gc::options::DEGREE_SUM_DYN:
                brancher = gc::make_degree_brancher<2, 3>(
                    s, g, cons->fg, vars, xvars, *cons, options, sum);
                brancher->use();
                break;
            case gc::options::DEGREE_PRODUCT_DYN:
                brancher = gc::make_degree_brancher<2, 3>(
                    s, g, cons->fg, vars, xvars, *cons, options, prod);
                brancher->use();
                break;
            case gc::options::DEGREE_UNION_DYN:
                brancher = std::make_unique<gc::DegreeUnionBrancher<2, 3>>(
                    s, g, cons->fg, vars, xvars, *cons, options);
                brancher->use();
                break;
            }
        }
  }

    gc_model(gc::graph<adjacency_struct>& ig, const gc::options& options,
        gc::statistics& statistics, std::pair<int, int> bounds,
        const std::vector<int>& debug_sol, const int k_core_threshold = -1)
        : options(options)
        , statistics(statistics)
        , degeneracy_sol{false}
        , dsatur_sol{false}
        , search_sol{false}
        , lb{bounds.first}
        , ub{bounds.second}
        , debug_sol(debug_sol)
        , solution(ig.capacity())
        , vertex_map(ig.capacity(), -1)
        , original(ig)
        , df(original)
        , toremove(0, ig.capacity(), 0)
        , reduction{preprocess(k_core_threshold)}
    {
        if (options.strategy != gc::options::BOUNDS and original.size() > 0
            and lb < ub) {

            std::cout << "create dense graph\n";

            final = gc::dense_graph(original, vertex_map);

            init_search(solver, final, original.nodes);
        }
    }

    template <class graph_struct>
    void get_solution(graph_struct& g, std::vector<int>& col)
    {
        int next{0};
        for (auto u : g.nodes) {
            for (auto v : g.partition[u]) {
                col[original.nodes[v]] = next;
            }
            ++next;
        }
    }

    // same as above expect that this is not a clique. "order" is the degeneracy
    // ordering
    int get_chordal_solution(std::vector<int>& col, std::vector<int>& order)
    {
        int maxc{0};
        std::vector<int> nodeset; //(0, g.capacity(), bitset::empt);
        gc::bitset colors(0, final.capacity(), gc::bitset::empt);

        for (auto i = order.rbegin(), iend = order.rend(); i != iend; ++i) {
            auto v = *i;
            colors.fill();
            for (auto u : nodeset) {
                if (!final.matrix[v].fast_contain(u))
                    continue;
                colors.fast_remove(col[original.nodes[u]]);
            }

            auto q{colors.min()};
            maxc = std::max(maxc, q);
            for (auto u : final.partition[v])
                col[original.nodes[u]] = q;
            nodeset.push_back(v);
        }
				
				return maxc + 1;
    }

    template <class graph_struct> int dsat_extend(graph_struct& g)
    {
        col.clear();

        auto actualub{ub};
        auto ncol{col.brelaz_color_guided(g, ub - 1, begin(g.nodes),
            begin(g.nodes) + col.core.size(), solution, 100, 12345 + ub)};

        if (ncol < ub) {
            actualub = reduction.extend_solution(solution, ub, true);
        } else
            actualub = ub;

        return actualub;
    }

    bool find_solution() { return find_solution(solver, final, lb, ub); }

    // color the residual graph
    //
    // stops at the first improving solution, return ub when there is none
    template <class graph_struct>
    bool find_solution(
        minicsp::Solver& s, graph_struct& g, int& lower_bound, int& upper_bound)
    {
        using minicsp::l_False;
        using minicsp::l_True;
        using minicsp::l_Undef;

        // int r{ub};

        minicsp::lbool sat{s.solveBudget()};

        assert(cons->ub >= upper_bound);
        cons->ub = upper_bound;

        if (sat == l_True) {
            upper_bound = g.nodes.size();
            get_solution(g, solution);

            // assert(solub < cons->ub);
            cons->sync_graph();

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[trace] SAT: " << cons->bestlb << ".."
                          << upper_bound;

            auto actualub{dsat_extend(original)};

            if (options.verbosity >= gc::options::YACKING)
                std::cout << ".." << actualub << std::endl;

            if (actualub < ub) {
                ub = actualub;
                statistics.notify_ub(ub);
            }

            return true;
        }
        lower_bound = upper_bound;
        return false;
    }

    // color the residual graph
    // model:lb and model:ub store the lower and upper bounds of the
    // RESIDUAL
    // graph
    //
    // store the solution in model::solution and extends it w.r.t. the IS
    // only
    bool solve(int ub_limit = -1)
    {
        return solve(solver, final, lb, ub, ub_limit);
    }

    template <class graph_struct>
    bool solve(minicsp::Solver& s, graph_struct& g, int& lower_bound,
        int& upper_bound, const int ub_limit = -1,
        const bool standard_extend = true)
    {
        using minicsp::l_False;
        using minicsp::l_True;
        using minicsp::l_Undef;

        minicsp::lbool sat{l_True};

        auto solution_found{false};
        while (sat != l_False and lower_bound < upper_bound) {

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[trace] solve in [" << cons->bestlb << ".."
                          << cons->ub << "[\n";

            sat = s.solveBudget();
						
            if (sat == l_True) {
                solution_found = true;
                search_sol = true;

                int solub = g.nodes.size();
								
                if (options.fillin) {

                    gc::degeneracy_finder<gc::dense_graph> df_leaf{g};
                    df_leaf.degeneracy_ordering();
                    int degeneracy{0};
                    for (auto v : df_leaf.order)
                        if (df_leaf.degrees[v] > degeneracy)
                            degeneracy = df_leaf.degrees[v];

                    solub = (degeneracy + 1);

                    assert(solub < cons->ub);

										int ncol{get_chordal_solution(solution, df_leaf.order)};
										
										assert(ncol == solub);

                } else {

                    get_solution(g, solution);
                }

                assert(solub < cons->ub);
                cons->sync_graph();

                if (options.verbosity >= gc::options::YACKING)
                    std::cout << "[trace] SAT: " << cons->bestlb << ".."
                              << solub;

                int actualub{upper_bound};
                int ISub{solub};
                if (standard_extend) {
                    // extends to the IS
                    ISub = reduction.extend_solution(
                        solution, upper_bound, false);
                    assert(ISub
                        <= upper_bound + (options.indset_constraints ? 1 : 0));

                    actualub = reduction.extend_solution(
                        solution, upper_bound, true);

                    if (options.verbosity >= gc::options::YACKING)
                        std::cout << ".." << ISub << ".." << actualub
                                  << std::endl;

                } else {

                    // std::cout << "SOL EXTEND UB = " << solub << std::endl;

                    col.clear();
                    auto ncol{col.brelaz_color_guided(original, upper_bound - 1,
                        begin(original.nodes),
                        begin(original.nodes) + col.core.size(), solution, 100,
                        12345 + upper_bound)};

                    if (options.verbosity >= gc::options::YACKING)
                        std::cout << ".." << ncol;

                    if (ncol < ub) {

                        actualub = reduction.extend_solution(
                            solution, upper_bound, true);

                        if (options.verbosity >= gc::options::YACKING)
                            std::cout << ".." << actualub << std::endl;

                        assert(ncol == actualub);

                        ub = actualub;
                    } else {
                        actualub = ncol;

                        if (options.verbosity >= gc::options::YACKING)
                            std::cout << ".." << actualub << std::endl;
                    }
                }

                // std::cout << "SOLVE UB = " << actualub << std::endl;
                statistics.notify_ub(actualub);
                if (options.verbosity >= gc::options::NORMAL)
                    statistics.display(std::cout);
                cons->ub = upper_bound = ISub;

                if (options.verbosity >= gc::options::YACKING)
                    std::cout << "bounds = [" << lower_bound << ".."
                              << upper_bound << "]" << std::endl;

                if (options.equalities) {
                    auto m{mineq(g.capacity(), cons->ub - 1)};
                    // mineq_constraint->set_lb(m);

                    if (options.verbosity >= gc::options::YACKING)
                        std::cout
                            << "[modeling] update implied constraint #eq >= "
                            << m << " / "
                            << (g.capacity() * (g.capacity() - 1) / 2)
                            << "\n\n";
                }

                if (options.xvars) {
                    for (auto v : xvars)
                        v.setmax(s, cons->ub - 2, minicsp::NO_REASON);
                }


                if (actualub < ub_limit)
                    break;

            } else if (sat == l_Undef) {
                if (options.verbosity >= gc::options::YACKING)
                    std::cout << "[trace] *** INTERRUPTED ***\n";
                break;
            } else {
                cons->bestlb = cons->ub;

                if (options.verbosity >= gc::options::YACKING)
                    std::cout << "[trace] UNSAT: " << cons->bestlb << ".."
                              << cons->ub << std::endl;

                statistics.notify_lb(cons->bestlb);

                if (options.verbosity >= gc::options::NORMAL)
                    statistics.display(std::cout);
                lower_bound = cons->bestlb;
            }
        }

        return solution_found;
    }

    void solve_with_dsatur()
    {
        C_Graphe G;
        DSATUR_ dsat_;

        G.nb_sommets = final.size();
        G.nb_aretes = original.count_edges();

        G.matrice_adjacence.resize(G.nb_sommets);
        G.sommets_voisins.resize(G.nb_sommets);

        G.adj = (bool*)malloc(G.nb_sommets * G.nb_sommets * sizeof(bool));
        G.sommets_voisins_bis
            = (int*)malloc(G.nb_sommets * G.nb_sommets * sizeof(int));
        G.degre = (int*)malloc(G.nb_sommets * sizeof(int));

        for (int i = 0; i < G.nb_sommets; i++) {
            G.matrice_adjacence[i].resize(G.nb_sommets);
            G.sommets_voisins[i].clear();
            G.degre[i] = 0;
        }
        for (int i = 0; i < G.nb_sommets; i++) {
            for (int j = 0; j < G.nb_sommets; j++) {
                G.matrice_adjacence[i][j] = 0;
                G.sommets_voisins_bis[i * G.nb_sommets + j] = 0;
                G.adj[i * G.nb_sommets + j] = false;
            }
        }

        for (auto u : final.nodes)
            for (auto v : final.matrix[u])
                if (G.matrice_adjacence[u][v] == 0) {
                    G.matrice_adjacence[u][v] = 1;
                    G.matrice_adjacence[v][u] = 1;

                    G.adj[u * G.nb_sommets + v] = true;
                    G.adj[v * G.nb_sommets + u] = true;

                    G.sommets_voisins[u].push_back(v);
                    G.sommets_voisins[v].push_back(u);

                    G.sommets_voisins_bis[(u * G.nb_sommets) + G.degre[u]] = v;
                    G.degre[u]++;
                    G.sommets_voisins_bis[(v * G.nb_sommets) + G.degre[v]] = u;
                    G.degre[v]++;
                }

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[trace] use external dsatur on [" << lb << ".." << ub
                      << "[\n";

        // std::vector<int> tmp_solution(G.nb_sommets);
        dsat_.print_progress = false;
        dsat_.DSATUR_algo(G, 10000, 2, lb, ub);

        if (ub > dsat_.UB) {
            search_sol = true;
            assert(original.size() == dsat_.init_n);
            assert(G.nb_sommets == dsat_.init_n);
            for (int v = 0; v < original.size(); ++v) {
                solution[original.nodes[v]] = dsat_.meilleure_coloration[v];
            }

            auto maxcol{*std::max_element(dsat_.meilleure_coloration,
                dsat_.meilleure_coloration + G.nb_sommets)};

            assert(maxcol == dsat_.UB - 1);
            for (auto a : original.nodes) {
                for (auto b : original.matrix[a]) {
                    assert(solution[a] != solution[b]);
                }
            }
        }

        lb = dsat_.LB;
        ub = dsat_.UB;

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[trace] update bounds [" << lb << ".." << ub << "]\n";
        }

        void finalize_solution(std::vector<pair<int, int>>& edges)
        {
            reduction.extend_solution(solution, ub, true);
            auto ncol{*std::max_element(begin(solution), end(solution)) + 1};
            std::cout << "[solution] " << ncol << "-coloring computed at "
                      << minicsp::cpuTime() << std::endl
                      << std::endl;

            if (options.printsolution) {
                for (int v = 0; v < original.capacity(); ++v)
                    std::cout << " " << std::setw(2) << solution[v];
                std::cout << std::endl;
            }

            if (options.checksolution) {
                for (auto e : edges) {
                    if (solution[e.first] == solution[e.second]) {
                        std::cout << "WRONG SOLUTION: " << e.first << " and "
                                  << e.second << " <- " << solution[e.first]
                                  << "\n";
                        break;
                    }
                }
            }
        }

        void print_stats()
        {
            if (statistics.best_lb >= statistics.best_ub)
                std::cout << "OPTIMUM " << statistics.best_ub << "\n";
            else
                std::cout << "Best bounds [" << statistics.best_lb << ", "
                          << statistics.best_ub << "]\n";

            minicsp::printStats(solver);
            statistics.display(std::cout);
            std::cout << std::endl;
        }
};


template <class input_format>
int color(gc::options& options, gc::graph<input_format>& g)
{
    options.describe(std::cout);

    if (options.verbosity >= gc::options::NORMAL)
        std::cout << "[reading] ";

    std::vector<std::pair<int, int>> edges;
    if (options.format == "snap")
        snap::read_graph(options.instance_file.c_str(),
            [&](int nv, int) { g = gc::graph<input_format>{nv}; },
            [&](int u, int v) {
                if (u != v) {
                    // num_edges += 1 - g.matrix[u].fast_contain(v);
                    g.add_edge(u, v);
                    // ++num_edges;
                }
            },
            [&](int, gc::weight) {});
    else if (options.format == "edg")
        edgeformat::read_graph(options.instance_file.c_str(),
            [&](int nv, int ne) {
                g = gc::graph<input_format>{nv};
                // num_edges = ne;
            },
            [&](int u, int v) { g.add_edge(u, v); }, [&](int, gc::weight) {});
    else
        dimacs::read_graph(options.instance_file.c_str(),
            [&](int nv, int) { g = gc::graph<input_format>{nv}; },
            [&](int u, int v) {
                if (u != v) {
                    g.add_edge(u - 1, v - 1);
                    // ++num_edges;
                    if (options.checksolution)
                        edges.push_back(std::pair<int, int>{u - 1, v - 1});
                }
            },
            [&](int, gc::weight) {});

    g.canonize();
    long num_edges{g.count_edges()};


    if (options.convert != "") {
        std::ofstream outfile(options.convert.c_str(), std::ios_base::out);
        if (options.convert.substr(options.convert.size() - 3, 3) == "edg") {
            outfile << g.size() << " " << num_edges << "\n";
            for (auto u : g.nodes) {
                for (auto v : g.matrix[u]) {
                    if (u < v)
                        outfile << u << " " << v << "\n";
                }
            }
        } else {
            gc::dyngraph dg(g);
            dg.print_dimacs(outfile);
        }
        return 1;
    }

    std::vector<int> sol;
    if (!options.solution_file.empty()) {
        std::cout << "[reading] Reading solution file\n";
        std::ifstream sifs(options.solution_file.c_str());
        int c{0};
        while (sifs) {
            sifs >> c;
            if (sifs)
                sol.push_back(c);
        }
    }

    if (options.verbosity >= gc::options::NORMAL) {
        g.describe(std::cout, num_edges);
        std::cout << " at " << minicsp::cpuTime() << std::endl;
    }

    gc::statistics statistics(g.capacity());
    if (options.preprocessing != gc::options::NO_PREPROCESSING)
        statistics.update_ub = false;

    if (options.verbosity >= gc::options::NORMAL) {
        statistics.describe(std::cout);
        std::cout << std::endl;
    }

    // std::cout << "MAIN (READ): " << g.size() << "(" << (int*)(&g) << ")"
    //           << std::endl;
    switch (options.strategy) {
    case gc::options::COLOR: {

        std::cout << "NOT IMPLEMENTED!\n";

    } break;
    case gc::options::BNB: {
        std::pair<int, int> bounds{0, g.size()};
        gc_model<input_format> model(g, options, statistics, bounds, sol);

        model.final.describe(std::cout);
        std::cout << std::endl;

        // model.solve(model.solver, model.final);
        model.solve();

        model.finalize_solution(edges);

        if (options.verbosity >= gc::options::QUIET)
            model.print_stats();
    } break;
    case gc::options::IDSATUR: {
			

        std::pair<int, int> bounds{0, g.size()};

        options.strategy = gc::options::BOUNDS; // so that we don't create the
        // dense graph yet
        options.ddsaturiter = 0;
        gc_model<input_format> model(g, options, statistics, bounds, sol);

        int limit{-1};

        while (model.lb < model.ub) {
            statistics.binds(NULL);
            gc::dense_graph g{model.dsatur_reduced()};

            minicsp::Solver s;

            model.init_search(s, g, model.original.nodes);

            statistics.binds(model.cons);
						

						int nub{model.ub}, nlb{model.lb};
            auto solution_found{false};
            if (options.idsaturlimit == 0) {
                solution_found = model.find_solution(s, g, nlb, nub);
            } else {
                solution_found = model.solve(s, g, nlb, nub, -1, false);
            }
						

            model.lb = std::max(nlb, model.lb);
            statistics.notify_lb(model.lb);

            // // we need to check lb < ub because the solution is already
            // extended in
            //             if (solution_found and model.lb < model.ub) {
            //                 // auto solub = model.cons->ub;
            //
            //                 if (options.verbosity >= gc::options::YACKING)
            //                     std::cout << "[trace] " << limit
            //                               << " SAT: " << model.cons->bestlb
            //                               << ".." << nub;
            //
            // 		std::cout << "DSAT EXTEND UB = " << nub << std::endl;
            //                 auto actualub{model.dsat_extend(model.original)};
            //
            //                 if (options.verbosity >= gc::options::YACKING)
            //                     std::cout << ".." << actualub << std::endl;
            //
            //                 if (actualub < model.ub) {
            //                     model.ub = actualub;
            //                     statistics.notify_ub(model.ub);
            //                 }
            //             }

            if (options.verbosity >= gc::options::NORMAL)
                statistics.display(std::cout);

            if (g.size() == model.original.size())
                break;

            if (--limit == 0)
                break;
        }

    } break;
    case gc::options::BOTTOMUP: {

        // gc::statistics& sub_statistics(g.capacity());

        std::pair<int, int> bounds{0, g.size()};

        options.strategy = gc::options::BOUNDS; // so that we don't create the
        // dense graph yet
        options.ddsaturiter = 0;
        gc_model<input_format> model(g, options, statistics, bounds, sol);

        // model.original.describe(std::cout);
        // std::cout << std::endl;

        int limit{-1};

        options.strategy = gc::options::BNB;
        while (model.lb < model.ub) {

            statistics.binds(NULL);
            // model.cons = NULL;

            // std::cout << "model.lb: " << model.lb << " / model.ub: " <<
            // model.ub << std::endl;
            gc::graph<input_format> g{model.dsatur_sparse_reduced()};
            // g.describe(std::cout);
            // std::cout << std::endl;

            // std::cout << "----\n";
            // std::cout << g.nodes << std::endl;
            // std::cout << g.nodeset << std::endl;
            // for( auto v : g.nodes ) {
            // 	std::cout << g.matrix[v] << std::endl;
            // }

            std::pair<int, int> bounds{model.lb, model.ub};

            statistics.update_ub = false;
            gc_model<input_format> tmp_model(
                g, options, statistics, bounds, sol);
            statistics.update_ub = true;

            // std::cout << "tmp_model.lb: " << tmp_model.lb << " /
            // tmp_model.ub: " << tmp_model.ub << std::endl;
            // std::cout << "model.lb: " << model.lb << " / model.ub: " <<
            // model.ub << std::endl;

            // minicsp::Solver s;
            //
            // model.init_search(s, g, model.original.nodes);
            //
            // statistics.binds(model.cons);

            auto solution_found{tmp_model.ub < model.ub};
            if (options.idsaturlimit == 0) {
                solution_found |= model.find_solution(tmp_model.solver,
                    tmp_model.final, tmp_model.lb, tmp_model.ub);
            } else {
                solution_found |= model.solve(tmp_model.solver, tmp_model.final,
                    tmp_model.lb, tmp_model.ub, -1, false);
                //     solution_found = model.find_solution();
                // } else {
                //     solution_found = model.solve(-1);
            }

            model.lb = std::max(tmp_model.lb, model.lb);
            statistics.notify_lb(model.lb);

            // we need to check lb < ub because the solution is already extended
            // in
            if (solution_found and model.lb < model.ub) {
                // auto solub = model.cons->ub;

                // std::cout << 11 << std::endl;

                if (options.verbosity >= gc::options::YACKING)
                    std::cout << "[trace] " << limit << " SAT: " << tmp_model.lb
                              << ".." << tmp_model.ub;

                // std::cout << 22 << std::endl;

                // std::cout << "DSAT EXTEND UB = " << model.ub << std::endl;
                auto actualub{model.dsat_extend(model.original)};

                if (options.verbosity >= gc::options::YACKING)
                    std::cout << ".." << actualub << std::endl;

                if (actualub < model.ub) {
                    model.ub = actualub;
                    statistics.notify_ub(model.ub);
                }
            }

            // exit(1);

            if (options.verbosity >= gc::options::NORMAL)
                statistics.display(std::cout);

            if (g.size() == model.original.size())
                break;

            if (--limit == 0)
                break;
        }

    } break;
    case gc::options::TOPDOWN: {

        std::pair<int, int> bounds{1, g.size()};

        options.strategy = gc::options::BOUNDS; // so that we don't create the
        // dense graph yet
        options.ddsaturiter = 0;
        gc_model<input_format> model(g, options, statistics, bounds, sol);

        // std::cout << model.lb << ".." << model.ub << std::endl;

        model.col.local_search(model.original, model.solution, statistics);





    } break;
    case gc::options::CLEVER: {

        std::vector<int> vmap(g.capacity(), -1);
        options.strategy = gc::options::BOUNDS; // so that we don't create the
        // dense graph yet
        gc_model<gc::vertices_vec> init_model(
            g, options, statistics, std::make_pair(0, g.size()), sol);

        options.strategy = gc::options::CLEVER;
        statistics.update_lb = false;

        // std::vector<std::pair<int, int>> tmp_edges;
        // for (auto u : init_model.g.nodes)
        //     for (auto v : init_model.g.matrix[u])
        //         if (u < v)
        //             tmp_edges.push_back(std::pair<int, int>{u, v});

        while (init_model.lb < init_model.ub) {

            if (options.verbosity >= gc::options::YACKING)
                std::cout << "[search] solve a tmp model with bounds ["
                          << init_model.lb << ".." << init_model.ub
                          << "[ focusing on the " << (init_model.ub - 1)
                          << "-core\n";

            vmap.clear();
            vmap.resize(g.capacity(), -1);

            gc::graph<gc::vertices_vec> gcopy(g, vmap);

            // print(gcopy);

            gc_model<gc::vertices_vec> tmp_model(gcopy, options, statistics,
                std::make_pair(init_model.lb, init_model.ub), sol,
                (init_model.ub - 1));
								

								// std::cout << "[" << tmp_model.lb << ".." << tmp_model.ub << "]\n";
								// assert(tmp_model.g.size() > 0 or tmp_model.ub > tmp_model.lb);

            if (tmp_model.ub > tmp_model.lb and tmp_model.final.size() > 0) {
                if (options.dsatur) {
                    tmp_model.solve_with_dsatur();
                } else {
                    tmp_model.solve(init_model.ub);
                    // tmp_model.solve(
                    //     tmp_model.solver, tmp_model.final, init_model.ub);
                }
            }

            if (options.verbosity >= gc::options::YACKING) {
                std::cout << "[search] tmp: [" << tmp_model.lb << "..";
                if (tmp_model.ub < init_model.ub)
                    std::cout << tmp_model.ub << "..";
            }

            assert(tmp_model.solution.size() == tmp_model.original.capacity());

            auto incumbent{init_model.ub};
            if (tmp_model.degeneracy_sol or tmp_model.dsatur_sol
                or tmp_model.search_sol) {
                incumbent = tmp_model.reduction.extend_solution(
                    tmp_model.solution, init_model.ub, true);
                // copy the tmp model solution into the init model
                for (int v = 0; v < tmp_model.original.capacity(); ++v)
                    init_model.solution[init_model.original.nodes[v]]
                        = tmp_model.solution[v];
		            assert(incumbent < init_model.ub);
		            init_model.ub = incumbent;
            }

            if (options.verbosity >= gc::options::YACKING)
                std::cout << init_model.ub << "]\n";

            init_model.lb = tmp_model.lb;

            // iub may not be equal to ilb even if the solver wasn't stopped:
            // coloring the removed vertices might have required extra colors
            // either way, [ilb, iub] are correct bounds
            statistics.notify_ub(init_model.ub);
            statistics.notify_lb(init_model.lb);

            if (options.verbosity >= gc::options::NORMAL)
                statistics.display(std::cout);

            statistics.unbinds();
            vmap.clear();
        }

        init_model.finalize_solution(edges);

    } break;
    case gc::options::BOUNDS: {
        std::pair<int, int> bounds{0, g.size()};
        gc_model<input_format> model(g, options, statistics, bounds, sol);
        // int is_ub{0};
        model.reduction.extend_solution(model.solution, bounds.second, true);
        auto ncol{
            *std::max_element(begin(model.solution), end(model.solution)) + 1};

        if (options.verbosity >= gc::options::YACKING)
            std::cout << "[solution] " << ncol << "-coloring computed at "
                      << minicsp::cpuTime() << std::endl
                      << std::endl;
        if (options.checksolution) {
            for (auto e : edges) {
                if (model.solution[e.first] == model.solution[e.second]) {
                    std::cout << "WRONG SOLUTION!!\n";
                }
            }
        }
    } break;

    case gc::options::TEST: {

        std::cout << "NOT IMPLEMENTED!\n";

    } break;
    }

    return 0;
}

int main(int argc, char* argv[])
{
    auto options = gc::parse(argc, argv);
    gc::graph<gc::vertices_vec> g;
    // gc::graph<gc::bitset> g;
    auto result = color(options, g);

    return result;
}
