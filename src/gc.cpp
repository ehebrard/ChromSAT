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

#include <minicsp/core/cons.cpp>
#include <minicsp/core/solver.hpp>
#include <minicsp/core/utils.hpp>

#include "./sota/Segundo/DSATUR/dsatur_algo.h"
#include "./sota/Segundo/DSATUR/graphe.h"

template <class graph_struct> void print(graph_struct& g)
{
    std::ofstream outfile("debug.col", std::ios_base::out);

    outfile << "p edge " << g.size() << " " << g.count_edges() << std::endl;
    for (auto u : g.nodes) {
        for (auto v : g.matrix[u]) {
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

// template< class adjacency_struct >
// struct graph_reduction {
//     const gc::graph<adjacency_struct>& g;
//     const gc::statistics& statistics;
//     std::vector<int> removed_vertices;
//     std::vector<int> dominator;
//     std::vector<gc::indset_constraint> constraints;
//     std::vector<vertex_status> status;
//     gc::bitset nodeset;
//     gc::bitset util_set;
//
//     explicit graph_reduction(
//         const gc::graph<adjacency_struct>& g, const gc::statistics&
//         statistics)
//         : g(g)
//         , statistics(statistics)
//         , status(g.capacity(), gc::vertex_status::in_graph)
//         , nodeset(0, g.capacity(), gc::bitset::empt)
//         , util_set(0, g.capacity(), gc::bitset::empt)
//     {
//     }
//
//     int extend_solution(std::vector<int>& col, const bool full = false)
//     {
//         int maxc{0};
//         nodeset.copy(g.nodeset);
//         for (auto v : g.nodes)
//             maxc = std::max(maxc, col[v]);
//
//         auto d{dominator.rbegin()};
//         for (auto i = removed_vertices.rbegin(), iend =
//         removed_vertices.rend();
//              i != iend; ++i) {
//             auto v = *i;
//
//             if (!full and status[v] != gc::vertex_status::indset_removed)
//                 break;
//
//             if (status[v] == gc::vertex_status::dominated_removed) {
//                 assert(nodeset.fast_contain(*d));
//                 col[v] = col[*d];
//                 ++d;
//             }
//
//             assert(status[v] != gc::vertex_status::in_graph);
//             util_set.clear();
//             for (auto u : g.matrix[v]) {
//                 if (!nodeset.fast_contain(u))
//                     continue;
//                 util_set.fast_add(col[u]);
//             }
//
//             for (int q = 0; q != g.capacity(); ++q) {
//                 if (util_set.fast_contain(q))
//                     continue;
//                 // assert(q <= statistics.best_ub);
//                 maxc = std::max(maxc, q);
//                 col[v] = q;
//                 break;
//             }
//             nodeset.fast_add(v);
//         }
//
//         return maxc + 1;
//     }
// };

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
    gc::graph_reduction<adjacency_struct> reduction;
    gc::dense_graph g;

    // gc::dyngraph dg;

    boost::optional<gc::fillin_info> fillin;

    minicsp::Solver s;
    gc::varmap vars;
    gc::cons_base* cons{NULL};
    std::vector<minicsp::cspvar> xvars;
    std::unique_ptr<gc::rewriter> rewriter;

    std::unique_ptr<gc::Brancher> brancher;

    minicsp::cons_pb* mineq_constraint{NULL};

    // std::vector<edge> var_map;

    gc::varmap create_chord_vars()
    {
        gc::minfill_buffer<gc::dense_graph> mb{g};
        mb.minfill();
        std::cout << "[modeling] "<< mb.fillin.size()
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

        std::cout << "[modeling] created " << s.nVars()
                  << " chord variables at " << minicsp::cpuTime() << "\n\n";

        return vars;
    }

    gc::varmap create_all_vars()
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

        std::cout << "[modeling] created " << s.nVars()
                  << " classic variables\n\n";

        if (options.equalities) {

            // std::cout << "UB = " << ub << std::endl;
            // int k = ub-1;
            // int a = g.capacity() / k;
            // int mineq = a * (g.capacity() - k * (a + 1) / 2);

            std::vector<int> w;
            w.resize(all_vars.size(), 1);
            auto m{mineq(g.capacity(), ub - 1)};
            mineq_constraint
                = new cons_pb(s, all_vars, w, m);
            std::cout << "[modeling] create implied constraint #eq >= " << m
                      << " / " << all_vars.size() << "\n\n";
        }

        return vars;
    }

    gc::varmap create_vars()
    {

        if (g.size() > 0 and options.ddsaturiter > 0 and lb < ub) {

            std::cout << "[modeling] launch dense dsatur ("
                      << options.ddsaturiter << " times) at "
                      << minicsp::cpuTime() << "\n";

            // we have a dense graph now, so let's compute dsatur the other way
            // just in case
            for (int i = 0; i < options.ddsaturiter; ++i) {
                auto sol{gc::brelaz_color(g, (options.ddsaturiter > 1))};
                int ncol{*max_element(begin(sol), end(sol)) + 1};

                if (ub > ncol) {
                    assert(g.size() == original.size());
                    for (int i = 0; i < original.size(); ++i) {
                        solution[original.nodes[i]] = sol[i];
                    }

                    auto actualncol = reduction.extend_solution(solution, true);

                    if (ub > actualncol) {
												dsatur_sol = true;
                        ub = actualncol;
                        statistics.notify_ub(ub);
                        statistics.display(std::cout);
                    }
                }
            }

            std::cout << "[preprocessing] finished at " << minicsp::cpuTime()
                      << "\n";
        }

        if (options.fillin)
            return create_chord_vars();
        else
            return create_all_vars();
    }

    void degeneracy_peeling(gc::graph<adjacency_struct>& g,
        gc::graph_reduction<adjacency_struct>& gr, const int k_core_threshold)
    {
        // graph_reduction<adjacency_struct> gr(g, statistics);
        // if (options.preprocessing == gc::options::NO_PREPROCESSING)
        //     return gr;

        std::cout << "[preprocessing] start peeling (" << k_core_threshold << ")\n";

        // histogram(g);

        // lb = bounds.first;
        // ub = bounds.second;

        // auto lb = std::max(lb, given_lb);

        bool lb_safe = true;
        auto threshold = lb;
        gc::clique_finder<adjacency_struct> cf{
            g, std::min(options.cliquelimit, g.size())};
        gc::mycielskan_subgraph_finder<adjacency_struct> mf(g, cf, false);
        gc::degeneracy_finder<gc::graph<adjacency_struct>> df{g};

        adjacency_struct toremove;
        toremove.initialise(0, g.capacity(), gc::bitset::empt);

        int stop = 10;
        do {
            if (g.size() == 0) {
                break;
            }

            cf.clear();
            df.clear();
            // mf.clear();

            std::cout << "[preprocessing] compute degeneracy ordering\n";
            int degeneracy = 0;
            df.degeneracy_ordering();

            std::vector<int>
                reverse; // TODO change that by using iterators in clique finder
            for (auto rit = df.order.rbegin(); rit != df.order.rend(); ++rit) {
                auto v{*rit};
                if (df.degrees[v] > degeneracy) {
                    degeneracy = df.degrees[v];
                }
                reverse.push_back(v);
                // gr.status[v] = gc::vertex_status::low_degree_removed;
            }

            if (ub > degeneracy + 1) {
                ub = degeneracy + 1;
                statistics.notify_ub(ub);
                statistics.display(std::cout);

                if (!degeneracy_sol) {
                    for (auto v : df.order) {
                        gr.removed_vertices.push_back(v);
                        gr.status[v] = gc::vertex_status::low_degree_removed;
                    }
                    gr.extend_solution(solution, true);
                    gr.removed_vertices.clear();
                    gr.status.resize(g.capacity(), gc::vertex_status::in_graph);
                    degeneracy_sol = true;
                }
            }

            std::cout << "[preprocessing] compute lower bound\n";

            auto plb = cf.find_cliques(reverse);

            if (options.boundalg != gc::options::CLIQUES) {
                cf.sort_cliques(plb);
                plb = mf.improve_cliques_larger_than(plb);
            }

            bool changes = false;
            if (lb < plb) {

                if (lb_safe) {
                    lb = plb;
                    statistics.notify_lb(lb);
                    statistics.display(std::cout);
                    changes = true;
                }
            }

            if (k_core_threshold > threshold) {
                threshold = k_core_threshold;
                lb_safe = false;
            }

            if (lb > threshold) {
                threshold = lb;
            }

            if (df.degrees[df.order[0]] < threshold) {

                for (auto v : df.order) {
                    if (df.degrees[v] >= threshold)
                        break;
                    toremove.add(v);
                    gr.removed_vertices.push_back(v);
                }

                std::cout << "[preprocessing] remove " << toremove.size()
                          << " low degree nodes (<" << threshold << ")\n";

                toremove.canonize();

                g.remove(toremove);
                for (auto u : toremove) {
                    gr.status[u] = gc::vertex_status::low_degree_removed;
                }

                toremove.clear();
                statistics.notify_removals(g.size());
                statistics.display(std::cout);
                changes = true;
            }
            // else {
            //
            //                 std::cout << "lowest degree node " << df.order[0]
            //                 << " ("
            //                           << df.degrees[df.order[0]] << ") is >=
            //                           " << threshold
            //                           << std::endl;
            //             }

            if (!changes)
                break;

        } while (stop-- > 0);

        // return gr;
    }

    void neighborhood_dominance(gc::graph<adjacency_struct>& g,
        gc::graph_reduction<adjacency_struct>& gr)
    {

        std::vector<int> nodes;
        for (auto u : g.nodes)
            nodes.push_back(u);

        gc::bitset removed(0, g.capacity() - 1, gc::bitset::empt);

        auto psize{g.size()};
        for (auto u : nodes)
            for (auto v : nodes)
                // check if u dominates v
                if (u != v and !g.matrix[u].fast_contain(v)
                    and g.nodeset.fast_contain(v)
                    and g.nodeset.fast_contain(u)) {
                    gr.util_set.copy(g.matrix[v]);
                    gr.util_set.setminus_with(g.matrix[u]);
                    if (!gr.util_set.intersect(g.nodeset)) {
                        // N(v) <= N(U)s

                        assert(!removed.fast_contain(v));

                        removed.add(v);
                        //
                        // std::cout << "\nrm " << v << " " << g.matrix[v] <<
                        // std::endl;
                        // std::cout << "bc " << u << " " << g.matrix[u] <<
                        // std::endl;

                        g.remove(v);
                        gr.removed_vertices.push_back(v);
                        gr.dominator.push_back(u);
                        gr.status[v] = gc::vertex_status::dominated_removed;
                    }
                }

        if (psize > g.size())
            std::cout << "[preprocessing] remove " << (psize - g.size())
                      << " dominated nodes\n";
    }

    // gc::graph_reduction<adjacency_struct> core_reduction(
    //     gc::graph<adjacency_struct>& g, std::pair<int, int> bounds,
    //     bool myciel = false)
    // {
    //     graph_reduction<adjacency_struct> gr(g, statistics);
    //     if (options.preprocessing == gc::options::NO_PREPROCESSING)
    //         return gr;
    //
    //     // std::cout << "CORE REDUCTION: " << g.size() << "(" << (int*)(&g)
    //     << ")"
    //     //           << std::endl;
    //
    //     lb = bounds.first;
    //     ub = bounds.second;
    //     int hlb{0};
    //     gc::clique_finder<adjacency_struct> cf(g);
    //     gc::mycielskan_subgraph_finder<adjacency_struct> mf(g, cf, false);
    //     gc::degeneracy_finder<gc::graph<adjacency_struct>> df{g};
    //
    //     gc::bitset forbidden(0, g.capacity(), gc::bitset::empt);
    //     gc::bitset util_set(0, g.capacity(), gc::bitset::empt);
    //     adjacency_struct removedv(0, g.capacity(), gc::bitset::empt);
    //     std::vector<int> toremove;
    //     bool removed{false};
    //     int niteration{0};
    //     do {
    //         ++niteration;
    //         removed = false;
    //         auto sol{gc::brelaz_color(g)};
    //         for (auto u : g.nodes)
    //             for (auto v : g.matrix[u])
    //                 assert(sol[u] != sol[v]);
    //         int hub{*max_element(begin(sol), end(sol)) + 1};
    //
    //         // df.degeneracy_ordering();
    //         //             int hub{*max_element(
    //         //                 begin(df.degrees), end(df.degrees))};
    //         if (ub < 0 || (hub < ub && hub >= lb)) {
    //             ub = hub;
    //             statistics.notify_ub(ub);
    //         }
    //
    //         hlb = cf.find_cliques(g.nodes);
    //         if (myciel)
    //             hlb = mf.improve_cliques_larger_than(lb);
    //
    //         if (hlb > lb) {
    //             lb = hlb;
    //             statistics.notify_lb(lb);
    //         }
    //         statistics.display(std::cout);
    //
    //         forbidden.clear();
    //         toremove.clear();
    //         for (auto u : g.nodes) {
    //             if (forbidden.fast_contain(u))
    //                 continue;
    //             util_set.copy(g.matrix[u]);
    //             util_set.intersect_with(g.nodeset);
    //             if (util_set.size() >= static_cast<size_t>(lb))
    //                 continue;
    //             removed = true;
    //             removedv.fast_add(u);
    //             // ++statistics.num_vertex_removals;
    //             toremove.push_back(u);
    //             gr.removed_vertices.push_back(u);
    //             gr.status[u] = gc::vertex_status::low_degree_removed;
    //             forbidden.union_with(g.matrix[u]);
    //         }
    //         for (auto u : toremove) {
    //             g.nodes.remove(u);
    //             g.nodeset.remove(u);
    //         }
    //     } while (removed);
    //     if (removedv.size() > 0) {
    //         for (auto v : g.nodes) {
    //             g.matrix[v].setminus_with(removedv);
    //             g.origmatrix[v].setminus_with(removedv);
    //         }
    //         statistics.notify_removals(g.size());
    //         statistics.display(std::cout);
    //     }
    //
    //     return gr;
    // }

    // template< class adjacency_struct >
    void find_is_constraints(gc::graph<adjacency_struct>& g,
        gc::graph_reduction<adjacency_struct>& gr)
    {
        gc::degeneracy_vc_solver<gc::graph<adjacency_struct>> vc(g);
				
				//
				// std::cout << "PRINT GRAPH\n";
				// vector<int> vmap(g.capacity(), -1);
				// gc::graph<adjacency_struct> pcopy(g, vmap);
				// print(pcopy);
				
				
				
        auto bs = vc.find_is();
        std::cout << "[preprocessing] extract IS constraint size = "
                  << bs.size() << "\n";
        // std::cout << bs << std::endl;
        for (auto v : bs) {
            gr.removed_vertices.push_back(v);
            gr.status[v] = gc::vertex_status::indset_removed;
            gr.constraints.emplace_back(gc::indset_constraint{g.matrix[v], v});
            g.nodes.remove(v);
            g.nodeset.remove(v);
        }
        for (auto v : g.nodes)
            g.matrix[v].intersect_with(g.nodeset);

        statistics.notify_removals(g.size());
        statistics.display(std::cout);
    }

    void sparse_upper_bound(gc::dyngraph& dg, const int maxiter)
    {
        gc::coloring col;

        int iter = 0;
        do {

            // sparse upper bound, do it always
            col.brelaz_color(dg, (options.sdsaturiter > 1 ? 1 : 0));
            auto ncol{*std::max_element(begin(col.color), end(col.color)) + 1};

            // // [TODO: WHEN GIVING A "FAKE" LB, WE CAN FIND A SMALLER
            // COLORING,
            // // 1/ WE SHOULD CHECK WHAT IS THE ACTUAL SIZE OF THAT COLORING 2/
            // // OTHERWISE WE CAN ONLY GUARANTEE LB]
            // if (ncol < lb) {
            //     ncol = lb;
            // } [DONE]

            if (ub > ncol) {
                for (int i = 0; i < original.size(); ++i) {
                    solution[original.nodes[i]] = col.color[i];
                }
                dsatur_sol = true;
								
								std::cout << "DSATUR SOLUTION\n";

                auto actualncol = reduction.extend_solution(solution, true);
                // std::cout << " ====> " << actualncol << std::endl;

                if (ub > actualncol) {
                    ub = actualncol;
                    statistics.notify_ub(ub);
                    statistics.display(std::cout);
                }
            }

            if (++iter >= maxiter)
                break;

            dg.undo();
            col.clear();
        } while (lb < ub);
    }

    gc::graph_reduction<adjacency_struct> preprocess(
        gc::graph<adjacency_struct>& g, const int k_core_threshold)
    {
        // auto gr{degeneracy_peeling(original)};

        gc::graph_reduction<adjacency_struct> gr(g, statistics, solution);
        if (options.preprocessing == gc::options::NO_PREPROCESSING)
            return gr;

        degeneracy_peeling(original, gr, k_core_threshold);

        if (options.preprocessing == gc::options::FULL)
            neighborhood_dominance(original, gr);

        int num_edges = 0;
        if (g.size() > 0 and lb < ub) {

            if (options.sdsaturiter > 0) {
                std::cout << "[preprocessing] launch sparse dsatur ("
                          << options.sdsaturiter << " times) at "
                          << minicsp::cpuTime() << "\n";
                gc::dyngraph dg(g);
                num_edges = dg.edges.size();
                sparse_upper_bound(dg, options.sdsaturiter);
            }

            if (g.size() > 0 and lb < ub
                and options.indset_constraints
                and options.strategy != gc::options::BOUNDS)
                find_is_constraints(g, gr);
        }

        std::cout << "[preprocessing] finished at " << minicsp::cpuTime()
                  << "\n\n[modeling] preprocessed graph: ";
        g.describe(std::cout, num_edges);
        std::cout << std::endl << std::endl;

        return gr;
    }

    void post_eqvar_debug_sol(Solver& s, const std::vector<int>& coloring)
    {
        std::vector<int> vertex_revmap(g.capacity());
        for (size_t i = 0; i != vertex_map.size(); ++i)
            if (vertex_map[i] >= 0)
                vertex_revmap[vertex_map[i]] = i;
        std::vector<int> sol(s.nVars());
        for (int i = 0; i != g.capacity(); ++i)
            for (int j = i + 1; j < g.capacity(); ++j)
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
        , reduction{preprocess(original, k_core_threshold)}
    {

        if (options.strategy != gc::options::BOUNDS and original.size() > 0
            and lb < ub) {

            // std::cout << "[modeling] create dense graph with "<<
            // original.size()
            //           << " nodes\n";

            g = gc::dense_graph(original, vertex_map);

            // // g.check_consistency();
            // for(int i=0; i<g.size(); ++i) {
            //      std::cout << original.matrix[i] << " / " << g.matrix[i] <<
            //      std::endl;
            // }
            // std::cout << std::endl;
            // original.describe(std::cout);
            // std::cout << std::endl;
            // g.describe(std::cout);
            // std::cout << std::endl;

            if (!options.dsatur) {
                vars = gc::varmap(create_vars());
                cons = gc::post_gc_constraint(s, g, fillin, vars,
                    reduction.constraints, vertex_map, options, statistics);

                double vm_usage;
                double resident_set;
                gc::process_mem_usage(vm_usage, resident_set);

                std::cout << "[modeling] created coloring constraint #nodes = "
                          << original.size() << ", #vars = " << s.nVars()
                          << ", memory = " << (long)resident_set << " \n";
            }

            if (!debug_sol.empty())
                post_eqvar_debug_sol(s, debug_sol);


            // g.tell_class();
            // cons->g.tell_class();
            // cons->cf.g.tell_class();
            //

            rewriter = std::make_unique<gc::rewriter>(s, g, cons, vars, xvars);

            setup_signal_handlers(&s);
            s.trace = options.trace;
            s.polarity_mode = options.polarity;
            s.verbosity = options.verbosity;

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
                    s.use_clause_callback(
                        [this](vec<minicsp::Lit>& clause, int btlvl) {
                            return rewriter->rewrite(clause, btlvl);
                        });
                }

                auto sum = [](int x, int y) { return x + y; };
                auto prod = [](int x, int y) { return x * y; };

                switch (options.branching) {
                case gc::options::VSIDS:
                    if (options.branching_low_degree) {
                        brancher = std::make_unique<gc::VSIDSBrancher>(
                            s, g, cons->fg, vars,xvars, *cons, options);
                        brancher->use();
                    } else
                        s.varbranch = minicsp::VAR_VSIDS;
                    break;
                case gc::options::VSIDS_GUIDED:
                    if (options.branching_low_degree) {
                        brancher = std::make_unique<gc::VSIDSBrancher>(
                            s, g, cons->fg, vars,xvars, *cons, options);
                        brancher->use();
                    } else
                        s.varbranch = minicsp::VAR_VSIDS;
                    s.phase_saving = false;
                    s.solution_phase_saving = true;
                    break;
                case gc::options::VSIDS_PHASED:
                    brancher = std::make_unique<gc::VSIDSPhaseBrancher>(
                        s, g, cons->fg, vars,xvars, *cons, options, -1, -1);
                    brancher->use();
                    break;
                case gc::options::VSIDS_CLIQUE:
                    brancher = std::make_unique<gc::VSIDSCliqueBrancher>(
                        s, g, cons->fg, vars,xvars, *cons, options);
                    break;
                case gc::options::VSIDS_COLORS_POSITIVE:
                    if (!options.xvars) {
                        std::cout << "VSIDS_COLORS_"
                                     "POSITIVE needs "
                                     "--xvars\n";
                        exit(1);
                    }
                    brancher = std::make_unique<gc::VSIDSColorBrancher>(
                        s, g, cons->fg, vars,xvars, *cons, options);
                    brancher->use();
                    break;
                case gc::options::BRELAZ:
                    brancher = std::make_unique<gc::BrelazBrancher>(
                        s, g, cons->fg, vars,xvars, *cons, options);
                    brancher->use();
                    break;
                case gc::options::PARTITION_SUM:
                    brancher = gc::make_partition_brancher<-1, -1>(
                        s, g, cons->fg, vars,xvars, *cons, options, sum);
                    brancher->use();
                    break;
                case gc::options::PARTITION_PRODUCT:
                    brancher = gc::make_partition_brancher<-1, -1>(
                        s, g, cons->fg, vars,xvars, *cons, options, prod);
                    brancher->use();
                    break;
                case gc::options::DEGREE_SUM:
                    brancher = gc::make_degree_brancher<-1, -1>(
                        s, g, cons->fg, vars,xvars, *cons, options, sum);
                    brancher->use();
                    break;
                case gc::options::DEGREE_PRODUCT:
                    brancher = gc::make_degree_brancher<-1, -1>(
                        s, g, cons->fg, vars,xvars, *cons, options, prod);
                    brancher->use();
                    break;
                case gc::options::DEGREE_UNION:
                    brancher
                        = std::make_unique<gc::DegreeUnionBrancher<-1, -1>>(
                            s, g, cons->fg, vars,xvars, *cons, options);
                    brancher->use();
                    break;
                case gc::options::PARTITION_SUM_DYN:
                    brancher = gc::make_partition_brancher<2, 3>(
                        s, g, cons->fg, vars,xvars, *cons, options, sum);
                    brancher->use();
                    break;
                case gc::options::PARTITION_PRODUCT_DYN:
                    brancher = gc::make_partition_brancher<2, 3>(
                        s, g, cons->fg, vars,xvars, *cons, options, prod);
                    brancher->use();
                    break;
                case gc::options::DEGREE_SUM_DYN:
                    brancher = gc::make_degree_brancher<2, 3>(
                        s, g, cons->fg, vars,xvars, *cons, options, sum);
                    brancher->use();
                    break;
                case gc::options::DEGREE_PRODUCT_DYN:
                    brancher = gc::make_degree_brancher<2, 3>(
                        s, g, cons->fg, vars,xvars, *cons, options, prod);
                    brancher->use();
                    break;
                case gc::options::DEGREE_UNION_DYN:
                    brancher = std::make_unique<gc::DegreeUnionBrancher<2, 3>>(
                        s, g, cons->fg, vars,xvars, *cons, options);
                    brancher->use();
                    break;
                }
            }
        }
    }

    void get_solution(std::vector<int>& col)
    {
        // std::vector<int> col(original.capacity(), -1);
        int next{0};

        for (auto u : g.nodes) {
            for (auto v : g.partition[u]) {
                col[original.nodes[v]] = next;

#ifdef DEBUG_IS
                std::cout << original.nodes[v] << " <- " << next << std::endl;
// assert(v == original.nodes[v]);
// std::cout << v << " -> " << original.nodes[v] << std::endl;
#endif
            }
            ++next;
        }
        // return col;
    }

    // color the residual graph
    // model:lb and model:ub store the lower and upper bounds of the RESIDUAL
    // graph
    //
    // store the solution in model::solution and extends it w.r.t. the IS only
    void solve()
    {
        using minicsp::l_False;
        using minicsp::l_True;
        using minicsp::l_Undef;

        minicsp::lbool sat{l_True};
        while (sat != l_False && lb < ub) {

            std::cout << "[trace] solve* in [" << cons->bestlb << ".."
                      << cons->ub << "[\n";

            sat = s.solveBudget();
            if (sat == l_True) {
                search_sol = true;
								
                get_solution(solution);

								
                int solub = g.nodes.size();
								if (options.fillin) {
					        gc::degeneracy_finder<gc::graph<adjacency_struct>> df{g};
									df.degeneracy_ordering();
									int degeneracy{0};
									for (auto v : df.order) 
											if (df.degrees[v] > degeneracy) 
					            				degeneracy = df.degrees[v];
															assert (cons->bestlb == degeneracy + 1);
									
								}
								
                assert(solub < cons->ub);
                cons->sync_graph();

                // extends to the IS
                int ISub = reduction.extend_solution(solution, false);
                assert(ISub <= ub + 1);

                std::cout << "[trace] SAT: " << cons->bestlb << ".." << solub
                          << ".." << ISub ;

								int actualub = reduction.extend_solution(solution, true);
								
								std::cout << ".." << actualub << std::endl;
								
                statistics.notify_ub(actualub);
                statistics.display(std::cout);
                cons->ub = ub = ISub;
                if (options.equalities) {
                    auto m{mineq(g.capacity(), cons->ub - 1)};
                    // mineq_constraint->set_lb(m);
                    std::cout
                        << "[modeling] update implied constraint #eq >= " << m
                        << " / " << (g.capacity() * (g.capacity() - 1) / 2)
                        << "\n\n";
                }

                if (options.xvars) {
                    for (auto v : xvars)
                        v.setmax(s, cons->ub - 2, minicsp::NO_REASON);
                }

            } else if (sat == l_Undef) {
                std::cout << "[trace] *** INTERRUPTED ***\n";
                break;
            } else {
                cons->bestlb = cons->ub;

                std::cout << "[trace] UNSAT: " << cons->bestlb << ".."
                          << cons->ub << std::endl;

                statistics.notify_lb(cons->bestlb);
                statistics.display(std::cout);
                lb = cons->bestlb;
            }
        }
    }

    void solve_with_dsatur()
    {
        C_Graphe G;
        DSATUR_ dsat_;

        G.nb_sommets = g.size();
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

        for (auto u : g.nodes)
            for (auto v : g.matrix[u])
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

        std::cout << "[trace] update bounds [" << lb << ".." << ub << "]\n";
    }

    void finalize_solution(std::vector<pair<int, int>>& edges)
    {
        reduction.extend_solution(solution, true);
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

template <class input_format>
int color(gc::options& options, gc::graph<input_format>& g)
{
    options.describe(std::cout);

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

    // if (options.format != "edg")
    g.canonize();
    long num_edges{g.count_edges()};
    // for (auto v : g.nodes) {
    //     num_edges += g.matrix[v].size();
    // }
    // num_edges /= 2;

    if (options.convert != "") {
        std::ofstream outfile(options.convert.c_str(), std::ios_base::out);

        // std::cout << options.convert.substr(options.convert.size()-3, 3) <<
        // std::endl;
        // exit(1);

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

    g.describe(std::cout, num_edges);
    std::cout << " at " << minicsp::cpuTime() << std::endl;
    // histogram(g);

    gc::statistics statistics(g.capacity());
    if (options.preprocessing != gc::options::NO_PREPROCESSING)
        statistics.update_ub = false;

    statistics.describe(std::cout);
    std::cout << std::endl;

    // std::cout << "MAIN (READ): " << g.size() << "(" << (int*)(&g) << ")"
    //           << std::endl;
    switch (options.strategy) {
    case gc::options::BNB: {
        std::pair<int, int> bounds{0, g.size()};
        if (options.preprocessing == gc::options::NO_PREPROCESSING) {
            std::cout << "compute init bounds\n";
            bounds = initial_bounds(
                g, statistics, options.boundalg != gc::options::CLIQUES);
        }

        gc_model<input_format> model(g, options, statistics, bounds, sol);
        model.solve();

        model.finalize_solution(edges);

        model.print_stats();
    } break;
    case gc::options::BOTTOMUP: {

        std::cout << "NOT IMPLEMENTED!\n";

        // statistics.update_ub = false;
        // auto bounds = initial_bounds(
        //     g, statistics, options.boundalg != gc::options::CLIQUES);
        // auto lb = bounds.first;
        // auto ub = bounds.second;
        // for (int i = lb; i < ub; ++i) {
        //     gc::graph<gc::bitset> gcopy{g};
        //     gc_model<gc::bitset> model(
        //         gcopy, options, statistics, std::make_pair(i, i + 1));
        //     auto ibounds = model.solve();
        //     auto ilb = ibounds.first;
        //     auto iub = ibounds.second;
        //     if (iub == i) {
        //         statistics.notify_ub(iub);
        //         statistics.display(std::cout);
        //         break;
        //     } else if (ilb != iub) {
        //         statistics.notify_lb(ilb);
        //         statistics.display(std::cout);
        //         std::cout << "INTERRUPTED\n";
        //     } else {
        //         statistics.notify_lb(ilb);
        //         statistics.display(std::cout);
        //     }
        //     statistics.unbinds();
        // }
    } break;
    case gc::options::TOPDOWN: {

        std::cout << "NOT IMPLEMENTED!\n";

        // std::vector<int> vmap(g.capacity());
        // options.strategy = gc::options::BOUNDS; // so that we don't create
        // the
        //                                         // dense graph yet
        // gc_model<gc::vertices_vec> init_model(
        //     g, options, statistics, std::make_pair(0, g.size()));
        //
        // options.strategy = gc::options::TOPDOWN;
        // statistics.update_lb = false;
        // auto lb = init_model.lb;
        // auto ub = init_model.ub;
        //
        // for (int i = ub - 1; i >= lb; --i) {
        //
        //     std::cout << "[search] solve a tmp model with bounds [" << i <<
        //     ".."
        //               << (i + 1) << "]\n";
        //
        //     vmap.resize(g.capacity());
        //     gc::graph<gc::vertices_vec> gcopy(g, vmap);
        //
        //     gc_model<gc::vertices_vec> tmp_model(
        //         gcopy, options, statistics, std::make_pair(i, i + 1));
        //
        //     auto ibounds = tmp_model.solve();
        //
        //     auto ilb = ibounds.first;
        //     auto iub = ibounds.second;
        //
        //     if (ilb == i + 1) {
        //         statistics.notify_lb(ilb);
        //         statistics.display(std::cout);
        //         break;
        //     } else if (ilb != iub) {
        //         statistics.notify_ub(iub);
        //         statistics.display(std::cout);
        //         std::cout << "INTERRUPTED\n";
        //     } else {
        //         statistics.notify_ub(iub);
        //         statistics.display(std::cout);
        //     }
        //     statistics.unbinds();
        //     vmap.clear();
        // }
    } break;
    case gc::options::CLEVER: {

        std::vector<int> vmap(g.capacity(), -1);
        options.strategy = gc::options::BOUNDS; // so that we don't create the
        // dense graph yet
        gc_model<gc::vertices_vec> init_model(
            g, options, statistics, std::make_pair(0, g.size()), sol);

        options.strategy = gc::options::CLEVER;
        statistics.update_lb = false;

        std::vector<std::pair<int, int>> tmp_edges;
        for (auto u : init_model.g.nodes)
            for (auto v : init_model.g.matrix[u])
                if (u < v)
                    tmp_edges.push_back(std::pair<int, int>{u, v});

        int limit = 5;
        while (init_model.lb < init_model.ub) {

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

            if (tmp_model.ub > tmp_model.lb) {
                if (options.dsatur) {
                    tmp_model.solve_with_dsatur();
                } else {
                    tmp_model.solve();
                }
            }

            std::cout << "[search] tmp: [" << tmp_model.lb << "..";
            if (tmp_model.ub < init_model.ub)
                std::cout << tmp_model.ub << "..";

            assert(tmp_model.solution.size() == tmp_model.original.capacity());

            auto incumbent{init_model.ub};
            if (tmp_model.degeneracy_sol or tmp_model.dsatur_sol
                or tmp_model.search_sol) {
                incumbent = tmp_model.reduction.extend_solution(
                    tmp_model.solution, true);
                // copy the tmp model solution into the init model
                for (int v = 0; v < tmp_model.original.capacity(); ++v)
                    init_model.solution[init_model.original.nodes[v]]
                        = tmp_model.solution[v];
		            assert(incumbent < init_model.ub);
		            init_model.ub = incumbent;
            }

						std::cout << init_model.ub << "]\n";
            init_model.lb = tmp_model.lb;

            // iub may not be equal to ilb even if the solver wasn't stopped:
            // coloring the removed vertices might have required extra colors
            // either way, [ilb, iub] are correct bounds
            statistics.notify_ub(init_model.ub);
            statistics.notify_lb(init_model.lb);
            statistics.display(std::cout);

            statistics.unbinds();
            vmap.clear();

            if (--limit == 0)
                break;
        }

        init_model.finalize_solution(edges);

    } break;
    case gc::options::BOUNDS: {
        std::pair<int, int> bounds{0, g.size()};
        gc_model<input_format> model(g, options, statistics, bounds, sol);
        // int is_ub{0};
        model.reduction.extend_solution(model.solution, true);
        auto ncol{
            *std::max_element(begin(model.solution), end(model.solution)) + 1};
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
