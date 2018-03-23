#include "prop.hpp"
#include "minicsp/core/utils.hpp"
#include "utils.hpp"

namespace gc
{

using namespace minicsp;

class gc_constraint : public minicsp::cons, public cons_base
{
private:
    Solver& s;
    graph& g;
    const std::vector<std::vector<Var>>& vars;
    const options& opt;
	 	statistics& stat;

    struct varinfo_t {
        int u{-1}, v{-1};
    };
    // indexed by varid
    std::vector<varinfo_t> varinfo;

    //
    backtrackable<int> lastdlvl;

    // the usual assortment of vecs, vectors and bitsets to avoid
    // reallocations
    vec<Lit> reason;
    std::vector<int> culprit;
    bitset diffuv, diffvu;

    bitset expl_N, expl_covered, expl_residue;

    neighbors_wrapper adjacency_list;

    std::vector<int> degeneracy_order;
    std::vector<int> heuristic;

    // bitset cur_neighbors;
    // std::vector<int> maxcliques;
    // std::vector<std::vector<int>> n_included;
    // long int n_prunings;

public:
    gc_constraint(Solver& solver, graph& pg,
        const std::vector<std::vector<Var>>& tvars, const options& opt, statistics& stat)
        : cons_base(pg)
        , s(solver)
        , g(pg)
        , vars(tvars)
        , opt(opt)
				, stat(stat)
        , lastdlvl(s)
        , expl_N(0, g.capacity() - 1, bitset::empt)
        , expl_covered(0, g.capacity() - 1, bitset::empt)
        , expl_residue(0, g.capacity() - 1, bitset::empt)
        , adjacency_list(g)
    {
        ub = g.capacity();
        assert(vars.size() == static_cast<size_t>(g.capacity()));
        for (int i = 0; i != g.capacity(); ++i) {
            assert(vars[i].size() == static_cast<size_t>(g.capacity()));
            for (int j = 0; j != g.capacity(); ++j) {
                if (j < i) {
                    assert(vars[i][j] == vars[j][i]);
                    continue;
                }
                if (j == i)
                    continue;
                if (vars[i][j] == minicsp::var_Undef)
                    continue;
                if (varinfo.size() <= static_cast<size_t>(vars[i][j]))
                    varinfo.resize(vars[i][j] + 1);
                varinfo[vars[i][j]] = {i, j};
                s.wake_on_lit(vars[i][j], this, nullptr);
                s.schedule_on_lit(vars[i][j], this);
            }
        }
        diffuv.initialise(0, g.capacity(), bitset::empt);
        diffvu.initialise(0, g.capacity(), bitset::empt);
        // cur_neighbors.initialise(0, g.capacity(), bitset::empt);

        DO_OR_THROW(propagate(s));
    }

    std::ostream& print(Solver&, std::ostream& os) const
    {
        os << "coloring";
        return os;
    }

    std::ostream& printstate(Solver&, std::ostream& os) const
    {
        os << "coloring";
        return os;
    }

    void sync_graph()
    {
        if (*lastdlvl < g.current_checkpoint()) {
            g.restore(*lastdlvl);
            *lastdlvl = g.current_checkpoint();
        }
        while (s.decisionLevel() > g.current_checkpoint()) {
            g.checkpoint();
            *lastdlvl = g.current_checkpoint();
        }
    }

    Clause* wake(Solver& s, Lit l)
    {
        sync_graph();

        auto info = varinfo[var(l)];
        auto u{g.rep_of[info.u]}, v{g.rep_of[info.v]};

        // helper: because u merged with v and v merged with x, x
        // merged with u
        auto merge_3way = [&](int u, int v, int x) -> Clause* {
            if (s.value(vars[u][x]) == l_True)
                return NO_REASON;
            if (u == v)
                return NO_REASON;
            reason.clear();
            reason.push(~Lit(vars[u][v]));
            reason.push(~Lit(vars[v][x]));
            if (g.origmatrix[x].fast_contain(u)) {
                return s.addInactiveClause(reason);
            }
            DO_OR_RETURN(s.enqueueFill(Lit(vars[u][x]), reason));
            return NO_REASON;
        };

        // helper: u is merged with v and v has an edge with x, so
        // u has an edge with x
        auto separate = [&](int u, int v, int x) -> Clause* {
            if (g.origmatrix[u].contain(x))
                return NO_REASON;
            if (u == v)
                return NO_REASON;
            assert(u != x && v != x);
            reason.clear();
            reason.push(~Lit(vars[u][v]));
            if (!g.origmatrix[v].fast_contain(x))
                reason.push(Lit(vars[v][x]));
            if (g.origmatrix[u].fast_contain(x))
                return s.addInactiveClause(reason);
            DO_OR_RETURN(s.enqueueFill(~Lit(vars[u][x]), reason));
            return NO_REASON;
        };

        if (!sign(l)) {
            // merging u and v

            if (u == v) {
                // already merged
                return NO_REASON;
            }

            // before merging, keep track of the differences in neighborhood
            diffuv.copy(g.matrix[u]);
            diffuv.setminus_with(g.matrix[v]);
            diffvu.copy(g.matrix[v]);
            diffvu.setminus_with(g.matrix[u]);

            // first, make sure u and v are merged
            DO_OR_RETURN(merge_3way(info.u, info.v, v));
            DO_OR_RETURN(merge_3way(v, info.u, u));

            // merge each partition with the other partition's
            // representative
            for (auto vp : g.partition[v]) {
                if (vp == v)
                    continue;
                DO_OR_RETURN(merge_3way(u, v, vp));
            }

            for (auto up : g.partition[u]) {
                if (up == u)
                    continue;
                DO_OR_RETURN(merge_3way(v, u, up));
            }

            // then merge all the rest, using the fact that each is
            // merged with the representative
            for (auto vp : g.partition[v]) {
                if (vp == v)
                    continue;

                for (auto up : g.partition[u]) {
                    if (up == u)
                        continue;
                    DO_OR_RETURN(merge_3way(vp, u, up));
                }
            }

            // now separate each vertex from potential new
            // neighbors. Here, we only use the original neighborhood
            // of the representative of each partition to determine
            // whether to generate a binary or ternary clause. We
            // could try to find a vertex in the partition that was
            // originally adjacent to the new neighbor, but it's too
            // much work (linear, even if it's mediated by a bitset)

            for (auto vp : g.partition[v]) {
                for (auto n : diffuv)
                    DO_OR_RETURN(separate(vp, u, n));
            }
            for (auto up : g.partition[u]) {
                for (auto n : diffvu)
                    DO_OR_RETURN(separate(up, v, n));
            }

            g.merge(u, v);
        } else {
            if (g.matrix[u].fast_contain(v))
                return NO_REASON;

            // separate u and v first
            DO_OR_RETURN(separate(u, info.u, info.v));
            DO_OR_RETURN(separate(v, info.v, u));

            // separate each vertex in each partition from the
            // representative of the other partition
            for (auto up : g.partition[u])
                DO_OR_RETURN(separate(up, u, v));
            for (auto vp : g.partition[v])
                DO_OR_RETURN(separate(vp, v, u));

            // and now all-to-all
            for (auto up : g.partition[u])
                for (auto vp : g.partition[v])
                    DO_OR_RETURN(separate(up, u, vp));

            g.separate(u, v);
        }

        return NO_REASON;
    }

    Clause* explain_naive_positive()
    {
        auto maxidx = std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques));
        culprit.clear();
        std::copy(begin(cf.cliques[maxidx]), end(cf.cliques[maxidx]),
            back_inserter(culprit));
        reason.clear();
        for (size_t i = 0; i != culprit.size() - 1; ++i)
            for (size_t j = i + 1; j != culprit.size(); ++j) {
                auto u = culprit[i], v = culprit[j];
                assert(g.rep_of[u] == u);
                assert(g.rep_of[v] == v);
                if (!g.origmatrix[u].fast_contain(v))
                    reason.push(Lit(vars[u][v]));
            }
        return s.addInactiveClause(reason);
    }

    // generate an explanation for u having neighborhood N. Appends
    // explanation to reason.
    void explain_N_naive(int u, const bitset& N, vec<Lit>& reason)
    {
        g.util_set.clear();
        for (auto v : g.partition[u])
            g.util_set.union_with(g.origmatrix[v]);
        g.util_set.intersect_with(N);
        expl_residue.copy(g.util_set);
        for (auto v : g.partition[u]) {
            g.util_set.setminus_with(g.origmatrix[v]);
            if (v != u)
                reason.push(~Lit(vars[u][v]));
            if (g.util_set.empty())
                break;
        }
        g.util_set.copy(g.matrix[u]);
        g.util_set.setminus_with(expl_residue);
        g.util_set.intersect_with(N);
        for (auto v : g.util_set) {
            reason.push(Lit(vars[u][v]));
        }
    }

    Clause* explain_naive()
    {
        auto maxidx = std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques));
        reason.clear();
        bitset& covered = expl_covered;
        covered.clear();
        for (auto u : cf.cliques[maxidx]) {
            expl_N.copy(cf.cliques[maxidx]);
            expl_N.setminus_with(covered);
            explain_N_naive(u, expl_N, reason);
            covered.fast_add(u);
        }
        return s.addInactiveClause(reason);
    }

    Clause* explain()
    {
        switch (opt.learning) {
        case options::NO_LEARNING:
            return INVALID_CLAUSE;
        case options::NAIVE_POSITIVE:
            return explain_naive_positive();
        case options::NAIVE:
            return explain_naive();
        default:
            assert(0);
            return INVALID_CLAUSE;
        }
    }

    Clause* propagate(Solver&) final
    {
        int lb{0};

        // recompute the degenracy order
        if (opt.ordering == options::DYNAMIC_DEGENERACY) {
            heuristic.clear();
            adjacency_list.get_degeneracy_order(heuristic);
            std::reverse(heuristic.begin(), heuristic.end());
            lb = cf.find_cliques(heuristic);
        } else if (opt.ordering == options::DEGENERACY
            or opt.ordering == options::INVERSE_DEGENERACY) {
            if (degeneracy_order.empty()) {
                adjacency_list.get_degeneracy_order(degeneracy_order);
                if (opt.ordering == options::INVERSE_DEGENERACY)
                    std::reverse(
                        degeneracy_order.begin(), degeneracy_order.end());
            }
            heuristic.clear();
            for (auto v : degeneracy_order)
                if (g.nodeset.fast_contain(v))
                    heuristic.push_back(v);
            lb = cf.find_cliques(heuristic);
        } else if (opt.ordering == options::PARTITION) {

            // sort by partition size
            heuristic.clear();
            for (auto v : g.nodes)
                heuristic.push_back(v);

            std::sort(heuristic.begin(), heuristic.end(),
                [&](const int x, const int y) {
                    return (g.partition[x].size() > g.partition[y].size());
                });

            lb = cf.find_cliques(heuristic);
        } else {
            // no ordering
            lb = cf.find_cliques(g.nodes);
        }

				if(opt.boundalg == options::FULLMYCIELSKI) {
						auto mlb{mf.full_myciel()};
						std::cout << mlb-lb << std::endl;
				} else if(opt.boundalg == options::MAXMYCIELSKI) {
						auto mlb{mf.improve_cliques_larger_than(lb)};
						std::cout << mlb-lb << std::endl;
				} else if(opt.boundalg == options::GREEDYMYCIELSKI) {
						auto mlb{mf.improve_greedy(lb-1)};
						std::cout << mlb-lb << std::endl;
				}


        if (s.decisionLevel() == 0 && lb > bestlb) {
            bestlb = lb;
            std::cout << "c new lower bound " << bestlb
                      << " time = " << minicsp::cpuTime()
                      << " conflicts = " << s.conflicts << std::endl;
        }
        // std::cout << "lb = " << lb << " num_cliques = " << cf.num_cliques
        //           << " |V| = " << g.nodes.size
        //           << " dlvl = " << s.decisionLevel() << "\n";
        if (cf.num_cliques == 1)
            assert(g.nodes.size() == cf.cliques[0].size());
        if (lb >= ub)
            return explain();
        return NO_REASON;
    }

    // Clause* prune_included_neighborhood(const int lb)
    // {
    //     maxcliques.clear();
    //     for (auto cl = 0; cl < cf.num_cliques; ++cl) {
    //         if (cf.clique_sz[cl] == lb) {
    //             maxcliques.push_back(cl);
    //         }
    //     }
    //     n_included.resize(maxcliques.size());
    //     for (unsigned i = 0; i < maxcliques.size(); ++i) {
    //         n_included[i].clear();
    //     }
    //     for (auto v : g.nodes) {
    //         for (unsigned i = 0; i < maxcliques.size(); ++i) {
    //             cur_neighbors.copy(g.matrix[v]);
    //             cur_neighbors.intersect_with(g.nodeset);
    //             if (cf.cliques[maxcliques[i]].intersect(cur_neighbors)) {
    //                 n_included[i].push_back(v);
    //             }
    //         }
    //     }
    //     for (unsigned i = 0; i < maxcliques.size(); ++i) {
    //         for (unsigned a = 0; a < n_included[i].size(); ++a) {
    //             for (unsigned b = a + 1; b < n_included[i].size(); ++b) {
    //                 if (!g.matrix[n_included[i][a]].fast_contain(
    //                         n_included[i][b])
    //                     && s.value(vars[n_included[i][a]][n_included[i][b]])
    //                         == l_Undef) {
    //                     cur_neighbors.copy(g.matrix[n_included[i][a]]);
    //                     cur_neighbors.union_with(g.matrix[n_included[i][b]]);
    //                     cur_neighbors.intersect_with(g.nodeset);
    //                     if (cur_neighbors.includes(cf.cliques[maxcliques[i]])) {
    //                         std::cout << "pruning " << n_included[i][a] << ","
    //                                   << n_included[i][b] << " ("
    //                                   << ++n_prunings << ")\n";
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     return NO_REASON;
    // }

    void check_consistency()
    {
        g.check_consistency();

        for (auto v : g.nodes) {
            assert(!g.partition[v].empty());
            assert(g.partition[v][0] == v);
            bitset bs(g.matrix[v]);
            for (auto u : g.partition[v]) {
                assert(g.rep_of[u] == v);
                bs.setminus_with(g.origmatrix[u]);
            }
            for (auto u : bs) {
                assert(vars[v][u] != var_Undef);
                if (s.value(vars[v][u]) != l_False)
                    std::cout << "Partition " << v
                              << " = "
                              // << print_container(g.partition[v])
                              << " has extra neighbor " << u << std::endl;
            }
            for (auto u : bs)
                assert(s.value(vars[v][u]) == l_False);
        }

        bool failed{false};
        for (auto& vec : vars)
            for (auto x : vec) {
                if (x == var_Undef)
                    continue;
                auto info = varinfo[x];
                if (s.value(x) == l_False) {
                    if (!g.matrix[g.rep_of[info.u]].contain(g.rep_of[info.v])) {
                        std::cout << info.u << "-" << info.v
                                  << " not connected, var is false\n";
                        failed = true;
                    }
                    if (!g.matrix[g.rep_of[info.v]].contain(g.rep_of[info.u])) {
                        std::cout << info.v << "-" << info.u
                                  << " not connected, var is false\n";
                        failed = true;
                    }
                }
                if (s.value(x) == l_True) {
                    if (g.rep_of[info.u] != g.rep_of[info.u]) {
                        std::cout << info.u << "-" << info.v
                                  << " not merged, var is true\n";
                        failed = true;
                    }
                }
            }
        if (failed)
            std::cout << "oops" << std::endl;
        assert(!failed);
    }
};

cons_base* post_gc_constraint(Solver& s, graph& g,
    const std::vector<std::vector<Var>>& vars, const options& opt, statistics& stat)
{
    auto cons = new gc_constraint(s, g, vars, opt, stat);
    s.addConstraint(cons);
    return cons;
}

} // namespace gc
