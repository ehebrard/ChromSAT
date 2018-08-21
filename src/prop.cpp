#include "prop.hpp"
#include "minicsp/core/utils.hpp"
#include "utils.hpp"
#include "statistics.hpp"
#include "mycielski.hpp"

#include <cmath>

#define NICE_TRACE

namespace gc
{

using namespace minicsp;

dense_graph cons_base::create_filled_graph(
    boost::optional<std::vector<std::pair<int, int>>> fillin)
{
    assert((fillin && opt.fillin) || (!fillin && !opt.fillin));
    if (opt.fillin) {
        dense_graph filled{g};
        for (auto e : *fillin)
            filled.add_edge(e.first, e.second);
        return filled;
    } else {
        dense_graph g;
        return g;
    }
}

class gc_constraint : public minicsp::cons, public cons_base
{
public:
    mycielskan_subgraph_finder<bitset> mf;

    const varmap& vars;
    std::vector<indset_constraint> isconses;

    statistics& stat;

    struct varinfo_t {
        int u{-1}, v{-1};
    };
    // indexed by varid
    std::vector<varinfo_t> varinfo;

    // for the adaptive bound policy: this is set to true on conflicts
    // and reset to false after the stronger policy runs
    bool run_expensive_bound{false};
    minicsp::Solver::clause_callback_t adaptive_callback;
    options::dual_policy bound_source;

    // the usual assortment of vecs, vectors and bitsets to avoid
    // reallocations
    vec<Lit> reason;
    std::vector<int> culprit;
    bitset diffuv, diffvu;

    bitset util_set;
    bitset neighborhood;
    bitset ordering_tmp, ordering_forbidden, ordering_removed_bs;

    // the representative chosen for each partition during explanation
    std::vector<int> expl_reps;

    std::vector<int> ordering_removed;

    std::vector<int> heuristic;

public:
    gc_constraint(Solver& solver, dense_graph& g,
        boost::optional<std::vector<std::pair<int, int>>> fillin,
        const varmap& tvars, const std::vector<indset_constraint>& isconses,
        const options& opt, statistics& stat)
        : cons_base(solver, opt, g, fillin)
        , mf(g, cf, opt.prune)
        , vars(tvars)
        , isconses(isconses)
        , stat(stat)
        , util_set(0, g.capacity() - 1, bitset::empt)
        , neighborhood(0, g.capacity() - 1, bitset::empt)
        , ordering_tmp(0, g.capacity() - 1, bitset::empt)
        , ordering_forbidden(0, g.capacity() - 1, bitset::empt)
        , ordering_removed_bs(0, g.capacity() - 1, bitset::empt)
    {
        stat.binds(this);
        ub = g.capacity();
        for (int i = 0; i != g.capacity(); ++i) {
            if (!g.nodes.contain(i))
                continue;
            for (int j = 0; j != g.capacity(); ++j) {
                if (!g.nodes.contain(j))
                    continue;
                if (j < i) {
                    assert(vars[i][j] == vars[j][i]);
                    continue;
                }
                if (j == i)
                    continue;
                auto ijvar = vars[i][j];
                if (ijvar == minicsp::var_Undef)
                    continue;
                if (varinfo.size() <= static_cast<size_t>(ijvar))
                    varinfo.resize(ijvar + 1);
                varinfo[ijvar] = {i, j};
                s.wake_on_lit(ijvar, this, nullptr);
                s.schedule_on_lit(ijvar, this);
            }
        }
        diffuv.initialise(0, g.capacity(), bitset::empt);
        diffvu.initialise(0, g.capacity(), bitset::empt);
        // cur_neighbors.initialise(0, g.capacity(), bitset::empt);

        if (opt.adaptive) {
            adaptive_callback
                = [this](const vec<Lit>&,
                      int) -> minicsp::Solver::clause_callback_result_t {
                this->run_expensive_bound = true;
                return minicsp::Solver::CCB_OK;
            };
            s.use_clause_callback(adaptive_callback);
        }

        // TODO
        // std::cout << "UB = " << ub << std::endl;
        // int k = ub-1;
        // int a = g.capacity() / k;
        //
        //
        //
        //                        mineq_constraint = new cons_pb(s, vars,
        // weights,
        // lb);
        //                        s.addConstraint(c);
        //
        // post_pb(Solver& s, std::vector<Var> const& vars,
        //               std::vector<int> const& weights, int lb);
        // int mineq = a * (g.capacity() - k * (a + 1) / 2);
        //                        mineq_constraint = new cons_pbvar(s, vbool, c,
        // rhs);
        //                        s.addConstraint(mineq_constraint);

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

    Clause* wake(Solver& s, Lit l)
    {
        sync_graph();

        auto info = varinfo[var(l)];
        auto u{g.rep_of[info.u]}, v{g.rep_of[info.v]};

        // helper: because u merged with v and v merged with x, x
        // merged with u
        auto merge_3way = [&](int u, int v, int x) -> Clause* {
            Var uxvar = vars[u][x];
            // if (opt.fillin && uxvar == var_Undef)
            //     return NO_REASON;
            if (s.value(uxvar) == l_True)
                return NO_REASON;
            if (u == v)
                return NO_REASON;
            if (!vars.contain(u, x))
                return NO_REASON;

#ifdef NICE_TRACE
            std::cout << " -> merge " << u << " with "
                      << " " << x << std::endl;
#endif

            reason.clear();
            reason.push(~Lit(vars[u][v]));
            reason.push(~Lit(vars[v][x]));
            if (g.origmatrix[x].fast_contain(u)) {
                return s.addInactiveClause(reason);
            }
            DO_OR_RETURN(s.enqueueFill(Lit(uxvar), reason));
            return NO_REASON;
        };

        // helper: u is merged with v and v has an edge with x, so
        // u has an edge with x
        auto separate = [&](int u, int v, int x) -> Clause* {
            Var uxvar = vars[u][x];
            if (opt.fillin && uxvar == var_Undef)
                return NO_REASON;
            if (g.origmatrix[u].fast_contain(x))
                return NO_REASON;
            if (u == v)
                return NO_REASON;
            if (!vars.contain(u, x))
                return NO_REASON;

#ifdef NICE_TRACE
            std::cout << " -> separate " << u << " from "
                      << " " << x << std::endl;
#endif

            assert(u != x && v != x);
            reason.clear();
            reason.push(~Lit(vars[u][v]));
            if (!g.origmatrix[v].fast_contain(x)) {
                reason.push(Lit(vars[v][x]));
            }
            DO_OR_RETURN(s.enqueueFill(~Lit(vars[u][x]), reason));
            return NO_REASON;
        };

        if (!sign(l)) {
            // merging u and v

#ifdef NICE_TRACE
            std::cout << "-merge " << u << " and " << v << std::endl;
#endif

            if (u == v) {
                // already merged
                return NO_REASON;
            }

            // before merging, keep track of the differences in neighborhood
            diffuv.copy(g.matrix[u]);
            diffuv.setminus_with(g.matrix[v]);
            diffvu.copy(g.matrix[v]);
            diffvu.setminus_with(g.matrix[u]);

            // merge each partition with the other partition's
            // representative
            for (auto vp : g.partition[v]) {
                if (vp == info.v)
                    continue;

                if (opt.fillin && !fg.matrix[vp].fast_contain(info.u))
                    continue;
                DO_OR_RETURN(merge_3way(info.u, info.v, vp));
            }

            for (auto up : g.partition[u]) {
                if (up == info.u)
                    continue;
                if (opt.fillin && !fg.matrix[up].fast_contain(info.v))
                    continue;
                DO_OR_RETURN(merge_3way(info.v, info.u, up));
            }

            // then merge all the rest, using the fact that each is
            // merged with the representative
            for (auto vp : g.partition[v]) {
                if (vp == info.v)
                    continue;

                for (auto up : g.partition[u]) {
                    if (up == info.u)
                        continue;
                    DO_OR_RETURN(merge_3way(vp, info.u, up));
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
                if (opt.fillin) {
                    util_set.copy(diffuv);
                    util_set.intersect_with(fg.matrix[vp]);
                    for (auto n : util_set)
                        DO_OR_RETURN(separate(vp, info.u, n));
                } else {
                    for (auto n : diffuv)
                        DO_OR_RETURN(separate(vp, info.u, n));
                }
            }
            for (auto up : g.partition[u]) {
                if (opt.fillin) {
                    util_set.copy(diffvu);
                    util_set.intersect_with(fg.matrix[up]);
                    for (auto n : util_set)
                        DO_OR_RETURN(separate(up, info.v, n));
                } else {
                    for (auto n : diffvu)
                        DO_OR_RETURN(separate(up, info.v, n));
                }
            }

            g.merge(u, v);
        } else {

#ifdef NICE_TRACE
            std::cout << "-separate " << u << " and " << v << std::endl;
#endif

            if (g.matrix[u].fast_contain(v))
                return NO_REASON;

            // separate each vertex in each partition from the
            // representative of the other partition
            for (auto up : g.partition[u])
                DO_OR_RETURN(separate(up, info.u, info.v));
            for (auto vp : g.partition[v])
                DO_OR_RETURN(separate(vp, info.v, info.u));

            // and now all-to-all
            for (auto up : g.partition[u])
                for (auto vp : g.partition[v])
                    DO_OR_RETURN(separate(up, info.u, vp));

            g.separate(u, v);
        }

        return NO_REASON;
    }

    // fills in reason with explanation for clq. Also generates
    // this->culprit and fills in this->expl_reps
    void explain_positive_clique(const bitset& clq)
    {
        culprit.clear();
        std::copy(begin(clq), end(clq), back_inserter(culprit));
        if (opt.learning == options::CHOOSE_POSITIVE) {
            // sort by increasing partition size
            std::sort(begin(culprit), end(culprit), [&](auto u, auto v) {
                return g.partition[u].size() < g.partition[v].size();
            });
        }

        assert(!culprit.empty());

        expl_reps.clear();
        expl_reps.resize(g.capacity(), -1);

        double maxactivity{0.0};
        Var maxvar = var_Undef;
        // returns true to stop
        auto bestmatch = [&](auto u, auto v, auto& ubag, auto& vbag, auto up) {
            if (opt.learning != options::CHOOSE_POSITIVE) {
                expl_reps[u] = g.rep_of[u];
                expl_reps[v] = g.rep_of[v];
                return true;
            }
            int bestu{-1}, bestv{-1};
            for (auto vp : vbag) {
                if (g.origmatrix[up].fast_contain(vp)) {
                    expl_reps[u] = up;
                    expl_reps[v] = vp;
                    maxvar = var_Undef;
                    return true;
                } else {
                    auto var = vars[up][vp];
                    auto varact = s.var_activity(var);
                    assert(s.value(var) == l_False);
                    if (varact > maxactivity) {
                        maxactivity = varact;
                        maxvar = var;
                        bestu = up;
                        bestv = vp;
                    }
                }
            }
            if (bestu >= 0) {
                expl_reps[u] = bestu;
                expl_reps[v] = bestv;
            }
            return false;
        };

        for (size_t i = 0; i != culprit.size() - 1; ++i)
            for (size_t j = i + 1; j != culprit.size(); ++j) {
                auto u = culprit[i], v = culprit[j];
                assert(g.rep_of[u] == u);
                assert(g.rep_of[v] == v);
                assert(g.matrix[u].fast_contain(v));
                assert(g.matrix[v].fast_contain(u));
                auto ubag = &g.partition[u];
                auto vbag = &g.partition[v];
                maxactivity = 0.0;
                maxvar = var_Undef;
                if (expl_reps[v] >= 0 && expl_reps[u] < 0) {
                    using std::swap;
                    swap(u, v);
                    swap(ubag, vbag);
                }
                if (expl_reps[v] < 0 && expl_reps[u] >= 0) {
                    // find a rep for vbag only
                    auto up = expl_reps[u];
                    bestmatch(u, v, *ubag, *vbag, up);
                } else if (expl_reps[v] < 0 && expl_reps[u] < 0) {
                    // find a rep for both vbag and ubag
                    for (auto up : *ubag) {
                        if (bestmatch(u, v, *ubag, *vbag, up))
                            break;
                    }
                } else {
                    auto ur = expl_reps[u];
                    auto vr = expl_reps[v];
                    maxvar = vars[ur][vr];
                    assert(maxvar == var_Undef || s.value(maxvar) == l_False);
                }
                if (maxvar != var_Undef) {
                    assert(s.value(maxvar) == l_False);
                    reason.push(Lit(maxvar));
                }
            }
    }

    Clause* explain_positive()
    {
        auto maxidx{std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques))};

        if (bound_source != options::CLIQUES && mf.explanation_clique != -1) {
            maxidx = mf.explanation_clique;
        }

        // explain the base clique
        reason.clear();
        explain_positive_clique(cf.cliques[maxidx]);

        if (bound_source != options::CLIQUES && mf.explanation_clique != -1) {
            // for (auto v : mf.explanation_subgraph.nodes) {
            for (auto i{culprit.size()};
                 i < mf.explanation_subgraph.nodes.size(); ++i) {
                auto v{mf.explanation_subgraph.nodes[i]};
                auto vr = expl_reps[v];
                if (vr < 0) {
                    expl_reps[v] = v;
                    vr = v;
                }
                neighborhood.copy(mf.explanation_subgraph.matrix[v]);
                for (auto u : neighborhood) {
                    auto ur = expl_reps[u];
                    if (ur < 0) {
                        expl_reps[u] = u;
                        ur = u;
                    }
                    if (mf.explanation_subgraph.nodes.index(v)
                        > mf.explanation_subgraph.nodes.index(u)) {
                        assert(g.rep_of[u] == u);
                        assert(g.rep_of[v] == v);
                        if (!g.origmatrix[ur].fast_contain(vr))
                            reason.push(Lit(vars[ur][vr]));
                    }
                }
            }
        }

        return s.addInactiveClause(reason);
    }

    Clause* explain()
    {
        switch (opt.learning) {
        case options::NO_LEARNING:
            return INVALID_CLAUSE;
        case options::NAIVE_POSITIVE:
        case options::CHOOSE_POSITIVE:
            return explain_positive();
        default:
            assert(0);
            return INVALID_CLAUSE;
        }
    }

    Clause* contractWhenNIncluded()
    {
        bool some_propagation = true;
        while (some_propagation) {
            some_propagation = false;
            for (auto u : g.nodes) {
                for (auto v : g.nodes) {
                    if (u != v && !g.origmatrix[u].fast_contain(v)
                        && s.value(vars[u][v]) == l_Undef) {
                        neighborhood.copy(g.matrix[v]);
                        neighborhood.setminus_with(g.matrix[u]);
                        if (!neighborhood.intersect(g.nodeset)) {
                            // N(v) <= N(U)s
                            some_propagation = true;
                            ++stat.num_neighborhood_contractions;

                            // explanation: the edges (ON(v) \ ON(u)) x u
                            neighborhood.copy(g.origmatrix[v]);
                            neighborhood.setminus_with(g.origmatrix[u]);
                            neighborhood.intersect_with(g.nodeset);

                            reason.clear();
                            for (auto w : neighborhood) {
                                reason.push(Lit(vars[u][w]));
                            }
                            DO_OR_RETURN(
                                s.enqueueFill(Lit(vars[u][v]), reason));
                        }
                    }
                }
            }
        }
        return NO_REASON;
    }

    // find an upper bound u for the maximum IS. n/ceil(u) is a lower
    // bound for the chromatic number
    Clause* propagate_is()
    {
        create_ordering();
        ccf.find_clique_cover(heuristic);
        if (std::ceil(g.nodes.size() / static_cast<double>(ccf.num_cliques))
            >= ub) {
            std::cout << "\tIS ub = " << ccf.num_cliques
                      << " |V| = " << g.nodes.size() << " coloring lb = "
                      << std::ceil(g.nodes.size()
                             / static_cast<double>(ccf.num_cliques))
                      << "\n";
            reason.clear();
            for (int i = 0; i != ccf.num_cliques; ++i)
                explain_positive_clique(ccf.cliques[i]);
            return s.addInactiveClause(reason);
        }
        return NO_REASON;
    }

    void create_ordering()
    {
        // recompute the degenracy order
        if (opt.ordering == options::DYNAMIC_DEGENERACY) {
            assert(false);
        } else if (opt.ordering == options::DEGENERACY
            or opt.ordering == options::INVERSE_DEGENERACY) {
            assert(false);
        } else if (opt.ordering == options::PARTITION) {

            if (opt.ordering_low_degree == options::PREPROCESSING_ORDERING) {
                bool removed_some{false};
                ordering_removed_bs.clear();
                ordering_removed.clear();
                do {
                    removed_some = false;
                    ordering_forbidden.clear();
                    for (auto v : g.nodes) {
                        if (ordering_forbidden.fast_contain(v)
                            || ordering_removed_bs.fast_contain(v))
                            continue;
                        ordering_tmp.copy(g.matrix[v]);
                        ordering_tmp.intersect_with(g.nodeset);
                        if (ordering_tmp.size()
                            >= static_cast<size_t>(std::max(bestlb, *lastlb)))
                            continue;
                        removed_some = true;
                        ordering_removed_bs.fast_add(v);
                        ordering_removed.push_back(v);
                        ordering_forbidden.union_with(g.matrix[v]);
                        ordering_forbidden.fast_add(v);
                    }
                } while (removed_some);
            } else if (opt.ordering_low_degree == options::DEGREE_ORDERING) {
                ordering_removed_bs.clear();
                ordering_removed.clear();
                for (auto v : g.nodes) {
                    ordering_tmp.copy(g.matrix[v]);
                    ordering_tmp.intersect_with(g.nodeset);
                    if (ordering_tmp.size()
                        >= static_cast<size_t>(std::max(bestlb, *lastlb)))
                        continue;
                    ordering_removed_bs.fast_add(v);
                    ordering_removed.push_back(v);
                }
            }

            // sort by partition size
            heuristic.clear();
            for (auto v : g.nodes)
                if (!ordering_removed_bs.fast_contain(v))
                    heuristic.push_back(v);

            std::sort(heuristic.begin(), heuristic.end(),
                [&](const int x, const int y) {
                    return (g.partition[x].size() > g.partition[y].size()
                        || (g.partition[x].size() == g.partition[y].size()
                               && x < y));
                });

            std::reverse(begin(ordering_removed), end(ordering_removed));
            for (auto v : ordering_removed)
                heuristic.push_back(v);

            assert(heuristic.size() == g.nodes.size());
        } else {
            // no ordering
            heuristic.clear();
            std::copy(
                g.nodes.begin(), g.nodes.end(), std::back_inserter(heuristic));
        }
    }

    Clause* propagate(Solver&) final
    {
        int lb{0};

        bound_source = options::CLIQUES;

        create_ordering();

        lb = cf.find_cliques(heuristic, opt.cliquelimit);

        if (lb < ub
            && (s.decisionLevel() == 0 || !opt.adaptive
                   || run_expensive_bound)) {
            run_expensive_bound = false;
            bound_source = opt.boundalg;
            auto mlb{lb};
            if (opt.boundalg == options::FULLMYCIELSKI) {
                mlb = mf.full_myciel(lb, ub, s, vars);
            } else if (opt.boundalg == options::MAXMYCIELSKI) {
                mlb = mf.improve_cliques_larger_than(lb, lb, ub, s, vars);
            } else if (opt.boundalg == options::GREEDYMYCIELSKI) {
                mlb = mf.improve_greedy(lb - 1, lb, ub, s, vars);
            }

            // std::cout << lb << " --> " << mlb << std::endl;
            stat.notify_bound_delta(lb, mlb);

            // if(stat.num_bound_delta%100 == 0)
            //              stat.display(std::cout);

            lb = mlb;
        }

        *lastlb = lb;

        bool use_global_bound = false;
        if (lb <= bestlb) {
            lb = bestlb;
            use_global_bound = true;
        }

        if (s.decisionLevel() == 0 && lb > bestlb) {
            bestlb = lb;
            stat.display(std::cout);
            // for (auto v : g.nodes) {
            //     if (g.matrix[v].size() < static_cast<size_t>(bestlb)) {
            //         ++stat.num_vertex_removals;
            //     }
            // }
        }
        if (cf.num_cliques == 1)
            assert(g.nodes.size() == cf.cliques[0].size());
        if (lb >= ub) {
            if (use_global_bound) {
                reason.clear();
                return s.addInactiveClause(reason);
            }
            return explain();
        }

        if (opt.indset_lb)
            DO_OR_RETURN(propagate_is());

        // check local constraints

#ifdef DEBUG_IS
        std::cout << "propag: " << lb << " " << actualub << std::endl;
#endif

#ifdef DEBUG_IS
        for (const auto& c : isconses) {
            std::cout << " cons " << c.vs << std::endl;
        }
#endif

        if (lb >= actualub - 1 and lb >= ub - 1) {
            for (int i = 0; i != cf.num_cliques; ++i) {

#ifdef DEBUG_IS
                std::cout << "clique " << cf.cliques[i] << " in " << g.nodeset
                          << std::endl;
#endif

                if (cf.clique_sz[i] < ub -1)
                    continue;

                for (const auto& c : isconses) {

#ifdef DEBUG_IS
                    std::cout << " cons " << c.vs << std::endl;
#endif

                    util_set.clear();
                    for (auto v : c.vs)
                        util_set.fast_add(g.rep_of[v]);
                    util_set.intersect_with(cf.cliques[i]);

#ifdef DEBUG_IS
                    std::cout << " inter = " << util_set << std::endl;
#endif

                    if (static_cast<int>(util_set.size()) >= actualub - 1) {
                        reason.clear();
                        explain_positive_clique(util_set);
                        return s.addInactiveClause(reason);
                    }
                }
            }
        }

        if (opt.dominance)
            return contractWhenNIncluded();

        return NO_REASON;
    }

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
        for (auto& vec : vars.vars)
            for (auto p : vec.second) {
                auto x = p.second;
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

void union_with_sets(
    const std::vector<int>& ind, const std::vector<bitset>& bv, bitset& result)
{
    for (auto i : ind)
        result.union_with(bv[i]);
}

bool intersect_vec_bs_p(const std::vector<int> vec, const bitset& bs)
{
    return std::any_of(
        begin(vec), end(vec), [&](int up) { return bs.fast_contain(up); });
};

// pick which partition will be our w
int pick_partition(const dense_graph& g, const bitset& clq, bitset& util_set,
    const std::vector<std::vector<int>>& partitions,
    const std::vector<int>& revmap)
{
    auto numcovered = [&](int v) {
        util_set.clear();
        union_with_sets(partitions[revmap[v]], g.origmatrix, util_set);
        // std::cout << "N(" << v << ") =" << util_set << "\n";
        return std::count_if(begin(clq), end(clq), [&](int u) {
            return u != v
                && intersect_vec_bs_p(partitions[revmap[u]], util_set);
        });
    };

    std::vector<int> nc(g.capacity());
    for (auto v : clq) {
        nc[v] = numcovered(v);
        // std::cout << "nc[" << v << "] = " << nc[v] << "\n";
    }
    return *std::min_element(begin(clq), end(clq), [&](int u, int v) {
        return std::pair<int, size_t>{-nc[u], partitions[revmap[u]].size()}
        < std::pair<int, size_t>{-nc[v], partitions[revmap[v]].size()};
    });
}

// pick a neighbor of w to remove from each partition of clq
void update_partitions(const dense_graph& g,
    std::vector<std::vector<int>>& partitions,
    std::vector<int>& covered_neighbors, const bitset& clq, int w,
    const std::vector<int>& revmap, bitset& util_set)
{
    util_set.clear();
    union_with_sets(partitions[revmap[w]], g.origmatrix, util_set);
    covered_neighbors.clear();
    for (auto u : clq) {
        auto i = std::find_if(begin(partitions[revmap[u]]),
            end(partitions[revmap[u]]),
            [&](int x) { return util_set.fast_contain(x); });
        if (i != end(partitions[revmap[u]])) {
            int x = *i;
            if (partitions[revmap[u]].size() > 1)
                partitions[revmap[u]].erase(i);
            covered_neighbors.push_back(x);
        }
    }
}

void remap_constraint(
    bitset& bs, bitset& util_set, const std::vector<int>& vertex_map)
{
    util_set.clear();
    for (auto v : bs)
        util_set.fast_add(vertex_map[v]);
    bs.copy(util_set);
}

cons_base* post_gc_constraint(Solver& s, dense_graph& g,
    boost::optional<std::vector<std::pair<int, int>>> fillin,
    const varmap& vars, const std::vector<indset_constraint>& isconses,
    const std::vector<int>& vertex_map, const options& opt, statistics& stat)
{
    // try {
    if (g.size() > 0) {
        auto cons = new gc_constraint(s, g, fillin, vars, isconses, opt, stat);
        s.addConstraint(cons);
        for (auto& c : cons->isconses)
            remap_constraint(c.vs, cons->util_set, vertex_map);
        return cons;
    }
    // } catch (minicsp::unsat &u) {
    return NULL;
}
} // namespace gc
