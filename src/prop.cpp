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

    mycielskan_subgraph_finder mf;

    const std::vector<std::vector<Var>>& vars;
    const options& opt;
    statistics& stat;

    struct varinfo_t {
        int u{-1}, v{-1};
    };
    // indexed by varid
    std::vector<varinfo_t> varinfo;

    // something which is kept in sync by the solver. if it diverges
    // from what the graph thinks, it means we have
    // backtracked/restarted, so we should resync to that point
    backtrackable<int> lastdlvl;

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
    bitset expl_N, expl_covered, expl_residue;
    bitset expl_clqcopy;
    bitset neighborhood;

    std::vector<int> myc_reason;
    // bitset count;
    std::vector<int> expl_clq;

    // the partitions involved in a clique
    std::vector<std::vector<int>> expl_partitions;
    // to which partition of a clique does a vertex belong
    std::vector<int> expl_revmap;
    // which vertex have we chosen as reprentative from each partition
    std::vector<int> expl_part_rep;
    // the vertices we keep from each partition as we construct the
    // clever explanation
    std::vector<std::vector<int>> expl_covered_neighbors;
    // the vertex chosen as w at each step
    std::vector<int> expl_myc_w;

    neighbors_wrapper adjacency_list;

    std::vector<int> degeneracy_order;
    std::vector<int> heuristic;

    // bitset cur_neighbors;
    // std::vector<int> maxcliques;
    // std::vector<std::vector<int>> n_included;
    // long int n_prunings;

public:
    gc_constraint(Solver& solver, graph& pg,
        const std::vector<std::vector<Var>>& tvars, const options& opt,
        statistics& stat)
        : cons_base(pg)
        , s(solver)
        , g(pg)
        , mf(g, cf, opt.prune)
        , vars(tvars)
        , opt(opt)
        , stat(stat)
        , lastdlvl(s)
        , util_set(0, g.capacity() - 1, bitset::empt)
        , expl_N(0, g.capacity() - 1, bitset::empt)
        , expl_covered(0, g.capacity() - 1, bitset::empt)
        , expl_residue(0, g.capacity() - 1, bitset::empt)
        , expl_clqcopy(0, g.capacity() - 1, bitset::empt)
        , neighborhood(0, g.capacity() - 1, bitset::empt)
        // , global_myciel_clique(-1)
        // , count(0, g.capacity()*g.capacity(), bitset::empt)
        , expl_partitions(g.capacity())
        , expl_revmap(g.capacity())
        , expl_part_rep(g.capacity())
        , expl_covered_neighbors(g.capacity())
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

        if (opt.adaptive) {
            adaptive_callback
                = [this](const vec<Lit>&,
                      int) -> minicsp::Solver::clause_callback_result_t {
                this->run_expensive_bound = true;
                return minicsp::Solver::CCB_OK;
            };
            s.use_clause_callback(adaptive_callback);
        }

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
        auto maxidx{std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques))};

        if (bound_source != options::CLIQUES && mf.explanation_clique != -1) {
            maxidx = mf.explanation_clique;
        }

        // explain the base clique
        reason.clear();

        culprit.clear();
        std::copy(begin(cf.cliques[maxidx]), end(cf.cliques[maxidx]),
            back_inserter(culprit));
        for (size_t i = 0; i != culprit.size() - 1; ++i)
            for (size_t j = i + 1; j != culprit.size(); ++j) {
                auto u = culprit[i], v = culprit[j];
                assert(g.rep_of[u] == u);
                assert(g.rep_of[v] == v);
                if (!g.origmatrix[u].fast_contain(v)) {
                    reason.push(Lit(vars[u][v]));
                }
            }

        if (bound_source != options::CLIQUES && mf.explanation_clique != -1) {
            // for (auto v : mf.explanation_subgraph.nodes) {
            for (auto i{culprit.size()};
                 i < mf.explanation_subgraph.nodes.size(); ++i) {
                auto v{mf.explanation_subgraph.nodes[i]};
                neighborhood.copy(mf.explanation_subgraph.matrix[v]);
                neighborhood.setminus_with(g.origmatrix[v]);
                // neighborhood.set_min(v);
                for (auto u : neighborhood) {
                    if (mf.explanation_subgraph.nodes.index(v)
                        > mf.explanation_subgraph.nodes.index(u)) {
                        assert(g.rep_of[u] == u);
                        assert(g.rep_of[v] == v);
                        reason.push(Lit(vars[u][v]));
                    }
                }
            }
        }

        return s.addInactiveClause(reason);
    }

    // generate an explanation for u having neighborhood N. Appends
    // explanation to reason.
    void explain_N_naive(int u, const bitset& N, vec<Lit>& reason)
    {
        util_set.clear();
        for (auto v : g.partition[u])
            util_set.union_with(g.origmatrix[v]);
        util_set.intersect_with(N);
        expl_residue.copy(util_set);
        for (auto v : g.partition[u]) {
            util_set.setminus_with(g.origmatrix[v]);
            if (v != u)
                reason.push(~Lit(vars[u][v]));
            if (util_set.empty())
                break;
        }
        util_set.copy(g.matrix[u]);
        util_set.setminus_with(expl_residue);
        util_set.intersect_with(N);
        for (auto v : util_set) {
            reason.push(Lit(vars[u][v]));
        }
    }

    // explain the clique using function e to explain each neighborhood
    template <typename F> Clause* explain_clique(F e)
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
            e(u, expl_N, reason);
            covered.fast_add(u);
        }
        return s.addInactiveClause(reason);
    }

    // Mycielskian explanations
    void explain_myc_w_edges(int w, const bitset& clq,
        const std::vector<std::vector<int>>& partitions,
        std::vector<int>& covered_neighbors);
    void explain_myc_u_edges(const bitset& clq,
        const std::vector<std::vector<int>>& partitions,
        std::vector<int>& covered_neighbors);
    void explain_myc_rec(bitset& clq, std::vector<std::vector<int>>& partitions,
        const std::vector<int>& revmap);
    Clause* explain_myc();
    void verify_myc_reason();

    Clause* explain()
    {
        switch (opt.learning) {
        case options::NO_LEARNING:
            return INVALID_CLAUSE;
        case options::NAIVE_POSITIVE:
            return explain_naive_positive();
        case options::NAIVE:
            return explain_clique([&](int u, auto&& N, auto&& reason) {
                explain_N_naive(u, N, reason);
            });
        case options::MYC_POSITIVE:
            return explain_myc();
            break;
        default:
            assert(0);
            return INVALID_CLAUSE;
        }
    }

    Clause* propagate(Solver&) final
    {
        int lb{0};

        bound_source = options::CLIQUES;

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
                    return (g.partition[x].size() > g.partition[y].size()
                        || (g.partition[x].size() == g.partition[y].size()
                               && x < y));
                });

            lb = cf.find_cliques(heuristic);
        } else {
            // no ordering
            lb = cf.find_cliques(g.nodes);
        }

        if (s.decisionLevel() == 0 || !opt.adaptive || run_expensive_bound) {
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
            stat.notify_bound_delta(mlb - lb);
            lb = mlb;
        }

        bool use_global_bound = false;
        if (lb <= bestlb) {
            lb = bestlb;
            use_global_bound = true;
        }

        if (s.decisionLevel() == 0 && lb > bestlb) {
            bestlb = lb;
            std::cout << "c new lower bound " << bestlb
                      << " time = " << minicsp::cpuTime()
                      << " conflicts = " << s.conflicts << std::endl;
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
int pick_partition(const graph& g, const bitset& clq, bitset& util_set,
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
void update_partitions(const graph& g,
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

int pick_vertex(const graph& g, const bitset& clq, int w,
    const std::vector<std::vector<int>>& partitions,
    const std::vector<int>& revmap, bitset& util_set)
{
    auto& part = partitions[revmap[w]];
    assert(!part.empty());
    return part[0];
}

void gc_constraint::explain_myc_w_edges(int w, const bitset& clq,
    const std::vector<std::vector<int>>& partitions,
    std::vector<int>& covered_neighbors)
{
    auto wp = pick_vertex(g, clq, w, partitions, expl_revmap, util_set);
    expl_part_rep[expl_revmap[w]] = wp;
    std::cout << "Picked " << wp << " as rep of partition " << w << "("
              << print_container<std::vector<int>>(partitions[expl_revmap[w]])
              << ")\n";
    util_set.clear();
    // w (as a partition) is adjacent to all vertices in
    // covered_neighbors. wp, which we have chosen as the vertex from
    // w which will represent the partition in the failed graph, may
    // not be adjacent to all of them.
    for (auto u : covered_neighbors) {
        if (!g.origmatrix[wp].fast_contain(u)) {
            std::cout << "1:reason += " << lit_printer(s, Lit(vars[wp][u]))
                      << "\n";
            reason.push(Lit(vars[wp][u]));
        }
        util_set.fast_add(g.rep_of[u]);
    }
    // now for all the other partitions, make sure that wp is
    // connected to whatever vertex remains
    for (auto u : clq) {
        if (expl_revmap[u] == expl_revmap[wp])
            continue;
        auto& part = partitions[expl_revmap[u]];
        assert(!part.empty());
        assert(u == g.rep_of[u]);
        if (util_set.fast_contain(g.rep_of[u]))
            continue;
        auto urep = expl_part_rep[expl_revmap[u]];
        if (urep == g.rep_of[wp])
            continue;
        assert(!g.origmatrix[wp].fast_contain(urep));
        std::cout << "2:reason += " << lit_printer(s, Lit(vars[wp][urep]))
                  << "\n";
        reason.push(Lit(vars[wp][urep]));
    }
}

void gc_constraint::explain_myc_u_edges(const bitset& clq,
    const std::vector<std::vector<int>>& partitions,
    std::vector<int>& covered_neighbors)
{
    for (auto u : covered_neighbors) {
        if (u == expl_part_rep[expl_revmap[u]])
            continue;
        std::cout << "Explaining N(" << u << ") \\subseteq N("
                  << expl_part_rep[expl_revmap[u]] << ")\n";
        for (auto v : clq) {
            auto vrep = expl_part_rep[expl_revmap[v]];
            if (expl_revmap[v] == expl_revmap[u])
                continue;
            if (g.origmatrix[u].fast_contain(vrep))
                continue;
            std::cout << "3:reason += " << lit_printer(s, Lit(vars[u][vrep]))
                      << "\n";
            reason.push(Lit(vars[u][vrep]));
        }
    }
}

void gc_constraint::explain_myc_rec(bitset& clq,
    std::vector<std::vector<int>>& partitions, const std::vector<int>& revmap)
{
    if (clq.size() == 0)
        return;

    for (auto u : clq)
        std::cout << "  partition[" << u << "] = "
                  << print_container<std::vector<int>>(partitions[revmap[u]])
                  << "\n";

    auto w = pick_partition(g, clq, util_set, partitions, revmap);
    clq.fast_remove(w);
    std::cout << "selected partition " << w << " as w\n";
    expl_myc_w.push_back(w);
    update_partitions(g, partitions, expl_covered_neighbors[revmap[w]], clq, w,
        revmap, util_set);
    std::cout << "kept neighbors "
              << print_container<std::vector<int>>(
                     expl_covered_neighbors[revmap[w]])
              << "\n";
    explain_myc_rec(clq, partitions, revmap);

    explain_myc_w_edges(w, clq, partitions, expl_covered_neighbors[revmap[w]]);
    explain_myc_u_edges(clq, partitions, expl_covered_neighbors[revmap[w]]);
    clq.fast_add(w);
    for (auto u : expl_covered_neighbors[revmap[w]])
        partitions[revmap[u]].push_back(u);
}

Clause* gc_constraint::explain_myc()
{
    auto maxidx = std::distance(begin(cf.clique_sz),
        std::max_element(
            begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques));
    auto& clq = expl_clq;
    clq.clear();
    std::copy(begin(cf.cliques[maxidx]), end(cf.cliques[maxidx]),
        std::back_inserter(clq));

    for (size_t i = 0; i != clq.size(); ++i) {
        auto v = clq[i];
        for (auto u : g.partition[v])
            expl_revmap[u] = i;
        expl_partitions[i] = g.partition[v];
    }

    reason.clear();
    expl_myc_w.clear();

    std::cout << "\nconflict " << s.conflicts + 1 << " clique "
              << print_container<std::vector<int>>(clq) << "\n";

    bitset& clqcopy = expl_clqcopy;
    clqcopy.copy(cf.cliques[maxidx]);
    explain_myc_rec(clqcopy, expl_partitions, expl_revmap);

    std::cout << "reason " << minicsp::print(s, &reason) << "\n\n"
              << std::flush;
    verify_myc_reason();
    return s.addInactiveClause(reason);
}

void gc_constraint::verify_myc_reason()
{
    for (Lit l : reason)
        assert(s.value(l) == l_False);

    assert(0);
}

cons_base* post_gc_constraint(Solver& s, graph& g,
    const std::vector<std::vector<Var>>& vars, const options& opt,
    statistics& stat)
{
    auto cons = new gc_constraint(s, g, vars, opt, stat);
    s.addConstraint(cons);
    return cons;
}
} // namespace gc
