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

    clique_finder cf;

public:
    gc_constraint(Solver& solver, graph& pg,
        const std::vector<std::vector<Var>>& tvars, const options& opt)
        : s(solver)
        , g(pg)
        , vars(tvars)
        , opt(opt)
        , lastdlvl(s)
        , cf(g)
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

        DO_OR_THROW(propagate(s));
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
            // merging u and v. Note that (u,v) may be (info.u,
            // info.v) or (info.v, info.u). u is the representative

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

    Clause *explain() {
        auto maxidx = std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques));
        culprit.clear();
        std::copy(begin(cf.cliques[maxidx]), end(cf.cliques[maxidx]),
            back_inserter(culprit));
        reason.clear();
        for(size_t i = 0; i != culprit.size()-1; ++i)
            for(size_t j = i+1; j != culprit.size(); ++j) {
                auto u = culprit[i], v = culprit[j];
                assert(g.rep_of[u]==u);
                assert(g.rep_of[v]==v);
                if (!g.origmatrix[u].fast_contain(v))
                    reason.push(Lit(vars[u][v]));
            }
        return s.addInactiveClause(reason);
    }

    Clause* propagate(Solver&) final
    {
        int lb = cf.find_cliques();
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
            assert(g.nodes.size == cf.cliques[0].size());
        if (lb >= ub) {
            if (opt.learning > options::NO_LEARNING)
                return explain();
            else
                return INVALID_CLAUSE;
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
                    std::cout << "Partition " << v << " = "
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
    const std::vector<std::vector<Var>>& vars, const options& opt)
{
    auto cons = new gc_constraint(s, g, vars, opt);
    s.addConstraint(cons);
    return cons;
}

} // namespace gc
