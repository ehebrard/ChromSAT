#ifndef __GC_PROP_HPP
#define __GC_PROP_HPP

#include "graph.hpp"
#include "minicsp/core/solver.hpp"
#include "options.hpp"
#include "varmap.hpp"

#include <boost/optional.hpp>

namespace gc
{

using dense_graph = graph<bitset>;

struct statistics;

// constraint on a subset of the vertices of the graph that get
// added during preprocessing. We place a tighter upper bound on
// these sets of vertices
struct indset_constraint {

    template <class adjacency_struct>
    indset_constraint(adjacency_struct& scope, int s)
        : source(s)
        , vs(scope.min(), scope.max(), bitset::empt)
    {
        vs.copy(scope);
    }

    // the vertices in the local constraint
    bitset vs;
    // where it came from (i.e., which vertex's neighborhood is
    // it). for debugging
    int source;
};

struct cons_base {
    minicsp::Solver& s;
    dense_graph& g;
    int ub;
    int bestlb{0};
    clique_finder<bitset> cf; // for cliques
    clique_finder<bitset> ccf; // for clique covers
    minicsp::backtrackable<int> lastlb;

    explicit cons_base(minicsp::Solver& s, dense_graph& g)
        : s(s)
        , g(g)
        , cf(g)
        , ccf(g)
        , lastlb(s)
        , lastdlvl(s)
    {
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

protected:
    // something which is kept in sync by the solver. if it diverges
    // from what the graph thinks, it means we have
    // backtracked/restarted, so we should resync to that point
    minicsp::backtrackable<int> lastdlvl;
};

cons_base* post_gc_constraint(minicsp::Solver& s, dense_graph& g,
    boost::optional<std::vector<std::pair<int, int>>> fillin, const varmap& vars,
    const std::vector<indset_constraint>& isconses, const options& opt,
    statistics& stat);

} // namespace gc

#endif
