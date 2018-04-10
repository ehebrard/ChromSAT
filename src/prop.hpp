#ifndef __GC_PROP_HPP
#define __GC_PROP_HPP

#include "graph.hpp"
#include "minicsp/core/solver.hpp"
#include "options.hpp"

namespace gc
{

struct statistics;

struct cons_base {
    minicsp::Solver& s;
    graph& g;
    int ub;
    int bestlb{0};
    clique_finder cf;
    minicsp::backtrackable<int> lastlb;

    explicit cons_base(minicsp::Solver& s, graph& g)
        : s(s)
        , g(g)
        , cf(g)
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

cons_base* post_gc_constraint(minicsp::Solver& s, graph& g,
    const std::vector<std::vector<minicsp::Var>>& vars, const options& opt,
    statistics& stat);

} // namespace gc

#endif
