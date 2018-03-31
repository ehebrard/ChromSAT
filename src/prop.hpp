#ifndef __GC_PROP_HPP
#define __GC_PROP_HPP

#include "graph.hpp"
#include "minicsp/core/solver.hpp"
#include "options.hpp"
#include "statistics.hpp"

namespace gc
{

struct cons_base {
    int ub;
    int bestlb{0};
    clique_finder cf;
		mycielskan_subgraph_finder mf;

    explicit cons_base(graph& g)
        : cf(g), mf(g, cf)
    {
    }
};

cons_base* post_gc_constraint(minicsp::Solver& s, graph& g,
    const std::vector<std::vector<minicsp::Var>>& vars, const options& opt, statistics& stat);

} // namespace gc

#endif
