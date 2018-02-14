#ifndef __GC_PROP_HPP
#define __GC_PROP_HPP

#include "graph.hpp"
#include "options.hpp"
#include "minicsp/core/solver.hpp"

namespace gc
{

struct cons_base {
    int ub;
    int bestlb{0};
};

cons_base* post_gc_constraint(minicsp::Solver& s, graph& g,
    const std::vector<std::vector<minicsp::Var>>& vars, const options& opt);

} // namespace gc

#endif
