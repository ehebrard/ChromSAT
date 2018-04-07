#ifndef __GC_OPTIONS_HPP
#define __GC_OPTIONS_HPP

#include <iostream>
#include <string>

namespace gc
{
struct options {
    // outputs a nice description of all options
    void describe(std::ostream&);

    // the actual options
    std::string cmdline; // for reference
    std::string instance_file;

    // minicsp options
    bool trace{false};

    enum learning_level { NO_LEARNING, NAIVE_POSITIVE, NAIVE, MYC_POSITIVE };
    learning_level learning{NAIVE_POSITIVE};

    bool xvars;

    int polarity;

    enum ordering_heuristic {
        NONE,
        DEGENERACY,
        INVERSE_DEGENERACY,
        PARTITION,
        DYNAMIC_DEGENERACY
    };
    ordering_heuristic ordering;

    enum dual_policy {
        CLIQUES,
        GREEDYMYCIELSKI,
        MAXMYCIELSKI,
        FULLMYCIELSKI,
    };
    dual_policy boundalg;

    bool prune;

    bool adaptive;

    enum branching_heuristic { VSIDS, BRELAZ, PARTITION_PRODUCT, PARTITION_SUM, DEGREE_PRODUCT, DEGREE_SUM, DEGREE_UNION, PARTITION_PRODUCT_DYN, PARTITION_SUM_DYN, DEGREE_PRODUCT_DYN, DEGREE_SUM_DYN, DEGREE_UNION_DYN };
    branching_heuristic branching;
		
		int cliquelimit;
};

options parse(int argc, char* argv[]);

} // namespace gc

#endif
