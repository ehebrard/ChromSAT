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

    enum learning_level {
        NO_LEARNING,
        NAIVE_POSITIVE,
        NAIVE
    };
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
        FULLMYCIELSKI
    };
    dual_policy boundalg;
};

options parse(int argc, char* argv[]);

} // namespace gc

#endif
