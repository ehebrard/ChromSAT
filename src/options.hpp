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
    std::string solution_file;

    // minicsp options
    bool trace{false};

    enum learning_level {
        NO_LEARNING,
        NAIVE_POSITIVE,
        CHOOSE_POSITIVE,
        DYNAMIC_REPS
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

    enum ordering_low_degree_sorting {
        NO_LD_ORDERING,
        PREPROCESSING_ORDERING,
        DEGREE_ORDERING
    };
    ordering_low_degree_sorting ordering_low_degree;

    enum dual_policy {
        CLIQUES,
        GREEDYMYCIELSKI,
        MAXMYCIELSKI,
        FULLMYCIELSKI,
    };
    dual_policy boundalg;

    bool prune;

    bool adaptive;

    enum branching_heuristic {
        VSIDS,
        BRELAZ,
        PARTITION_PRODUCT,
        PARTITION_SUM,
        DEGREE_PRODUCT,
        DEGREE_SUM,
        DEGREE_UNION,
        PARTITION_PRODUCT_DYN,
        PARTITION_SUM_DYN,
        DEGREE_PRODUCT_DYN,
        DEGREE_SUM_DYN,
        DEGREE_UNION_DYN,
        VSIDS_PHASED,
        VSIDS_GUIDED,
        VSIDS_CLIQUE,
        VSIDS_COLORS_POSITIVE
        // ,BRELAZ_GUIDED
    };
    branching_heuristic branching;

    bool branching_low_degree;

    int cliquelimit;

    enum solution_strategy {
        BNB,
        BOTTOMUP,
        TOPDOWN,
        BOUNDS,
        CLEVER,
        COLOR,
        TEST,
        IDSATUR
    };
    solution_strategy strategy;

    enum preprocessing_types { NO_PREPROCESSING, LOW_DEGREE, FULL };
    preprocessing_types preprocessing;

    bool indset_constraints;

    enum fillin_type { FILLIN_NONE = 0, FILLIN_PROPAGATE, FILLIN_DECOMPOSE };
    int fillin;

    bool dominance;

    bool indset_lb;

    std::string format;

    int verbosity;

    bool checksolution;
    bool printsolution;

    int sdsaturiter;

    int ddsaturiter;

    std::string convert;

    bool equalities;

    bool dsatur;

    int memlimit;

    int myciellimit;

    bool ubfocus;

    bool triangle_up;

    int samplebase;

    int probewidth;
		
    enum core_type {
        ALL,
        WITNESS,
        LOWER
    };
    core_type core;
};

options parse(int argc, char* argv[]);

} // namespace gc

#endif
