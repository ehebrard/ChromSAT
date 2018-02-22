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

    enum learning_level { NO_LEARNING, NAIVE_EXPL };
    learning_level learning{NAIVE_EXPL};

    bool xvars;
		
		int polarity;
};

options parse(int argc, char* argv[]);

} // namespace gc

#endif
