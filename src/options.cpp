#include "options.hpp"

#include <tclap/CmdLine.h>

#include <numeric>

namespace gc
{

struct argbase {
    virtual ~argbase() {}
    virtual void assign() = 0;
};

template <typename Opt, typename ClapArg, typename E = void>
struct arg : public argbase {
    ClapArg carg;
    Opt& opt;

    template <typename... T>
    arg(TCLAP::CmdLine& cmd, Opt& opt, T&&... args)
        : carg(std::forward<T>(args)...)
        , opt(opt)
    {
        cmd.add(carg);
    }

    virtual void assign() override { opt = carg.getValue(); }
};

template <typename Opt, typename ClapArg>
struct arg<Opt, ClapArg, typename std::enable_if<std::is_enum<Opt>{}>::type>
    : public argbase {
    ClapArg carg;
    Opt& opt;

    template <typename... T>
    arg(TCLAP::CmdLine& cmd, Opt& opt, T&&... args)
        : carg(std::forward<T>(args)...)
        , opt(opt)
    {
        cmd.add(carg);
    }

    virtual void assign() override
    {
        opt = static_cast<typename std::remove_reference<Opt>::type>(
            carg.getValue());
    }
};

struct cmdline {
    TCLAP::CmdLine cmd;
    std::vector<std::unique_ptr<argbase>> args;

    cmdline(const std::string& message, const char delimiter = ' ',
        const std::string& version = "none", bool helpAndVersion = true)
        : cmd(message, delimiter, version, helpAndVersion)
    {
    }

    template <typename ClapArg, typename Opt, typename... T>
    void add(Opt& opt, T&&... clapargs)
    {
        args.emplace_back(std::move(std::make_unique<arg<Opt, ClapArg>>(
            cmd, opt, std::forward<T>(clapargs)...)));
    }

    void parse(int argc, char* argv[])
    {
        cmd.parse(argc, argv);
        for (auto& arg : args)
            arg->assign();
    }
};

options parse(int argc, char* argv[])
{
    using namespace TCLAP;
    using namespace std::string_literals;
    cmdline cmd("GraphColoring", ' ');

    options opt;
    opt.cmdline = std::accumulate(argv, argv + argc, ""s,
        [&](std::string acc, const char* arg) { return acc + " " + arg; });

    cmd.add<UnlabeledValueArg<std::string>>(opt.instance_file, "file",
        "instance file", true, "data/DIMACS_cliques/brock200_1.clq", "string");
    cmd.add<UnlabeledValueArg<std::string>>(
        opt.solution_file, "solfile", "solution file", false, "", "string");
    cmd.add<SwitchArg>(opt.trace, "", "trace", "enable minicsp tracing", false);
    cmd.add<ValueArg<int>>(opt.learning, "", "learning",
        "CDCLeaning & explanation level [0-1]", false, 1, "int");
    cmd.add<SwitchArg>(
        opt.xvars, "", "xvars", "add x (color) variables to the model", false);
    cmd.add<ValueArg<int>>(
        opt.polarity, "", "polarity", "polarity policy", false, 0, "int");
    cmd.add<ValueArg<int>>(opt.ordering, "", "ordering",
        "clique finding heuristic [0-4]", false, 3, "int");
    cmd.add<ValueArg<int>>(opt.ordering_low_degree, "", "ord-low-degree",
        "Use low degree information to improve clique ordering", false, 0, "int");
    cmd.add<ValueArg<int>>(opt.boundalg, "", "bound",
        "lower bound algorithm [0-3]", false, 2, "int");
    cmd.add<SwitchArg>(opt.prune, "", "prune", "enable pruning", false);
    cmd.add<SwitchArg>(opt.adaptive, "", "adaptive",
        "Switch between CLIQUES and declared bound policy dynamically", true);
    cmd.add<ValueArg<int>>(opt.branching, "", "branching",
        "Variable branching heuristic [0-14]", false, 1, "int");
    cmd.add<SwitchArg>(opt.branching_low_degree, "", "branch-low-degree",
        "Use low degree information to improve branching", false);
    cmd.add<ValueArg<int>>(opt.cliquelimit, "", "cliquelimit",
        "Maximum number of cliques in the lower bound algorithm during "
        "preprocessing",
        false, 1000, "int");
    cmd.add<ValueArg<int>>(opt.strategy, "", "strategy",
        "Solution strategy [0=BNB-1=bottom-up-2=top-down-3=preprocessing "
        "only-4=top-down as in (Verma et al.)-5=color-6=test]",
        false, 0, "int");
    cmd.add<ValueArg<int>>(opt.preprocessing, "", "preprocessing",
        "Level of preprocessing [0=none-1=low-degree-2=low-degree (sparse)]", false, 1, "int");
    cmd.add<SwitchArg>(opt.dominance, "", "dominance",
        "enable neighborhood-based dominance", false);
    cmd.add<SwitchArg>(opt.indset_constraints, "", "indset",
        "reduce by independent set constraints", false);
    cmd.add<ValueArg<int>>(opt.fillin, "", "fillin",
        "compute fill-in of graph (0:no,1:prop,2:decomp", false, 0, "int");
    cmd.add<SwitchArg>(opt.indset_lb, "", "indset-lb",
        "compute IS-based lower bound", false);
    cmd.add<ValueArg<std::string>>(opt.format, "", "format",
        "File format", false, "dimacs", "string");
    cmd.add<ValueArg<int>>(
        opt.verbosity, "", "verbosity", "Verbosity level", false, 0, "int");
    cmd.add<SwitchArg>(
        opt.checksolution, "", "checksolution", "checks the coloring", false);
    cmd.add<SwitchArg>(
        opt.printsolution, "", "printsolution", "prints the coloring", false);
    cmd.add<ValueArg<int>>(opt.sdsaturiter, "", "sdsaturiter",
        "# of sparse dsatur iterations", false, 1, "int");
    cmd.add<ValueArg<int>>(opt.ddsaturiter, "", "ddsaturiter",
        "# of dense dsatur iterations", false, 1, "int");
    cmd.add<ValueArg<std::string>>(opt.convert, "", "convert",
        "output in <file> using dimacs format", false, "", "string");
    cmd.add<SwitchArg>(opt.equalities, "", "equalities",
        "add implied number of equalitiy constraints", false);
    cmd.add<SwitchArg>(
        opt.dsatur, "", "dsatur", "use dsatur instead of minicsp", false);
    cmd.add<SwitchArg>(opt.ubfocus, "", "ubfocus",
        "do not try hard to get lb's [default=false]", false);
    cmd.add<ValueArg<int>>(opt.memlimit, "", "memlimit",
        "does not compute dense graph above this limit and downgrade reasoning "
        "on large instances (default: no limit)",
        false, -1, "int");
    cmd.add<ValueArg<int>>(opt.myciellimit, "", "myciellimit",
        "does not use myciel on graphs larger than the limit (default: no limit)",
        false, -1, "int");
    cmd.add<SwitchArg>(opt.triangle_up, "", "triangleup",
        "do partition-aware unit propagation", false);
    cmd.add<ValueArg<int>>(opt.samplebase, "", "samplebase",
        "size of the sampling base for max clique (default = 512)", false, 512,
        "int");
    cmd.add<ValueArg<int>>(opt.probewidth, "", "probewidth",
        "size of the max probe width (default = 64)", false, 64, "int");
    cmd.add<ValueArg<int>>(opt.core, "", "core",
        "Core type []",
        false, 0, "int");

    cmd.parse(argc, argv);
    return opt;
}

void options::describe(std::ostream& os)
{
    os << "[options] GC configuration\n";
    os << "[options] cmdline = " << cmdline << "\n";
    os << "[options] Instance file = " << instance_file << "\n";
    os << "[options] Clause learning = " << learning << "\n";
    os << "[options] Polarity policy = " << polarity << "\n";
    os << "[options] Clique ordering = " << ordering << "\n";
    os << "[options]  ... low degree = " << ordering_low_degree << "\n";
    os << "[options] Color variables = " << (xvars ? "present" : "absent") << "\n";
    os << "[options] Bound policy    = " << boundalg << "\n";
    os << "[options] Adaptive bounds = " << adaptive << "\n";
    os << "[options] Preprocessing   = " << (preprocessing == 0 ? "none" : preprocessing == 1 ?  "dense" : "sparse") << "\n";
    os << "[options] IS constraints  = " << (int)indset_constraints << "\n";
    os << "[options] fillin          = "
       << (fillin == options::FILLIN_NONE
                  ? "no"
                  : (fillin == options::FILLIN_PROPAGATE ? "yes (propagate)"
                                                         : "yes (decompose)"))
       << "\n";
    os << "[options] Triangle UP     = " << triangle_up << "\n";
    os << "[options] Branching strat = ";
    switch (branching) {
    case gc::options::VSIDS:
        os << "VSIDS\n";
        break;
    case gc::options::BRELAZ:
        os << "Brelaz\n";
        break;
    case gc::options::PARTITION_PRODUCT:
        os << "Largest edge partition size (x)\n";
        break;
    case gc::options::PARTITION_SUM:
        os << "Largest edge partition size (+)\n";
        break;
    case gc::options::DEGREE_PRODUCT:
        os << "Largest edge degree (x)\n";
        break;
    case gc::options::DEGREE_SUM:
        os << "Largest edge degree (+)\n";
        break;
    case gc::options::DEGREE_UNION:
        os << "Largest edge neighborhood\n";
        break;
    case gc::options::PARTITION_PRODUCT_DYN:
        os << "Largest edge partition size (x) (select)\n";
        break;
    case gc::options::PARTITION_SUM_DYN:
        os << "Largest edge partition size (+) (select)\n";
        break;
    case gc::options::DEGREE_PRODUCT_DYN:
        os << "Largest edge degree (x) (select)\n";
        break;
    case gc::options::DEGREE_SUM_DYN:
        os << "Largest edge degree (+) (select)\n";
        break;
    case gc::options::DEGREE_UNION_DYN:
        os << "Largest edge neighborhood (select)\n";
        break;
    case gc::options::VSIDS_PHASED:
        os << "VSIDS with partition phase heuristic\n";
        break;
    case gc::options::VSIDS_GUIDED:
        os << "VSIDS with solution phase saving\n";
        break;
    case gc::options::VSIDS_CLIQUE:
        os << "VSIDS restricted to merging with large cliques\n";
        break;
    case gc::options::VSIDS_COLORS_POSITIVE:
        os << "VSIDS restricted to assignments to color variables\n";
        break;
    }
    os << "[options]  ... low degree = " << branching_low_degree << "\n";
    os << "[options] Strategy        = ";
    switch (strategy) {
    case BNB:
        os << "branch and bound";
        break;
    case BOTTOMUP:
        os << "bottom-up";
        break;
    case TOPDOWN:
        os << "top-down";
        break;
    case BOUNDS:
        os << "bounds";
        break;
    case CLEVER:
        os << "(Verma et al.)";
        break;
    case COLOR:
        os << "heuristic coloring";
        break;
    case TEST:
        os << "test";
        break;
    case IDSATUR:
        os << "iterated dsatur";
        break;
    }
    os << "\n";
    os << std::endl;
}

} // namespace gc
