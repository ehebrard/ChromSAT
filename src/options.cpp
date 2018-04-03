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
    cmd.add<SwitchArg>(opt.trace, "", "trace", "enable minicsp tracing", false);
    cmd.add<ValueArg<int>>(opt.learning, "", "learning",
        "CDCLeaning & explanation level [0-1]", false, 1, "int");
    cmd.add<SwitchArg>(
        opt.xvars, "", "xvars", "add x (color) variables to the model", false);
    cmd.add<ValueArg<int>>(
        opt.polarity, "", "polarity", "polarity policy", false, 0, "int");
    cmd.add<ValueArg<int>>(opt.ordering, "", "ordering",
        "clique finding heuristic [0-4]", false, 3, "int");
    cmd.add<ValueArg<int>>(opt.boundalg, "", "bound",
        "lower bound algorithm [0-3]", false, 0, "int");
    cmd.add<SwitchArg>(opt.prune, "", "prune", "enable pruning", false);
    cmd.add<SwitchArg>(opt.adaptive, "", "adaptive",
        "Switch between CLIQUES and declared bound policy dynamically", false);
    cmd.add<ValueArg<int>>(opt.branching, "", "branching",
        "Variable branching heuristic [0-1]", false, 0, "int");

    cmd.parse(argc, argv);
    return opt;
}

void options::describe(std::ostream& os)
{
    os << "GC configuration\n";
    os << "cmdline = " << cmdline << "\n";
    os << "Instance file = " << instance_file << "\n";
    os << "Clause learning = " << learning << "\n";
    os << "Polarity policy = " << polarity << "\n";
    os << "Clique ordering = " << ordering << "\n";
    os << "Color variables = " << (xvars ? "present" : "absent") << "\n";
    os << "Bound policy    = " << boundalg << "\n";
    os << "Adaptive bounds = " << adaptive << "\n";
    os << std::endl;
}

} // namespace gc
