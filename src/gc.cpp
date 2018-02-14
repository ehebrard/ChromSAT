#include <iostream>

#include "graph.hpp"
#include "dimacs.hpp"
#include "prop.hpp"
#include "options.hpp"

#include <minicsp/core/solver.hpp>
#include <minicsp/core/utils.hpp>

struct gc_model
{
    gc::graph& g;
    minicsp::Solver& s;
    const gc::options& options;
    std::vector<std::vector<minicsp::Var>> vars;
    gc::cons_base *cons;

    gc_model(gc::graph& g, minicsp::Solver& s, const gc::options& options)
        : g(g)
        , s(s)
        , options(options)
    {
        vars.resize(g.capacity());
        for (size_t i = 0; i != vars.size(); ++i) {
            vars[i].resize(g.capacity());
            vars[i][i] = minicsp::var_Undef;
            for (size_t j = 0; j != i; ++j)
                vars[i][j] = vars[j][i];
            for (size_t j = i + 1; j != vars.size(); ++j)
                if (g.matrix[i].fast_contain(j))
                    vars[i][j] = minicsp::var_Undef;
                else {
                    vars[i][j] = s.newVar();
                    if (options.trace) {
                        using namespace std::string_literals;
                        using std::to_string;
                        auto n = "e"s + to_string(i) + "-"s + to_string(j);
                        s.setVarName(vars[i][j], n);
                    }
                }
        }

        cons = gc::post_gc_constraint(s, g, vars, options);
    }

    void solve()
    {
        using minicsp::l_False;
        using minicsp::l_True;
        using minicsp::l_Undef;

        minicsp::lbool sat{l_True};
        while (sat != l_False) {
            s.solveBudget();
            if (sat == l_True) {
                std::cout << "c new UB " << g.nodes.size
                          << " time = " << minicsp::cpuTime()
                          << " conflicts = " << s.conflicts << std::endl;
                assert(g.nodes.size < static_cast<size_t>(cons->ub));
                cons->ub = g.nodes.size;
            } else if (sat == l_Undef) {
                std::cout << "*** INTERRUPTED ***\n";
                return;
            }
        }
    }

    void print_stats()
    {
        if (cons->bestlb >= cons->ub)
            std::cout << "OPTIMUM " << cons->ub << "\n";
        else
            std::cout << "Best bounds [" << cons->bestlb << ", " << cons->ub
                      << "]\n";
        minicsp::printStats(s);
    }
};

int main(int argc, char *argv[])
{
    auto options = gc::parse(argc, argv);
    options.describe(std::cout);

    gc::graph g;
    dimacs::read_graph(options.instance_file.c_str(),
        [&](int nv, int) { g = gc::graph{nv}; },
        [&](int u, int v) { g.add_edge(u - 1, v - 1); },
        [&](int, gc::weight) {});
    g.describe(std::cout);

    minicsp::Solver s;
    setup_signal_handlers(&s);
    s.trace = options.trace;
    if (options.learning == gc::options::NO_LEARNING)
        s.learning = false;

    gc_model model(g, s, options);
    model.solve();
    model.print_stats();
}
