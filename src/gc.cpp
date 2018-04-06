#include <iostream>

#include "dimacs.hpp"
#include "graph.hpp"
#include "options.hpp"
#include "prop.hpp"
#include "rewriter.hpp"
#include "statistics.hpp"
#include "utils.hpp"

#include <minicsp/core/cons.hpp>
#include <minicsp/core/solver.hpp>
#include <minicsp/core/utils.hpp>

struct gc_model {
    gc::graph& g;
    minicsp::Solver& s;
    const gc::options& options;
    gc::statistics& statistics;
    std::vector<std::vector<minicsp::Var>> vars;
    gc::cons_base* cons;
    std::vector<minicsp::cspvar> xvars;
    gc::rewriter rewriter;

    std::vector<std::vector<minicsp::Var>> create_vars()
    {
        std::vector<std::vector<minicsp::Var>> vars;
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
        return vars;
    }

    gc_model(gc::graph& g, minicsp::Solver& s, const gc::options& options,
        gc::statistics& statistics, std::pair<int, int> bounds)
        : g(g)
        , s(s)
        , options(options)
        , statistics(statistics)
        , vars(create_vars())
        , cons(gc::post_gc_constraint(s, g, vars, options, statistics))
        , rewriter(s, g, *cons, vars, xvars)
    {
        auto [lb, ub] = bounds;
        cons->bestlb = std::max(lb, cons->bestlb);
        cons->ub = std::min(ub, cons->ub);

        if (options.xvars) {
            xvars = s.newCSPVarArray(g.capacity(), 0, cons->ub - 2);
            for (size_t i = 0; i != xvars.size(); ++i) {
                for (size_t j = i + 1; j != xvars.size(); ++j) {
                    if (g.matrix[i].fast_contain(j))
                        minicsp::post_neq(s, xvars[i], xvars[j], 0);
                    else
                        minicsp::post_eq_re(
                            s, xvars[i], xvars[j], 0, minicsp::Lit(vars[i][j]));
                }
            }

            // rewrite clauses to not use x literals
            s.use_clause_callback([this](vec<minicsp::Lit>& clause, int btlvl) {
                return rewriter.rewrite(clause, btlvl);
            });
        }

        switch (options.branching) {
        case gc::options::VSIDS:
            s.varbranch = minicsp::VAR_VSIDS;
            break;
        case gc::options::BRELAZ:
            if (!options.xvars) {
                std::cout << "Cannot use Brelaz ordering without xvars\n";
                exit(1);
            }
            s.varbranch = minicsp::VAR_DOM;
            s.valbranch = minicsp::VAL_LEX;
            break;
        }
    }

    void solve()
    {
        using minicsp::l_False;
        using minicsp::l_True;
        using minicsp::l_Undef;

        minicsp::lbool sat{l_True};
        while (sat != l_False) {
            sat = s.solveBudget();
            if (sat == l_True) {
                std::cout << "c new UB " << g.nodes.size()
                          << " time = " << minicsp::cpuTime()
                          << " conflicts = " << s.conflicts << std::endl;
                assert(g.nodes.size() < static_cast<size_t>(cons->ub));
                cons->ub = g.nodes.size();
                if (options.xvars) {
                    for (auto v : xvars)
                        v.setmax(s, cons->ub - 2, minicsp::NO_REASON);
                }
            } else if (sat == l_Undef) {
                std::cout << "*** INTERRUPTED ***\n";
                return;
            } else {
                cons->bestlb = cons->ub;
                std::cout << "c new lower bound " << cons->ub
                          << " time = " << minicsp::cpuTime()
                          << " conflicts = " << s.conflicts << std::endl;
                std::cout << "UNSAT\n";
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
        statistics.describe(std::cout);
    }
};

std::pair<int, int> preprocess(gc::graph& g)
{
    gc::clique_finder cf{g};
    int lb{cf.find_cliques(g.nodes)};
    auto sol{gc::brelaz_color(g)};
    for (auto u : g.nodes)
        for (auto v : g.matrix[u])
            assert(sol[u] != sol[v]);

    int ub{*max_element(begin(sol), end(sol)) + 1};
    std::cout << "c new UB " << ub << " time = " << minicsp::cpuTime()
              << " conflicts = 0" << std::endl;
    return std::pair<int, int>{lb, ub};
}

int main(int argc, char* argv[])
{
    auto options = gc::parse(argc, argv);
    options.describe(std::cout);

    gc::statistics statistics;

    gc::graph g;
    dimacs::read_graph(options.instance_file.c_str(),
        [&](int nv, int) { g = gc::graph{nv}; },
        [&](int u, int v) {
            if (u != v)
                g.add_edge(u - 1, v - 1);
        },
        [&](int, gc::weight) {});
    g.describe(std::cout);

    minicsp::Solver s;
    setup_signal_handlers(&s);
    s.trace = options.trace;
    s.polarity_mode = options.polarity;

    auto [lb, ub] = preprocess(g);

    if (options.learning == gc::options::NO_LEARNING)
        s.learning = false;

    gc_model model(g, s, options, statistics, std::pair<int, int>(lb, ub));
    model.solve();
    model.print_stats();
}
