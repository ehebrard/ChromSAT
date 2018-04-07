#include <iostream>

#include "brancher.hpp"
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
    minicsp::Solver s;
    const gc::options& options;
    gc::statistics& statistics;
    std::vector<std::vector<minicsp::Var>> vars;
    gc::cons_base* cons;
    std::vector<minicsp::cspvar> xvars;
    gc::rewriter rewriter;

    std::unique_ptr<gc::Brancher> brancher;
		
		// std::vector<edge> var_map;

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

    gc_model(gc::graph& g, const gc::options& options,
        gc::statistics& statistics, std::pair<int, int> bounds)
        : g(g)
        , options(options)
        , statistics(statistics)
        , vars(create_vars())
        , cons(gc::post_gc_constraint(s, g, vars, options, statistics))
        , rewriter(s, g, *cons, vars, xvars)
    {
        setup_signal_handlers(&s);
        s.trace = options.trace;
        s.polarity_mode = options.polarity;

        if (options.learning == gc::options::NO_LEARNING)
            s.learning = false;

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
            brancher = std::make_unique<gc::BrelazBrancher>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::PARTITION_SUM:
            brancher = std::make_unique<gc::PartitionBrancher<-1,-1,[](int x, int y) { return x+y; }>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::PARTITION_PRODUCT:
            brancher = std::make_unique<gc::PartitionBrancher<-1,-1,[](int x, int y) { return x*y; }>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::DEGREE_SUM:
            brancher = std::make_unique<gc::DegreeBrancher<-1,-1,[](int x, int y) { return x+y; }>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::DEGREE_PRODUCT:
            brancher = std::make_unique<gc::DegreeBrancher<-1,-1,[](int x, int y) { return x*y; }>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::DEGREE_UNION:
            brancher = std::make_unique<gc::DegreeUnionBrancher<-1,-1>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::PARTITION_SUM_DYN:
            brancher = std::make_unique<gc::PartitionBrancher<2,3,[](int x, int y) { return x+y; }>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::PARTITION_PRODUCT_DYN:
            brancher = std::make_unique<gc::PartitionBrancher<2,3,[](int x, int y) { return x*y; }>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::DEGREE_SUM_DYN:
            brancher = std::make_unique<gc::DegreeBrancher<2,3,[](int x, int y) { return x+y; }>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::DEGREE_PRODUCT_DYN:
            brancher = std::make_unique<gc::DegreeBrancher<2,3,[](int x, int y) { return x*y; }>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        case gc::options::DEGREE_UNION_DYN:
            brancher = std::make_unique<gc::DegreeUnionBrancher<2,3>>(
                s, g, vars, xvars, *cons);
            brancher->use();
            break;
        }
    }

    std::pair<int, int> solve()
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
                          << " conflicts = " << s.conflicts 
													<< " delta = " << statistics.get_bound_increase() << std::endl;
                assert(g.nodes.size() < static_cast<size_t>(cons->ub));
                cons->ub = g.nodes.size();
                if (options.xvars) {
                    for (auto v : xvars)
                        v.setmax(s, cons->ub - 2, minicsp::NO_REASON);
                }
            } else if (sat == l_Undef) {
                std::cout << "*** INTERRUPTED ***\n";
                break;
            } else {
                cons->bestlb = cons->ub;
                std::cout << "c new lower bound " << cons->ub
                          << " time = " << minicsp::cpuTime()
                          << " conflicts = " << s.conflicts 
													<< " delta = " << statistics.get_bound_increase() << std::endl;
                std::cout << "UNSAT\n";
            }
        }
        return std::make_pair(cons->bestlb, cons->ub);
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
		gc::mycielskan_subgraph_finder mf(g, cf, false);
    int lb{cf.find_cliques(g.nodes)};
		lb = mf.improve_cliques_larger_than(lb-1);
    auto sol{gc::brelaz_color(g)};
    for (auto u : g.nodes)
        for (auto v : g.matrix[u])
            assert(sol[u] != sol[v]);

    int ub{*max_element(begin(sol), end(sol)) + 1};
    std::cout << "c new UB " << ub << " time = " << minicsp::cpuTime()
              << " conflicts = 0" 
							<< " delta = 0" << std::endl;
    std::cout << "c new lower bound " << lb << " time = " << minicsp::cpuTime()
              << " conflicts = 0" 
							<< " delta = 0" << std::endl;

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

    auto [lb, ub] = preprocess(g);

    switch(options.strategy) {
		    case gc::options::BNB: {
		        gc_model model(g, options, statistics, std::make_pair(lb, ub));
		        model.solve();
		        model.print_stats();
		        break;
		    }
		    case gc::options::BOTTOMUP: {
		        for (int i = lb; i < ub; ++i) {
		            gc::graph gcopy{g};
		            gc_model model(
		                gcopy, options, statistics, std::make_pair(i, i + 1));
		            auto [ilb, iub] = model.solve();
		            if (iub == i) {
		                std::cout << "OPTIMUM " << ub << "\n";
		                break;
		            } else if (ilb != iub) {
		                std::cout << "best bounds [" << i << "," << ub << "\n";
		                std::cout << "INTERRUPTED\n";
		            } else {
		                std::cout << "c new lower bound " << ilb
		                          << " time = " << minicsp::cpuTime() 
															<< " conflicts = 0 delta = 0"<< "\n";
		            }
		        }
		    }
		}
}

