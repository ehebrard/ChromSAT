#ifndef GC_BRANCHER_HPP
#define GC_BRANCHER_HPP

#include "graph.hpp"
#include "prop.hpp"

namespace gc
{

struct Brancher {
    minicsp::Solver& s;
    graph& g;
    const std::vector<std::vector<minicsp::Var>>& evars;
    const std::vector<minicsp::cspvar>& xvars;
    cons_base& constraint;

    int64_t numdecisions{0}, numchoices{0};

    Brancher(minicsp::Solver& s, graph& g,
        const std::vector<std::vector<minicsp::Var>>& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint)
        : s(s)
        , g(g)
        , evars(evars)
        , xvars(xvars)
        , constraint(constraint)
    {
    }
    virtual ~Brancher() {}

    void use()
    {
        s.varbranch = minicsp::VAR_USER;
        s.user_brancher = [this](std::vector<minicsp::Lit>& cand) {
            constraint.sync_graph();

            select_candidates(cand);
            ++numdecisions;
            numchoices += cand.size();
        };
    }

    virtual void select_candidates(std::vector<minicsp::Lit>& cand) = 0;
};

struct BrelazBrancher : public Brancher {
    using Brancher::Brancher;
    std::vector<int> mindom;
    bitset util_set;

    BrelazBrancher(minicsp::Solver& s, graph& g,
        const std::vector<std::vector<minicsp::Var>>& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint)
        : Brancher(s, g, evars, xvars, constraint)
        , util_set(0, g.capacity() - 1, bitset::empt)
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        int mind{-1};
        mindom.clear();
        for (auto v : g.nodes) {
            auto x = xvars[v];
            auto xd = x.domsize(s);
            if (xd == 1)
                continue;
            if (mindom.empty() || xd < mind) {
                mindom.clear();
                mind = x.domsize(s);
                mindom.push_back(v);
            } else if (xd == mind)
                mindom.push_back(v);
        }

        if (mindom.empty())
            return;

        int tiedv = mindom.size();
        int maxdv{-1};
        int maxd{-1};
        if (mindom.size() > 1) {
            for (auto v : mindom) {
                util_set.copy(g.matrix[v]);
                util_set.intersect_with(g.nodeset);
                int deg = util_set.size();
                if (deg > maxd) {
                    maxd = deg;
                    maxdv = v;
                }
            }
        } else
            maxdv = mindom[0];
        auto x = xvars[maxdv];
        // std::cout << "Branching on vertex " << maxdv
        //           << " domsize = " << x.domsize(s) << " degree = " << maxd
        //           << " domsize ties = " << tiedv << "\n";
        cand.clear();
        cand.push_back(x.e_eq(s, x.min(s)));
    }
};

} // namespace gc

#endif
