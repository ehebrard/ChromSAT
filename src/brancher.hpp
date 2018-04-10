#ifndef GC_BRANCHER_HPP
#define GC_BRANCHER_HPP

#include "graph.hpp"
#include "prop.hpp"
#include "options.hpp"

namespace gc
{

struct Brancher {
    minicsp::Solver& s;
    graph& g;
    const std::vector<std::vector<minicsp::Var>>& evars;
    const std::vector<minicsp::cspvar>& xvars;
    cons_base& constraint;
    const options& opt;

    int64_t numdecisions{0}, numchoices{0};

    Brancher(minicsp::Solver& s, graph& g,
        const std::vector<std::vector<minicsp::Var>>& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : s(s)
        , g(g)
        , evars(evars)
        , xvars(xvars)
        , constraint(constraint)
        , opt(opt)
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
    // a maximal clique
    std::vector<int> clique;
    bitset util_set;

    // vertices that are ignored because they have low degree
    std::vector<int> low_degree;

    BrelazBrancher(minicsp::Solver& s, graph& g,
        const std::vector<std::vector<minicsp::Var>>& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : Brancher(s, g, evars, xvars, constraint, opt)
        , util_set(0, g.capacity() - 1, bitset::empt)
    {
        auto& cf = constraint.cf;
        auto maxidx = std::distance(begin(cf.clique_sz),
            std::max_element(
                begin(cf.clique_sz), begin(cf.clique_sz) + cf.num_cliques));
        auto& clq = cf.cliques[maxidx];
        for (auto v : clq)
            clique.push_back(v);
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        for (auto v : clique) {
            auto x = xvars[v];
            if (x.domsize(s) != 1) {
                cand.clear();
                cand.push_back(x.e_eq(s, x.min(s)));
                return;
            }
        }

        int mind{-1};
        mindom.clear();
        low_degree.clear();
        for (auto v : g.nodes) {
            auto x = xvars[v];
            auto xd = x.domsize(s);
            if (xd == 1)
                continue;
            if (opt.branching_low_degree) {
                util_set.copy(g.matrix[v]);
                util_set.intersect_with(g.nodeset);
                if (util_set.size() < static_cast<size_t>(*constraint.lastlb)) {
                    low_degree.push_back(v);
                    continue;
                }
            }
            if (mindom.empty() || xd < mind) {
                mindom.clear();
                mind = x.domsize(s);
                mindom.push_back(v);
            } else if (xd == mind)
                mindom.push_back(v);
        }

        if (mindom.empty()) {
            if (!low_degree.empty()) {
                auto x = xvars[low_degree.back()];
                cand.clear();
                cand.push_back(x.e_eq(s, x.min(s)));
                return;
            } else
                return;
        }

        // int tiedv = mindom.size();
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

template<int N, int D>
struct EdgeBrancher : public Brancher {
    using Brancher::Brancher;
    std::vector<edge> e_cand;
    std::vector<int> nodes;
    bitset neighbors_u;
    bitset neighbors_v;
    bitset counter;
    int max_tied;

    EdgeBrancher(minicsp::Solver& s, graph& g,
        const std::vector<std::vector<minicsp::Var>>& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : Brancher(s, g, evars, xvars, constraint, opt)
        , neighbors_u(0, g.capacity() - 1, bitset::empt)
        , neighbors_v(0, g.capacity() - 1, bitset::empt)
        , counter(0, g.capacity() - 1, bitset::empt)
        , max_tied(g.capacity() * g.capacity())
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand) = 0;

    void select_branch(std::vector<minicsp::Lit>& cand)
    {
        for (auto e : e_cand) {
            if (N < 0)
                cand.push_back(minicsp::Lit(evars[e.first][e.second]));
            else {
                auto u{e.first};
                auto v{e.second};

                counter.copy(g.matrix[u]);
                counter.intersect_with(g.matrix[v]);
                counter.intersect_with(g.nodeset);

                auto inter_size = counter.size();

                counter.copy(g.matrix[u]);
                counter.union_with(g.matrix[v]);
                counter.intersect_with(g.nodeset);

                auto union_size = counter.size();

                if (inter_size * D > union_size * N) {
                    cand.push_back(minicsp::Lit(evars[u][v]));
                } else {
                    cand.push_back(~minicsp::Lit(evars[u][v]));
                }
            }
        }
    }
};

template <int N, int D, typename Op>
struct PartitionBrancher : public EdgeBrancher<N, D> {
    Op op;

    PartitionBrancher(minicsp::Solver& s, graph& g,
        const std::vector<std::vector<minicsp::Var>>& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt, Op op)
        : EdgeBrancher<N, D>(s, g, evars, xvars, constraint, opt)
        , op(op)
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        int best_crit{0};

        this->nodes.clear();
        std::copy(begin(this->g.nodeset), end(this->g.nodeset),
            back_inserter(this->nodes));

        this->e_cand.clear();
        for (auto u : this->nodes) {
            auto u_part_size{this->g.partition[u].size()};
            for (auto v : this->nodes) {
                if (u == v || this->g.matrix[u].fast_contain(v))
                    continue;

                auto criterion{op(u_part_size, this->g.partition[v].size())};
                if (criterion >= best_crit) {
                    if (criterion == best_crit)
                        this->e_cand.clear();
                    else
                        best_crit = criterion;
                    if (this->e_cand.size()
                        < static_cast<size_t>(this->max_tied))
                        this->e_cand.push_back(edge{u, v});
                }
            }
        }

        this->select_branch(cand);
    }
};

template <int N, int D, typename Op>
std::unique_ptr<PartitionBrancher<N, D, Op>> make_partition_brancher(
    minicsp::Solver& s, graph& g,
    const std::vector<std::vector<minicsp::Var>>& evars,
    const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
    const options& opt, Op op)
{
    return std::make_unique<PartitionBrancher<N, D, Op>>(
        s, g, evars, xvars, constraint, opt, op);
}

template <int N, int D, typename Op>
struct DegreeBrancher : public EdgeBrancher<N, D> {
    Op op;

    DegreeBrancher(minicsp::Solver& s, graph& g,
        const std::vector<std::vector<minicsp::Var>>& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt, Op op)
        : EdgeBrancher<N, D>(s, g, evars, xvars, constraint, opt)
        , op(op)
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        int best_crit{0};

        this->nodes.clear();
        std::copy(begin(this->g.nodeset), end(this->g.nodeset),
            back_inserter(this->nodes));

        this->e_cand.clear();
        for (auto u : this->nodes) {
            this->neighbors_u.copy(this->g.matrix[u]);
            this->neighbors_u.intersect_with(this->g.nodeset);
            auto u_degree{this->neighbors_u.size()};
            for (auto v : this->nodes) {
                if (u == v || this->g.matrix[u].fast_contain(v))
                    continue;

                this->neighbors_v.copy(this->g.matrix[v]);
                this->neighbors_v.intersect_with(this->g.nodeset);

                auto criterion{op(this->neighbors_v.size(), u_degree)};
                if (criterion >= best_crit) {
                    if (criterion == best_crit)
                        this->e_cand.clear();
                    else
                        best_crit = criterion;
                    if (this->e_cand.size()
                        < static_cast<size_t>(this->max_tied))
                        this->e_cand.push_back(edge{u, v});
                }
            }
        }

        this->select_branch(cand);
    }
};

template <int N, int D, typename Op>
std::unique_ptr<DegreeBrancher<N, D, Op>> make_degree_brancher(
    minicsp::Solver& s, graph& g,
    const std::vector<std::vector<minicsp::Var>>& evars,
    const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
    const options& opt, Op op)
{
    return std::make_unique<DegreeBrancher<N, D, Op>>(
        s, g, evars, xvars, constraint, opt, op);
}

template<int N, int D>
struct DegreeUnionBrancher : public EdgeBrancher<N, D> {
    // using EdgeBrancher::EdgeBrancher;

    DegreeUnionBrancher(minicsp::Solver& s, graph& g,
        const std::vector<std::vector<minicsp::Var>>& evars,
        const std::vector<minicsp::cspvar>& xvars, cons_base& constraint,
        const options& opt)
        : EdgeBrancher<N, D>(s, g, evars, xvars, constraint, opt)
    {
    }

    void select_candidates(std::vector<minicsp::Lit>& cand)
    {
        int best_crit{0};

        this->nodes.clear();
        std::copy(begin(this->g.nodeset), end(this->g.nodeset),
            back_inserter(this->nodes));

        this->e_cand.clear();
        for (auto u : this->nodes) {
            this->neighbors_u.copy(this->g.matrix[u]);
            for (auto v : this->nodes) {
                if (u == v || this->g.matrix[u].fast_contain(v))
                    continue;
                this->counter.copy(this->g.matrix[v]);
                this->counter.intersect_with(this->g.nodeset);
                int criterion = this->counter.size();
                if (criterion >= best_crit) {
                    if (criterion == best_crit)
                        this->e_cand.clear();
                    else
                        best_crit = criterion;
                    if (this->e_cand.size()
                        < static_cast<size_t>(this->max_tied))
                        this->e_cand.push_back(edge{u, v});
                }
            }
        }

        this->select_branch(cand);
    }
};

} // namespace gc

#endif
