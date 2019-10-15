#ifndef __CG_ALGORITHM_HH
#define __CG_ALGORITHM_HH

#include <limits>
#include <vector>

#include "bi_graph.hpp"
#include "cliquer.hpp"
#include "dsatur.hpp"
#include "intstack.hpp"
#include "statistics.hpp"

namespace gc
{

template <class graph_struct> class coloring_algorithm;

template <class graph_struct> class heuristic
{

private:
    coloring_algorithm<graph_struct>& env;

    degeneracy_finder<graph_struct> df;
    std::vector<int> degeneracy;
    std::vector<int> dg_rank;

public:
    // dsatur brelaz;

    heuristic(coloring_algorithm<graph_struct>& env)
        : env(env)
        , df(env.G)
    {
        // N = env.G.nodes.capacity();
        degeneracy.resize(env.N);
        // max_clique.initialise(0, N - 1, gc::bitset::empt);
        // brelaz.use_recolor = false;
    }

    int degeneracy_coloring(std::vector<int>& coloring)
    {

        /*** compute degeneracy and k-cores, and initialise ub ***/
        df.degeneracy_ordering();
        int ub{df.degeneracy + 1};

        /*** compute the coloring corresponding to the degeneracy ordering
        * ***/
        // brelaz.use_recolor = false;
        auto bub{
            env.brelaz.greedy(env.G, rbegin(df.order), rend(df.order), ub)};

        if (ub > bub)
            ub = bub;

        degeneracy = df.degrees;
        coloring.resize(env.N, 0);
        dg_rank.resize(env.N);
        int r{0};
        for (auto d{begin(df.order)}; d < end(df.order); ++d) {
            auto v{*d};
            coloring[*d] = env.brelaz.color[*d];
            dg_rank[*d] = r++;
            if (d != begin(df.order))
                degeneracy[v] = std::max(degeneracy[*(d - 1)], degeneracy[v]);
        }

        return ub;
    }

    int dsatur_coloring(const int ub, std::vector<int>& coloring)
    {
        env.l.get_largest_clique();

        double tbefore = minicsp::cpuTime();
        auto global_criterion = [&](int better, int worse) {
            // degree [DYNAMIC]
            // if (G.matrix[better].size() > G.matrix[worse].size())
            // if (dg_rank[better] > dg_rank[worse])
            //             return true;

            // 1/ start with the max clique
            auto better_in_maxclique{env.l.max_clique.contain(better)};
            auto worse_in_maxclique{env.l.max_clique.contain(worse)};
            if (better_in_maxclique and !worse_in_maxclique)
                return true;
            if (worse_in_maxclique != better_in_maxclique)
                return false;

            if (degeneracy[better] > degeneracy[worse]
                or (degeneracy[better] == degeneracy[worse]
                       and env.G.matrix[better].size()
                           > env.G.matrix[worse].size()))
                return true;
            return false;
        };

        env.brelaz.clear();
        /*** compute a better coloring using dsatur with complex tie breaking
         * ***/
        auto nub{env.brelaz.brelaz_color_score(
            env.G, ub - 1, global_criterion, env.G.size(), 12345)};

        if (nub < ub) {

            coloring.resize(env.N, 0);

            for (auto v{0}; v < env.N; ++v) {
                coloring[v] = env.brelaz.color[v];
            }
        } else
            nub = ub;
        env.stats.notify_dsatur_time(minicsp::cpuTime() - tbefore);

        return nub;
    }
};

template <class graph_struct> class lower_bound
{

private:
    coloring_algorithm<graph_struct>& env;

    std::vector<int> largest_cliques;

    gc::cliquer<graph_struct> cq;
    gc::bi_graph B;

public:
    int clique_sz;
    gc::bitset max_clique;

    lower_bound(coloring_algorithm<graph_struct>& env)
        : env(env)
        , cq(env.G)
        , clique_sz{0}
        , max_clique(0, env.G.capacity() - 1, gc::bitset::empt)
    {
    }

    template <class random_it> int clique(random_it beg, random_it end)
    {
        double tbefore = minicsp::cpuTime();
        cq.clear();
        clique_sz = cq.find_cliques(beg, end, env.options.cliquelimit);

        env.stats.notify_nclique(cq.num_cliques);
        env.stats.notify_clique_time(minicsp::cpuTime() - tbefore);

        return clique_sz;
    }

    void get_largest_clique()
    {
        if (clique_sz == 0)
            return;

        largest_cliques.clear();

        for (auto i{0}; i < cq.cliques.size(); ++i) {
            if (cq.cliques[i].size() >= clique_sz) {
                largest_cliques.push_back(i);
            }
        }

        max_clique.clear();
        max_clique.union_with(cq.cliques[largest_cliques[0]]);
    }

    int matching()
    {
        double tbefore = minicsp::cpuTime();
        get_largest_clique();

        int matching_bound = clique_sz;

        if (largest_cliques.size() > 1) {

            for (auto i{0}; i < largest_cliques.size(); ++i) {
                auto c1{largest_cliques[i]};

                for (auto j{i + 1}; j < largest_cliques.size(); ++j) {
                    auto c2{largest_cliques[j]};

                    B.get_from_cliques(env.G, cq.cliques[c1], cq.cliques[c2]);
                    auto mm{B.hopcroftKarp()};

                    auto mm_bound = (cq.cliques[c1].size()
                        + cq.cliques[c2].size() - B.I - mm);

                    if (mm_bound > matching_bound) {
                        matching_bound = mm_bound;
                    }
                }
            }

            env.stats.notify_bound_delta(clique_sz, matching_bound);
        }
        env.stats.notify_matching_time(minicsp::cpuTime() - tbefore);

        return matching_bound;
    }
};

template <class graph_struct> class selector
{

private:
    coloring_algorithm<graph_struct>& env;

public:
    selector(coloring_algorithm<graph_struct>& env)
        : env(env)
    {
    }

    arc select()
    {

        // std::cout << env.l.max_clique << " " << env.l.clique_sz << "/" <<
        // env.l.max_clique.size() << std::endl;

        int c = 0;
        for (auto vp{begin(env.brelaz.order)};
             vp != begin(env.brelaz.order) + env.l.clique_sz; ++vp) {

            // std::cout << *vp << " " << c << std::endl;

            assert(env.brelaz.color[*vp] == c++);
        }

        arc e;
        for (auto vp{begin(env.brelaz.order) + env.l.clique_sz};
             vp != end(env.brelaz.order); ++vp) {
            auto v{*vp};
            if (env.brelaz.color[v] < env.l.clique_sz) {
                e = arc{v, env.brelaz.order[env.brelaz.color[v]]};
                break;
            }
        }
        return e;
    }

    void make_choice()
    {
        ++env.depth;
        env.G.trail.push_back(env.N);
        arc e{select()};
        env.G.contract(e[0], e[1]);
    }
};

template <class graph_struct> class coloring_algorithm
{

public:
    graph_struct& G;
    gc::statistics& stats;
    gc::options& options;

    int UB;
    int LB;
    int depth;
    int N;
    dsatur brelaz;

    heuristic<graph_struct> h;
    lower_bound<graph_struct> l;
    selector<graph_struct> s;

public:
    std::vector<int> coloring;

    coloring_algorithm(graph_struct& g, statistics& stats, gc::options& options)
        : G(g)
        , stats(stats)
        , options(options)
        , UB{static_cast<int>(g.size())}
        , LB{1 + (G.num_edges > 0)}
        , depth{0}
        , N{static_cast<int>(G.nodes.capacity())}
        , h(*this)
        , l(*this)
        , s(*this)
    {
        stats.notify_lb(LB);
        stats.notify_ub(UB);
        brelaz.use_recolor = false;
    }

    void update_lb(const int lb)
    {
        if (depth == 0 and lb > LB) {
            LB = lb;
            stats.notify_lb(lb);
            stats.custom_force_display(std::cout);
        }
    }

    void update_ub(const int ub)
    {
        if (ub < UB) {
            UB = ub;
            stats.notify_ub(ub);
            stats.custom_force_display(std::cout);
        }
    }

    // void contract(arc e) {
    // 	trail.push_back(limit);
    // 	G.contract(e[0], e[1]);
    // }

    void backtrack()
    {
        --depth;
        arc e{G.backtrack(N)};
        G.addition(e[0], e[1]);
        ++stats.total_conflicts;
    }

    void find_coloring()
    {
        int period = std::pow(10, 6 - options.verbosity);

        update_ub(h.degeneracy_coloring(coloring));

        update_ub(h.dsatur_coloring(UB, coloring));

        while (LB < UB) {

            if (options.verbosity > 1
                && stats.notify_iteration(depth) % period == 0)
                stats.custom_force_display(std::cout);

            auto lb{l.clique(begin(brelaz.order), end(brelaz.order))};
            update_lb(lb);

            update_ub(h.dsatur_coloring(UB+1, coloring));

            if (UB > lb) {

                s.make_choice();

            } else if (depth > 0) {

                backtrack();

            } else
                break;
        }
				
				stats.custom_force_display(std::cout);
    }
}; 


} // namespace gc

#endif // ALGORITHM
