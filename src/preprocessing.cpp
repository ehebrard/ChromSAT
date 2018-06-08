#include <iostream>

#include "dimacs.hpp"
#include "graph.hpp"
#include "mycielski.hpp"
#include "options.hpp"
#include "rewriter.hpp"
#include "statistics.hpp"
#include "utils.hpp"
#include "vertices_vec.hpp"


template< class adjacency_struct >
struct graph_reduction {
    const gc::graph<adjacency_struct>& g;
    std::vector<int> removed_vertices;
    gc::bitset nodeset;
    gc::bitset util_set;

    explicit graph_reduction(const gc::graph<adjacency_struct>& g)
        : g(g)
        , nodeset(0, g.capacity(), gc::bitset::empt)
        , util_set(0, g.capacity(), gc::bitset::empt)
    {
    }

    int extend_solution(std::vector<int>& col)
    {
        int maxc{0};
        nodeset.copy(g.nodeset);
        for (auto v : g.nodes)
            maxc = std::max(maxc, col[v]);
        for (auto i = removed_vertices.rbegin(), iend = removed_vertices.rend();
             i != iend; ++i) {
            auto v = *i;
            util_set.clear();
            for (auto u : g.matrix[v]) {
                if (!nodeset.fast_contain(u))
                    continue;
                util_set.fast_add(col[u]);
            }
            for (int q = 0; q != g.capacity(); ++q) {
                if (util_set.fast_contain(q))
                    continue;
                maxc = std::max(maxc, q);
                col[v] = q;
                break;
            }
            nodeset.fast_add(v);
        }
        return maxc + 1;
    }
};


template< class adjacency_struct >
struct gc_preprocessor {
    int lb{0}, ub{-1};
    const gc::options& options;
    gc::statistics& statistics;

    graph_reduction<adjacency_struct> reduction;
    gc::graph<adjacency_struct>& g;

    graph_reduction<adjacency_struct> preprocess(
        gc::graph<adjacency_struct>& g, std::pair<int, int> bounds, bool myciel = false)
    {
        graph_reduction gr(g);
        if (options.preprocessing == gc::options::NO_PREPROCESSING)
            return gr;

        lb = bounds.first;
        ub = bounds.second;
        int hlb{0};
        gc::clique_finder<adjacency_struct> cf{g};
        gc::mycielskan_subgraph_finder<adjacency_struct> mf(g, cf, false);

        adjacency_struct forbidden(0, g.capacity(), gc::bitset::empt);
        adjacency_struct util_set(0, g.capacity(), gc::bitset::empt);
        adjacency_struct removedv(0, g.capacity(), gc::bitset::empt);
        std::vector<int> toremove;
        bool removed{false};
        int niteration{0};
        do {
            ++niteration;
            removed = false;
            auto sol{gc::brelaz_color(g)};
            for (auto u : g.nodes)
                for (auto v : g.matrix[u])
                    assert(sol[u] != sol[v]);
            int hub{*max_element(begin(sol), end(sol)) + 1};
            if (ub < 0 || (hub < ub && hub >= lb)) {
                ub = hub;
                statistics.notify_ub(ub);
            }

            hlb = cf.find_cliques(g.nodes);
            if (myciel)
                hlb = mf.improve_cliques_larger_than(lb);

            if (hlb > lb) {
                lb = hlb;
                statistics.notify_lb(lb);
            }
            statistics.display(std::cout);

            forbidden.clear();
            toremove.clear();
            for (auto u : g.nodes) {
                if (forbidden.fast_contain(u))
                    continue;
                util_set.copy(g.matrix[u]);
                util_set.intersect_with(g.nodeset);
                if (util_set.size() >= static_cast<size_t>(lb))
                    continue;
                removed = true;
                removedv.fast_add(u);
                ++statistics.num_vertex_removals;
                toremove.push_back(u);
                gr.removed_vertices.push_back(u);
                forbidden.union_with(g.matrix[u]);
            }
            for (auto u : toremove) {
                g.nodes.remove(u);
                g.nodeset.remove(u);
            }
        } while (removed);
        if (removedv.size() > 0) {
            for (auto v : g.nodes) {
                g.matrix[v].setminus_with(removedv);
                g.origmatrix[v].setminus_with(removedv);
            }
        }
        return gr;
    }

    gc_preprocessor(gc::graph<adjacency_struct>& g, const gc::options& options,
        gc::statistics& statistics, std::pair<int, int> bounds)
        : options(options)
        , statistics(statistics)
        , reduction(preprocess(g, bounds))
        , g(g)
    {
        auto plb = bounds.first;
        auto pub = bounds.second;

        lb = std::max(lb, plb);
        if (ub < 0 || pub < ub)
            ub = pub;
    }

    void print_stats()
    {
        statistics.display(std::cout);
        std::cout << std::endl;
    }
};

template< class adjacency_struct >
std::pair<int, int> initial_bounds(
    const gc::graph<adjacency_struct>& g, gc::statistics& stat, bool myciel = false)
{
    // gc::degeneracy_finder df{g};
    // df.degeneracy_ordering();
    // for( auto u : df.order ) {
    // 	std::cout << u << "(" << g.matrix[u].size() << ") ";
    // }
    // std::cout << std::endl;

    gc::clique_finder<adjacency_struct> cf{g};
    gc::mycielskan_subgraph_finder<adjacency_struct> mf(g, cf, false);
    int lb{cf.find_cliques(g.nodes)};

    if (myciel)
        lb = mf.improve_cliques_larger_than(lb - 1);

    auto sol{gc::brelaz_color(g)};
    for (auto u : g.nodes)
        for (auto v : g.matrix[u])
            assert(sol[u] != sol[v]);

    int ub{*max_element(begin(sol), end(sol)) + 1};
    stat.notify_lb(lb);
    stat.notify_ub(ub);
    stat.display(std::cout);
    return std::make_pair(lb, ub);
}


int main(int argc, char* argv[])
{
    auto options = gc::parse(argc, argv);
    options.describe(std::cout);
		
    
		gc::graph<gc::vertices_vec> g;
    dimacs::read_graph(options.instance_file.c_str(),
        [&](int nv, int) { g = gc::graph<gc::vertices_vec>{nv}; },
        [&](int u, int v) {
            if (u != v)
                g.add_edge(u - 1, v - 1);
        },
        [&](int, gc::weight) {});
    g.describe(std::cout);


		gc::statistics statistics(g.capacity());
		
    std::pair<int, int> bounds{0, g.capacity()};
		bounds = initial_bounds(g, statistics, options.boundalg != gc::options::CLIQUES);

		gc_preprocessor<gc::vertices_vec> p(g, options, statistics, bounds);

		// graph_reduction<gc::vertices_vec> gr =


}
