#include <iostream>

#include "dimacs.hpp"
#include "graph.hpp"
#include "mycielski.hpp"
#include "options.hpp"
#include "rewriter.hpp"
#include "statistics.hpp"
#include "utils.hpp"
#include "vertices_vec.hpp"


template<class adjacency_struct>
void histogram(gc::graph<adjacency_struct>& g)
{
    std::vector<int> degrees;
    for (auto v : g.nodes) {
        degrees.push_back(g.matrix[v].size());
    }
    std::sort(begin(degrees), end(degrees));
    int psize = degrees.size() / 10;
    int pstart{0};
    while (static_cast<size_t>(pstart) < degrees.size()) {
        int pend
            = std::min(static_cast<size_t>(pstart + psize), degrees.size() - 1);
        while (static_cast<size_t>(pend) < degrees.size() - 1
            && degrees[pend + 1] == degrees[pend])
            ++pend;
        std::cout << (pend - pstart + 1) << " vertices: degree "
                  << degrees[pstart] << "-" << degrees[pend] << "\n";
        pstart = pend + 1;
    }
}

// store the pruning done on the original graph in order to trace them back (e.g. for printing a coloring)
template< class adjacency_struct >
struct graph_reduction {
    const gc::graph<adjacency_struct>& residual;
    std::vector<int> removed_vertices;
    gc::bitset nodeset;
    gc::bitset util_set;

    explicit graph_reduction(const gc::graph<adjacency_struct>& g)
        : residual(g)
        , nodeset(0, g.capacity(), gc::bitset::empt)
        , util_set(0, g.capacity(), gc::bitset::empt)
    {
    }

    int extend_solution(std::vector<int>& col, const gc::graph<adjacency_struct>& g)
    {
        int maxc{0};
        nodeset.copy(residual.nodeset);
        for (auto v : residual.nodes)
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
        // if (options.preprocessing == gc::options::NO_PREPROCESSING)
        //     return gr;

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

						// std::cout << "iteration " << niteration << std::endl;
						
            removed = false;
						
						// compute an upper bound
            auto sol{gc::brelaz_color(g)};
            for (auto u : g.nodes)
                for (auto v : g.matrix[u])
                    assert(sol[u] != sol[v]);
            int hub{*max_element(begin(sol), end(sol)) + 1};
            if (ub < 0 || (hub < ub && hub >= lb)) {
                ub = hub;
                statistics.notify_ub(ub);
            }

						// compute a lower bound
            hlb = cf.find_cliques(g.nodes);
            if (myciel)
                hlb = mf.improve_cliques_larger_than(lb);
            if (hlb > lb) {
                lb = hlb;
                statistics.notify_lb(lb);
            }
            statistics.display(std::cout);

						// remove nodes whose neighborhood is not of size at least lb [they can take any of the lb-|N(v)| colors not taken by their neighbors]
            forbidden.clear();
            toremove.clear();
            for (auto u : g.nodes) {
                if (forbidden.fast_contain(u))
                    continue;
								
								// TODO: REPLACE BY A "NEIGHBORHOOD_SIZE" METHOD
                if (g.matrix[u].size() >= static_cast<size_t>(lb))
                    continue;
								
								// v 
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
		            for (auto v : g.matrix[u]) {
		                g.matrix[v].remove(u);
		                g.origmatrix[v].remove(u);
		            }
            }
        } while (removed);
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
	
	
	
	  auto sol{gc::brelaz_color(g)};
	  for (auto u : g.nodes)
	      for (auto v : g.matrix[u])
	          assert(sol[u] != sol[v]);
	  int ub{*max_element(begin(sol), end(sol)) + 1};
		stat.notify_ub(ub);
		stat.display(std::cout);
	

    gc::clique_finder<adjacency_struct> cf{g};
    gc::mycielskan_subgraph_finder<adjacency_struct> mf(g, cf, false);
    int lb{cf.find_cliques(g.nodes)};

    if (myciel)
        lb = mf.improve_cliques_larger_than(lb - 1);



    stat.notify_lb(lb);
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
		g.sort();

		gc::graph<gc::vertices_vec> h(g);
		

		gc::statistics statistics(g.capacity());
		
		
		std::cout << "\noriginal graph =\n";
		histogram(g);

		gc::degeneracy_finder<gc::vertices_vec> df(g);
		df.degeneracy_ordering();
		df.display_ordering();
		
   		std::pair<int, int> bounds{0, g.capacity()};
		// bounds = initial_bounds(g, statistics, options.boundalg != gc::options::CLIQUES);

		gc_preprocessor<gc::vertices_vec> p(g, options, statistics, bounds);
		

		if(g.size() == 0) {
			
			std::cout << "\nresidual graph is empty\n";
			
			std::vector<int> sol(g.capacity(), -1);
			auto ncolors = p.reduction.extend_solution(sol,h);
			std::cout << "coloring = "<< ncolors << std::endl;
			
		  for (auto u : h.nodes)
		      for (auto v : h.matrix[u])
		          assert(sol[u] != sol[v]);
			
		} else if(g.size() == h.size()) {
			
			std::cout << "\nno reduction\n";
			
		} else {
			
			std::cout << "\nreduced graph =\n";
			histogram(g);
			
		}


}
