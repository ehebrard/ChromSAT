#include <iostream>
#include <iomanip>

#include "dimacs.hpp"
#include "graph.hpp"
#include "mycielski.hpp"
#include "options.hpp"
#include "rewriter.hpp"
#include "statistics.hpp"
#include "utils.hpp"
#include "vertices_vec.hpp"
#include "interval_list.hpp"
#include "sparse_dynamic_graph.hpp"

#include <minicsp/core/utils.hpp>


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


template< class adjacency_struct >
void preprocess(gc::options& options, gc::graph<adjacency_struct>& g) {
	
  dimacs::read_graph(options.instance_file.c_str(),
      [&](int nv, int) { g = gc::graph<adjacency_struct>{nv}; },
      [&](int u, int v) {
          if (u != v)
              g.add_edge(u - 1, v - 1);
      },
      [&](int, gc::weight) {});
  g.describe(std::cout);
	g.canonize();
	
	gc::statistics statistics(g.capacity());
	
	
	gc::degeneracy_finder<adjacency_struct> df(g);	
	df.degeneracy_ordering();
	
	for(auto v : df.order) {
		std::cout << std::setw(4) << v << " " << std::setw(3) << df.degrees[v] << " " << std::setw(3) << g.matrix[v].size() << std::endl;
	}
	
	gc::graph<adjacency_struct> h(g);
	
  std::pair<int, int> bounds{0, g.capacity()};
	
	
	// bounds = initial_bounds(g, statistics, options.boundalg != gc::options::CLIQUES);
	// gc::clique_finder<adjacency_struct> cf{g};
	// // bitset restriction(0, h.capacity(), bitset::empt);
	// for(auto vp = rbegin(df.order); vp != rend(df.order); ++vp) {
	// 	for(auto up = vp+1; up != rend(df.order); ++up) {
	// 		if(!h.matrix[*vp].fast_contain(*up)) {
	// 			g.nodeset.copy(h.matrix[*vp]);
	// 			g.nodeset.union_with(h.matrix[*up]);
	//
	//
	// 			if(g.nodeset.size() >= bounds.first) {
	// 	    			g.nodes.clear();
	// 				for(auto x : g.nodeset) {
	// 					g.nodes.add(x);
	// 				}
	//
	//
	// 				// std::cout << g.nodeset << std::endl;
	// 				auto maxclique = cf.find_cliques(g.nodes);
	// 				if(maxclique >= bounds.first) {
	// 					std::cout << *vp << " " << *up << std::endl;
	// 				} else {
	// 					std::cout << ".";
	// 					std::cout.flush();
	// 				}
	// 				cf.clear();
	// 			}
	// 		}
	// 	}
	// }
	
	
	
	std::cout << "\noriginal graph (" << g.size() << ") =\n";
	histogram(g);
	
  // std::pair<int, int> bounds{0, g.capacity()};
	// bounds = initial_bounds(g, statistics, options.boundalg != gc::options::CLIQUES);

	gc_preprocessor<adjacency_struct> p(g, options, statistics, bounds);
	

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
		
		std::cout << "\nreduced graph (" << g.size() << ") =\n";
		histogram(g);

		gc::degeneracy_finder<adjacency_struct> df(g);
		df.degeneracy_ordering();
		df.display_ordering();

	}

}



template< class adjacency_struct >
void upper_bound(gc::options& options, gc::graph<adjacency_struct>& g) {
	
	
	double tnow, t = minicsp::cpuTime();
	std::cout << "read file into sorted vectors...";
	std::cout.flush();
	
	int num_edges = 0;
  dimacs::read_graph(options.instance_file.c_str(),
      [&](int nv, int) { g = gc::graph<adjacency_struct>{nv}; },
      [&](int u, int v) {
          if (u != v) {
              g.add_edge(u - 1, v - 1);
							++num_edges;
					}
      },
      [&](int, gc::weight) {});
  // g.describe(std::cout);
	g.canonize();
	
	// gc::statistics statistics(g.capacity());
			
	tnow = minicsp::cpuTime();
	std::cout << (tnow - t) << std::endl << g.size() << " nodes, " << num_edges << " edges\n";
	std::cout << "compute brelaz coloring...";
	std::cout.flush();
	t = tnow;
	
	
  auto sol{gc::brelaz_color(g)};
  // for (auto u : g.nodes)
  //     for (auto v : g.matrix[u])
  //         assert(sol[u] != sol[v]);
  

	tnow = minicsp::cpuTime();
	std::cout << (tnow - t) << std::endl;
	t = tnow;
	
	std::cout << "num_colors = " << *max_element(begin(sol), end(sol)) + 1 << std::endl;
	


}


void upper_bound(gc::options& options, gc::dyngraph& g) 
{
		
  // dimacs::read_graph(options.instance_file.c_str(),
  //     [&](int nv, int) { g = gc::dyngraph{nv}; },
  //     [&](int u, int v) {
  //         if (u != v) {
  //             g.add_edge(u - 1, v - 1);
  // 					}
  //     },
  //     [&](int, gc::weight) {});
	
	double tnow, t = minicsp::cpuTime();
	std::cout << "read file into sorted vectors...";
	std::cout.flush();
	
	
	gc::graph<gc::vertices_vec> sg;
  dimacs::read_graph(options.instance_file.c_str(),
      [&](int nv, int) { sg = gc::graph<gc::vertices_vec>{nv}; },
      [&](int u, int v) {
          if (u != v) 
              sg.add_edge(u - 1, v - 1);
      },
      [&](int, gc::weight) {});
  // sg.describe(std::cout);
	sg.canonize();
	
	tnow = minicsp::cpuTime();
	std::cout << (tnow - t) << std::endl;
	std::cout << "copy into sparse graph...";
	std::cout.flush();
	t = tnow;

	g = gc::dyngraph{sg.capacity()};
	for(auto v : sg.nodes) 
		for(auto u : sg.matrix[v]) 
			if(u > v) 
				g.add_edge(v,u);
			
		
	tnow = minicsp::cpuTime();
	std::cout << (tnow - t) << std::endl << g.size() << " nodes, " << g.num_edges << " edges\n";
	std::cout << "compute brelaz coloring...";
	std::cout.flush();
	t = tnow;
	
	
	gc::coloring col;
	col.brelaz_color(g);
	
	tnow = minicsp::cpuTime();
	std::cout << (tnow - t) << std::endl;
	t = tnow;
	
	std::cout << "num_colors = " << *std::max_element(begin(col.color), end(col.color)) + 1 << std::endl;

}



int main(int argc, char* argv[])
{
    auto options = gc::parse(argc, argv);
    options.describe(std::cout);
		
		if(options.preprocessing == gc::options::SPARSE) {
				gc::graph<gc::vertices_vec> g;
				preprocess(options, g);
		} else if(options.preprocessing == gc::options::LOW_DEGREE) {
				gc::graph<gc::bitset> g;	
				preprocess(options, g);
		} else {
			std::cout << "\nsparse graph\n";
				gc::dyngraph g;
				upper_bound(options, g);
				
			std::cout << "\nbitset graph\n";
				gc::graph<gc::bitset> h;	
				upper_bound(options, h);
		}
		
}
