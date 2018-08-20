#ifndef __GC_LEXBFS_HPP
#define __GC_LEXBFS_HPP

#include "graph.hpp"

/* Inspired from Habib & al. 
	"Lex-BFS and partition refinement ..." (2000) */

// TODO: lexBFS, chordality test, coloring
// No degree consideration here
// How are the iterators used ?

namespace gc
{

template <class graph_struct> struct lex_BFS_ordering {

    const graph_struct& orig;

    // result
    graph_struct g;
    std::vector<int> order;


    // buffers
    std::vector<std::list<int>::iterator> iterators;
    bitset ordered;
    std::list<std::list<int>> buckets;
    bitset util_set;
	std::vector<int> label;


    lex_BFS_ordering(const graph_struct& g)
		: orig(g)
		, g(g)
		, iterators(g.capacity())
		, label(g.capacity())
		, ordered(0, g.capacity() - 1, bitset::empt)
		, util_set(0, g.capacity() - 1, bitset::empt)
    {
    }
	
	// From Paul, Viennot 1997 "algos autour de lex-BFS"
	// Simplicial elimination order
	void simplicial_EO()
	{
		int n = g.size();
		for (auto v : g.nodes) {
			label[v] = -1;
		}
		
		for (int i = n; i > 0; --i) {
			// Choose v with label max
			auto it = max_element(std::begin(label), std::end(label));
			int v = int(it - label.begin());
			order.push_back(v);
			ordered.fast_add(v);
			label[v] = -2; // means v is ordered

			std::cout << "Labels in order : ";
			for (auto u : label){
				std::cout << u << " ";}
			std::cout << std::endl;	
			
			for (auto u : g.matrix[v]) {
				if (label[u] > -2) {
					label[u] += i;
				}
			}
		}
	}

//	void partition_lex_BFS()
//	{
//		int n = g.capacity();

//		for (auto v : g.nodes) {
//			buckets[0].push_front(v);
//			iterators[v] = buckets[0].begin();
//			label[v] = -1;
//		}

//		
//		while(true) {
//        	size_t i{0};
//            for (; i != buckets.size(); ++i) // go backward instead ?
//                if (buckets[i].size() == 1)
//                    break;
//            if (i == buckets.size()+1) 
//                break;
//			
//			// Pick vertex v from the last class made of non-visited vertices
//			auto v = buckets[i].back();
//			order.push_back(v);
//            buckets[i].pop_back();
//            ordered.fast_add(v);
//			--n;

//			for (int j = 0; j <= i; ++j) {
//				for (auto u : buckets[j]) {
//					util_set.copy(g.matrix[v]);
//					util_set.intersect_with(g.nodeset);
//				}
//				
//			}
//					
//		}
//	}




};


} // namespace gc

#endif
