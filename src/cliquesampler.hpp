#include <iostream>

#include "graph.hpp"

#ifndef __CLIQUESAMPLER_HPP
#define __CLIQUESAMPLER_HPP

namespace gc
{

template <class adjacency_struct> struct clique_sampler {

    adjacency_struct start_set;
    adjacency_struct cand_set;

    template <class viterator>
    int find_clique(basic_graph<adjacency_struct>& g, const int l,
        viterator first, viterator last, const int basewidth=64)
    {
        int lb_old{l}, lb{l};
				
				start_set.clear();
				
				for(auto vp{first}; vp!=last; ++vp) {
						start_set.push_back(*vp);
						std::cout << " " << (*vp);
				}
				std::cout << std::endl;
				
				return lb;
    }
};

} // namespace gc

#endif // __CLIQUESAMPLER_HPP
