#ifndef __CG_CA_GRAPH_HH
#define __CG_CA_GRAPH_HH

#include "intstack.hpp"

#include <vector>

namespace gc
{


	template <class adjacency_struct> class ca_graph
	{
	public:
	    intstack nodes;
	    std::vector<std::vector<int>> neighbors;

			// struct for merge, merged nodes point to their parent, remaining nodes point to themselves
			std::vector<int> parent;
			
			// struct for merge, rank of the tree whose root is v
			// struct for merge, size of the tree whose root is v

	    ca_graph() {}
	    explicit ca_graph(int nv)
	        : neighbors(nv)
	    {
	        nodes.reserve(nv);
	        nodes.fill();
	    }
			
			// find does not do the "flatten" operation to make backtracking easier
			const find(const int v);
			
	};



} // namespace gc

#endif
