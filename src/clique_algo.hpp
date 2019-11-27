
#include "coloring_heuristic.hpp"
#include "sorted_graph.hpp"
#include "sparse_set.hpp"
#include "options.hpp"
#include "graph.hpp"
#include "statistics.hpp"


#ifndef __CLIQUE_ALGO_HPP
#define __CLIQUE_ALGO_HPP

namespace gc
{

	class clique_algo {
		
	public:
		
		int ub;
		int lb;
		
		clique_algo(sorted_graph& g) : G(g), df(G) {
			
	    df.degeneracy_ordering();

	    ub = df.degeneracy + 1;
			
			lb = 0;
			
			current.reserve(G.capacity());
			current.fill();
			commit();
			
		}
		
		void initialise_lower_bound();
		
		int clique_size() const;
		
		void get_branches();
	
		bool dead_end();
		
		bool root_node();
		
		int bound();
	
		void branch();
		
		void search();
		
		void backtrack();
		
		void commit();
		
		void undo();
		
		void add_to_clique(const int v);
		
		
	private:
		
		std::vector<int> decision_stack;
		
		sparse_set current;
		
		std::vector<size_t> dtrail;
		
		sorted_graph& G;
		
		coloring_heuristic brelaz;
		
		degeneracy_finder<sorted_graph> df;
		
		std::vector<int> coloring;
		
		long long int num_fails{0};
	
	};

} // namespace gc

#endif // __CLIQUE_ALGO_HPP


