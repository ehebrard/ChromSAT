#ifndef __GC_STATISTICS_HPP
#define __GC_STATISTICS_HPP

#include <iostream>
#include "prop.hpp"
// #include <minicsp/core/solver.hpp>


using namespace minicsp;


namespace gc
{
struct statistics {
	
		statistics(const int size) {
			
				changed = true;
				cons= NULL;
			
				// total_time = 0;
				total_conflicts = 0;
				best_lb = 0;
				best_ub = size;
				
				update_lb = true;
				update_ub = true;
							
			
				total_bound_1 = 0;
				total_bound_2 = 0;
				num_neighborhood_contractions = 0;
				num_vertex_removals = 0;
		}
	
    // outputs a nice description of all statistics
    void describe(std::ostream&);

		void binds(cons_base* c);
		void unbinds();

		cons_base *cons;
		
		bool changed;
		
		int best_lb; 
		void notify_lb(const int l);
		int best_ub;
		void notify_ub(const int u);
			

		bool update_lb;
		bool update_ub;

		// double total_time;
		uint64_t total_conflicts;


    // the actual statistics
		uint64_t num_neighborhood_contractions;
		
		uint64_t num_vertex_removals;
			
    uint64_t total_bound_1;
		uint64_t total_bound_2; 
		void notify_bound_delta(const int b1, const int b2);
		double get_bound_increase() const;
		
};

} // namespace gc

#endif
