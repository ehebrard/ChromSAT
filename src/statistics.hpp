#ifndef __GC_STATISTICS_HPP
#define __GC_STATISTICS_HPP

#include <iostream>

namespace gc
{
struct statistics {
	
		statistics() {
				total_bound_1 = 0;
				total_bound_2 = 0;
				num_neighborhood_contractions = 0;
		}
	
    // outputs a nice description of all statistics
    void describe(std::ostream&) const;

    // the actual statistics
		long long int num_neighborhood_contractions;
		
		
    long long int total_bound_1;
		long long int total_bound_2; 
		void notify_bound_delta(const int b1, const int b2);
		double get_bound_increase() const;
		
};

} // namespace gc

#endif
