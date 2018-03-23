#ifndef __GC_STATISTICS_HPP
#define __GC_STATISTICS_HPP

#include <iostream>

namespace gc
{
struct statistics {
	
		statistics() {
				total_bound_delta = 0;
				num_bound_delta = 0;
		}
	
    // outputs a nice description of all statistics
    void describe(std::ostream&) const;

    // the actual statistics
    long long int total_bound_delta; 
		long long int num_bound_delta;		
		void notify_bound_delta(const int d);
		
};

} // namespace gc

#endif
