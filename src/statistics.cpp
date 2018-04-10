
#include <iomanip>

#include "statistics.hpp"
#include <minicsp/core/utils.hpp>


namespace gc
{
		
void statistics::notify_bound_delta(const int b1, const int b2)
{
		total_bound_1 += b1;
		total_bound_2 += b2;
}

void statistics::notify_lb(const int l) 
{
		if(best_lb < l) {
				best_lb = l;
				changed = true;
		}
}

void statistics::notify_ub(const int u) 
{
		if(best_ub > u) {
				best_ub = u;
				changed = true;
		}
}

double statistics::get_bound_increase() const {
		if(total_bound_1)
				return (double)total_bound_2/(double)total_bound_1;
		return 0;
}

void statistics::describe(std::ostream& os)
{	
		if(update_lb && cons && best_lb < cons->bestlb) {
				changed = true;
				best_lb = cons->bestlb;
		}
		if(update_ub && cons && best_ub > cons->ub) {
				changed = true;
				best_ub = cons->ub;
		}
		
		
		
		if(changed) {
				os.setf(std::ios_base::fixed, std::ios_base::floatfield);
		    os << "d lb = " << std::setw(4) << std::left << best_lb
					 << "| ub = " << std::setw(4) << std::left << best_ub
					 << "| time = " << std::setw(10) << std::left << std::setprecision(4) << minicsp::cpuTime()
					 << "| conflicts = " << std::setw(10) << std::left << (cons ? total_conflicts + cons->s.conflicts : total_conflicts)
					 << "| delta = " << std::setw(8) << std::left << std::setprecision(4) << get_bound_increase() 
					 << "| #dom = " << std::setw(10) << std::left << num_neighborhood_contractions
					 << "| #rem = " << std::setw(10) << std::left << num_vertex_removals
		    	 << std::endl;
		}
		
		changed = false;
}

void statistics::binds( gc::cons_base* c ) {
		cons = c;
}

void statistics::unbinds() {
		total_conflicts += cons->s.conflicts;
		cons = NULL;
}

} // namespace gc
