#include "statistics.hpp"


namespace gc
{
		
void statistics::notify_bound_delta(const int d)
{
	total_bound_delta += d;
	++num_bound_delta;
}

void statistics::describe(std::ostream& os) const
{
    os << "GC statistics\n";
    os << "avg bound delta = " << (double)total_bound_delta/(double)num_bound_delta << "\n";
    os << std::endl;
}

} // namespace gc
