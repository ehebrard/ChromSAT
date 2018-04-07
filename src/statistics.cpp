#include "statistics.hpp"


namespace gc
{
		
void statistics::notify_bound_delta(const int b1, const int b2)
{
	total_bound_1 += b1;
	total_bound_2 += b2;
}

double statistics::get_bound_increase() const {
		return (double)total_bound_2/(double)total_bound_1;
}

void statistics::describe(std::ostream& os) const
{
    os << "GC statistics\n";
    os << "avg bound increase = " << get_bound_increase() << "\n";
    os << std::endl;
}

} // namespace gc
