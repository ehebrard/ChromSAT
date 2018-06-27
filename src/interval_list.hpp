

#include <assert.h>
#include <iomanip>
#include <vector>




#ifndef __INTERVAL_LIST_HPP
#define __INTERVAL_LIST_HPP





namespace gc
{

class interval_list
{

public:
	
	static const int BIG = 0x1000000;
	
	std::vector<int> inf;
	std::vector<int> sup;
	std::vector<int> next;
	std::vector<int> freed;
	
	int size;
	
	
	interval_list();
	
	bool add(const int x);
	int get() const;
	void clear();
	
	std::ostream& display(std::ostream& os) const;

};

std::ostream& operator<<(std::ostream& os, const interval_list& x);

std::ostream& operator<<(std::ostream& os, const interval_list* x);

} // namespace gc

#endif // __INTERVAL_LIST_HPP
