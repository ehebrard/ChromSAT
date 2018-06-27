
#include <iostream>
#include "interval_list.hpp"



gc::interval_list::interval_list() { 
		inf.push_back(-0x1000000);
		inf.push_back(0x1000000);
		sup.push_back(-0x1000000);
		sup.push_back(0x1000000);
		next.push_back(1);
		next.push_back(-1);
	
		size = 0; 
}


bool gc::interval_list::add(const int x) {
	int prev=0, pos=next[0], idx;
	while(pos > 1 && sup[pos] < x) {
			prev = pos;
			pos = next[pos];
	}
	if(inf[pos] > x) {	
		  if(sup[prev] == inf[pos] - 2) { // merge
				sup[prev] = sup[pos];
				next[prev] = next[pos];
				freed.push_back(pos);
			} else if(sup[prev] == x - 1) {
				sup[prev] = x;
			} else if(inf[pos] == x + 1) {
				inf[pos] = x;
			} else { // new singleton interval
				if(freed.size() == 0) {
						inf.push_back(x);
						sup.push_back(x);
						idx = next.size();
						next.push_back(pos);
				} else {
						idx = freed.back();
						freed.pop_back();
						inf[idx] = x;
						sup[idx] = x; 
						next[idx] = pos;
						next[prev] = idx;
				}
				next[prev] = idx;
			}
		
			++size;
			return true;
	}
	
	return false;
}


int gc::interval_list::get() const {
		int c = sup[next[0]];
		return (c == 0x1000000 ? 0 : c + 1); 
}

void gc::interval_list::clear() {
	inf.resize(2);
	sup.resize(2);
	next.resize(2);
	next[0] = 1;
	size = 0;
}

/*!@name Printing*/
//@{
std::ostream& gc::interval_list::display(std::ostream& os) const
{
	
		os << "{";
		int prev=0, pos=next[0];
		while(pos > 1) {
				if(inf[pos] == sup[pos])
						os << " " << inf[pos];
				else
					os << " [" << inf[pos] << ".." << sup[pos] << "]";
				prev = pos;
				pos = next[pos];
		}
	
    os << " }";
		
		// os << std::endl;
		// for( auto s : inf ) {
		// 	os << std::setw(2) << s << " ";
		// }
		// os << std::endl;
		// for( auto s : sup ) {
		// 	os << std::setw(2) << s << " ";
		// }
		// os << std::endl;
		// for( auto s : next ) {
		// 	os << std::setw(2) << s << " ";
		// }
		// os << std::endl;
		
		
    return os;
}
//@}


std::ostream& gc::operator<<(std::ostream& os, const gc::interval_list& x)
{
    return x.display(os);
}

std::ostream& gc::operator<<(std::ostream& os, const gc::interval_list* x)
{
    return (x ? x->display(os) : os);
}
