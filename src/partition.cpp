
#include "partition.hpp"

    /*!@name Constructors*/
    //@{
    partition::partition(const size_t n = 0, const size_t m = 1) {
    	resize(m,n);
    }
		
		size_t partition::size() { return bag.size(); }

    void partition::resize(const size_t m, const size_t n) {
			bag.resize(m);
    	index_.resize(n, NOTIN);
    }

    /*!@name List Manipulation*/
    //@{
    void partition::move(const int elt, const int from, const int to) {
    	list_[from][index_[elt]] = list[from].back();
			_index[list[from].back()] = _index[elt];
			list[from].pop_back();
			add_elt(elt, to);
    }
		void partition::add_elt(const int elt, const int to) {
			index_[elt] = bag[to].size();
			bag[to].push_back(elt);
		}
		void partition::swap(const int a, const int b) {
			std::swap(bag[a].list_, bag[b].list_);
		}
		void partition::remove(const int a) {
			std::swap(bag[a].list_, bag.back().list_);
			bag.pop_back();
		}
    //@}

    /*!@name Miscellaneous*/
    //@{
    std::ostream& display(std::ostream& os) const {
    	for(auto it{begin(bag)}; it!=end(bag); ++it) {
    		os << (it - begin(bag)) << ":";
				for(auto jt{begin(*it)}; jt!=end(*it); ++jt) {
					os << " " << *jt;
				}
				os << std::endl;
    	}
			return os;
    }

		std::ostream& operator<<(std::ostream& os, const partition& x) {
			return x.display(os);
		}
