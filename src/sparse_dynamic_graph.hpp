

#include <assert.h>
#include <iomanip>
#include <vector>

#include "bitset.hpp"
#include "intstack.hpp"
#include "interval_list.hpp"


#ifndef __DYNGRAPH_HPP
#define __DYNGRAPH_HPP


namespace gc
{
	
	using bitset = BitSet;
	
	template <class T>
	struct Pair
	{

	public:
	    T data[2];

	    Pair(T a, T b)
	    {
	        data[0] = a;
	        data[1] = b;
	    }

	    T operator[](const int i) const { return data[i]; }
	    T& operator[](const int i) { return data[i]; }
			
	    /*!@name Printing*/
	    //@{
	    std::ostream& display(std::ostream& os) const
	    {
	        os << "<" << data[0] << " " << data[1] << ">";
	        return os;
	    }
	    //@}
	};
	
	template <class T>
	std::ostream& operator<<(std::ostream& os, const Pair<T>& x) {
		return x.display(os);
	}

	template <class T>
	std::ostream& operator<<(std::ostream& os, const Pair<T>* x) {
		return x->display(os);
	}
	
	using Edge = Pair<int>; 
	


template <typename T>
std::ostream& vecdisplay(const std::vector<T> vec, std::ostream& os)
{

    os << "[";
    for (typename std::vector<T>::const_iterator it = vec.begin();
         it != vec.end(); ++it) {
        os << " " << *it;
    }
    os << " ]";

    return os;
}



class dyngraph
{

public:
    /** nodes */
    // the maximum capacity of the graph (number of nodes)
    size_t capacity;
    // the list of nodes of the graph, in no particular order
    intstack node;
    // the nodes in a bitset, to help bitwise operations
    gc::bitset nodeset;

    /** edges */
    // the current number of edges (the list of edges is not updated)
    int num_edges;
    // the (original) edges of the graph
    std::vector<Edge> edges;
    // ranks stores the indices of the nodes in eachother's neighbor list so
    // that we can remove the edge efficiently
    std::vector<Edge> ranks;

    /** neighborhood */
    // the list of neighbors of each node
    std::vector<std::vector<int>> neighbor;
    // the indices of the incident edges
    std::vector<std::vector<int>> nb_index;


    dyngraph() {}
    dyngraph(const int n);
    dyngraph(const dyngraph& g);

    // helpers
    int size() const;
    bool null() const;
    bool empty() const;
    bool full() const;
    bool has_node(int x) const;
    int degree(const int x) const { return neighbor[x].size(); }
		double get_density() const;

    // remove every edge and node
    void clear();

    // rename the vertices so that they are sorted by non-increasing degree
    void sort(bool non_decreasing=true);

    // node addition/removal
    void declare_node(const int x);
    void add_node(const int x);

    void rem_node(const int x);

    // edge addition/removal
    int add_edge(const int x, const int y);
    void rem_edge(const int i);

    // printing
    std::ostream& display(std::ostream& os) const;
    void print_dimacs(std::ostream& os) const;
		
		// debug
		void verify(const char* msg);
		
		
		void maximal_matching(std::vector<int>& matching, int& nmatch, std::vector<int>& indexlist);
		// void brelaz_color(coloring& col);
		// void degeneracy();
		
		
		
	private:
		void rem_edge(const int x, const int y, const int i);

};

std::ostream& operator<<(std::ostream& os, const dyngraph& x);

std::ostream& operator<<(std::ostream& os, const dyngraph* x);


struct coloring {

		std::vector<int> color;
		std::vector<int> order;
		std::vector<int> rank;
		std::vector<int> first;
		std::vector<interval_list> satur;
		
		void brelaz_color(dyngraph& g);
		
		void remove(const int y, const int d);
		
		void clear() {
				first.clear();
				for(auto v : order) 
						satur[v].clear();
		}

};


} // namespace gc

#endif // __DYNGRAPH_HPP
