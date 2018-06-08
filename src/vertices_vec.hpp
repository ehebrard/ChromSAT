#ifndef __CG_VERTICES_VEC_HH
#define __CG_VERTICES_VEC_HH

#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>



// class Bitset<unsigned int, float>;

class bitset;

namespace gc
{

//class BitSet;

// eligible for adjacency_struct template
class vertices_vec
{
public:
	std::vector<int> vertices;
	std::vector<int> buffer;

	// Constructors
	vertices_vec();
	vertices_vec(vertices_vec const& v);
	vertices_vec(const int lb, const int ub, const int p);
	~vertices_vec();

	void copy(vertices_vec const& v) {vertices = v.vertices;}
	void copy(bitset const& elts); // {vertices.clear(); for(auto v : elts) vertices.push_back(v);}

	// Std methods for vector
	bool empty() const {return vertices.empty();}
	size_t size() const {return vertices.size();}
	void clear() {vertices.clear();}
	void push_back(int i) {vertices.push_back(i);}
	void reserve(size_t n) {vertices.reserve(n);}
	
	// Iterators
	std::vector<int>::iterator begin() {return vertices.begin();}
	std::vector<int>::iterator end() {return vertices.end();}
	std::vector<int>::reverse_iterator rbegin() {return vertices.rbegin();}
	std::vector<int>::reverse_iterator rend() {return vertices.rend();}
	
	std::vector<int>::const_iterator begin() const {return vertices.begin();}
	std::vector<int>::const_iterator end() const {return vertices.end();}
	std::vector<int>::const_reverse_iterator crbegin() const {return vertices.crbegin();}
	std::vector<int>::const_reverse_iterator crend() const {return vertices.crend();}
	
  void initialise(const int lb, const int ub, const int p);

	inline void add(const int v) {vertices.push_back(v);}
	inline void fast_add(const int v) {vertices.push_back(v);}
  bool fast_contain(const int elt) const;

	void sort_by_number() {std::sort (vertices.begin(), vertices.end());}
	bool is_sorted() const {return std::is_sorted(vertices.begin(), vertices.end());}
	

	// Modifies this
	void intersect_with(vertices_vec const& v);
	void safe_intersect_with(vertices_vec & v); // Sorts v and vertices (!) 

	
	void union_with(vertices_vec const& v);
	void safe_union_with(vertices_vec & v);
	void union_with(int i); // Sorts vertices (!)

	void difference_with(vertices_vec const& v);
	void safe_difference_with(vertices_vec & v); // Sorts v and vertices (!)
	void setminus_with(vertices_vec & v); // Sorts v and vertices (!) 
	
	void add_interval(const int l, const int u);
	
	void remove(const int e);
	
	bool includes(vertices_vec const& v) const;
	bool intersect(vertices_vec const& v) const;
	
	
	int min() const;
	int max() const;

	// Return a third vertices_vec
/*
	vertices_vec intersect_with(vertices_vec const& v);
	
	vertices_vec union_with(vertices_vec const& v) const;
	vertices_vec safe_union_with(vertices_vec & v) const;
	vertices_vec union_with(int i) const;
*/	
	void display() const;
	
};

class neighbors_vec: public vertices_vec
{
	public:
		int degree();
};

} // namespace gc
#endif
