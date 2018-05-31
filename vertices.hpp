#ifndef __CG_VERTICES_VEC_HH
#define __CG_VERTICES_VEC_HH

#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>

namespace gc
{

// eligible for adjacency_struct template
class vertices_vec
{
public:
	std::vector<int> vertices;
	std::vector<int> buffer;

	// Constructors
	vertices_vec();
	vertices_vec(vertices_vec const& v);
	~vertices_vec();

	void copy_from(vertices_vec const& v) {vertices = v.vertices;}

	// Std methods for vector
	bool empty() const {return vertices.empty();}
	size_t size() const {return vertices.size();}
	void push_back(int const& i) {vertices.push_back(i);}
	void reserve(size_t n) {vertices.reserve(n);}
	
	// Iterators
	std::vector<int>::iterator begin() {return vertices.begin();}
	std::vector<int>::iterator end() {return vertices.end();}
	std::vector<int>::const_iterator begin() const {return vertices.begin();}
	std::vector<int>::const_iterator end() const {return vertices.end();}
	

	void sort_by_number() {std::sort (vertices.begin(), vertices.end());}
	bool is_sorted() const {return std::is_sorted(vertices.begin(), vertices.end());}
	

	// Modifies this
	void intersect_with(vertices_vec const& v);
	void safe_intersect_with(vertices_vec & v); // Sorts v (!) 

	
	void union_with(vertices_vec const& v);
	void safe_union_with(vertices_vec & v);
	void union_with(int i); // Safed

	void difference_with(vertices_vec const& v);
	void safe_difference_with(vertices_vec & v); // Sorts v (!) 

	void display() const;
	
};

class neighbors_vec: public vertices_vec
{
	public:
		int degree();
};

} // namespace gc
#endif