#include "vertices_vec.hpp"

// TODO
/* Copying : clear vector before ? use  = , std::copy, or else ?
/* 
*/ 


namespace gc
{


/********************************************************/
//					VERTICES_VEC						//
/********************************************************/
vertices_vec::vertices_vec(){}
vertices_vec::vertices_vec(vertices_vec const& v):vertices(v.vertices){}

vertices_vec::~vertices_vec(){}

// Intersect with (Assuming both vectors are sorted)
void vertices_vec::intersect_with(vertices_vec const& v)
{	
	// Put intersection in the buffer
	buffer.clear();
	std::set_intersection(vertices.begin(), vertices.end(), v.vertices.begin(), v.vertices.end(), std::back_inserter(buffer));
	
	// Copy the buffer in vertices
	//vertices.clear();
	vertices = buffer;	
}

void vertices_vec::safe_intersect_with(vertices_vec & v)
{
	if (!v.is_sorted()) v.sort_by_number();
	if (!is_sorted()) sort_by_number();
	intersect_with(v);
}

// Union with another vertices_vec
void vertices_vec::union_with(vertices_vec const& v)
{
	// Put union in the buffer
	buffer.clear();
	std::set_union(vertices.begin(), vertices.end(), v.vertices.begin(), v.vertices.end(), std::back_inserter(buffer));
	
	// Copy the buffer in vertices
	//vertices.clear();
	vertices = buffer;	
}

void vertices_vec::safe_union_with(vertices_vec & v)
{
	if (!v.is_sorted()) v.sort_by_number();
	if (!is_sorted()) sort_by_number();
	union_with(v);
}

// Union with a single vertex
void vertices_vec::union_with(int i)
{	
	if (!is_sorted()) sort_by_number();
	vertices.insert(std::upper_bound(vertices.begin(), vertices.end(), i), i);
}


// Difference with
void vertices_vec::difference_with(vertices_vec const& v)
{
	// Put difference in the buffer
	buffer.clear();
	std::set_difference(vertices.begin(), vertices.end(), v.vertices.begin(), v.vertices.end(), 
                        std::inserter(buffer, buffer.begin()));
	// Copy the buffer in vertices
	//vertices.clear();
	vertices = buffer;
}

void vertices_vec::safe_difference_with(vertices_vec & v)
{
	if (!v.is_sorted()) v.sort_by_number();
	if (!is_sorted()) sort_by_number();
	difference_with(v);
}

/*
vertices_vec vertices_vec::union_with(vertices_vec const& v)
{
	vertices_vec union_vec(this);
	union_vec.union_with(v);
	return union_vec;
}

vertices_vec vertices_vec::safe_union_with(vertices_vec & v)
{
	vertices_vec union_vec(this);
	union_vec.safe_union_with(v);
	return union_vec;
}

vertices_vec vertices_vec::union_with(int i)
{
	vertices_vec union_vec(this);
	union_vec.union_with(i);
	return union_vec;
}
*/

void vertices_vec::display() const {
	std::cout << "{ ";
	for(std::vector<int>::const_iterator it = vertices.begin(); it!= vertices.end(); ++it){
	std::cout << *it << " ";
	}
	std::cout << "}"<< std::endl;
}

/********************************************************/
//					NEIGHBORS_VEC						//
/********************************************************/

int neighbors_vec::degree()
{
	return vertices.size();
}
} // namespace gc


