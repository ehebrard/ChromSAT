#ifndef __CG_SPARSE_GRAPH_HH
#define __CG_SPARSE_GRAPH_HH

#include "bitset.hpp"
#include "intstack.hpp"

#include <algorithm>
#include <list>
#include <vector>


// #define _DEBUG_CLIQUE

namespace gc
{

using weight = int64_t;
using bitset = BitSet;

/// UNUSED 
using edge = std::pair<int,int>;

// NEW STUCTURE (use IntStack instead ?)
using vertices_vec = std::vector<int>;


/****************************************************************************
					CLASS FOR SPARSE GRAPH

using adjacency vector describing neighbors of each vertex
*****************************************************************************/

class sparse_graph
{
public:
	vertices_vec nodeset;
    IntStack nodes;
	std::vector<vertices_vec> adjacency;
	// we keep a copy of the original adjacency because we modify it
    // when we do merge/separate
	std::vector<vertices_vec> origadjacency;
	
	// Constructors
    sparse_graph() {}
	//~sparse_graph() {clear();}
	
	// With nv vertices, empty but reserved adjacency
    explicit sparse_graph(int nv)
    {
        nodes.reserve(nv);
        nodes.fill(); // Set nodes.size_ at nodes.list_.size()...
		adjacency.reserve(nv); // Vector of size nv of empty vector, no neighboors yet
		origadjacency.reserve(nv);
		vertices_vec temp;
		for (int i = 0; i < nv; ++i){
			adjacency.push_back(temp);
			origadjacency.push_back(temp);
		}
    }
	// What are those ?
    sparse_graph(sparse_graph&) = default;
    sparse_graph(sparse_graph&&) = default;
    sparse_graph& operator=(const sparse_graph&) = default;
    sparse_graph& operator=(sparse_graph&&) = default;
	
	int capacity() const { return adjacency.size(); } // adjacency.size() instead of matrix.size()... same thing ?

	// Adding each other in their neighbors
	void add_edge(const int u, const int v)
	{
		adjacency[u].push_back(v);
		adjacency[v].push_back(u);
		origadjacency[u].push_back(v);
		origadjacency[v].push_back(u);
	}

	// if u is not already in adjacency[v]	
	void safe_add_edge(const int u, const int v)
	{	 
		std::vector<int>::iterator it;
		it = find(adjacency[v].begin(), adjacency[v].end(), u);
  		if (it == adjacency[v].end()) add_edge(u,v); 
	}

	// TO STUDY
    void add_node(const int v)
    {
        nodes.add(v);
    }
	
    void remove_node(int v)
    {
        nodes.remove(v);
    }

	// Argument ?
	/*
    void add_clique(const bitset& C)
    {
        for (auto v : C) {
            add_node(v);
            matrix[v].union_with(C);
            matrix[v].remove(v);
        }
    }
	*/

	// Ok
    void clear()
    {
        for (auto v : nodes) {
            adjacency[v].clear();
        }
		adjacency.clear();
        nodes.clear();
    }

	void describe(std::ostream& os) const;

	std::ostream& display_adjacency(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const sparse_graph& x);


/****************************************************************/
//					SPARSE CLIQUE FINDER						//
/****************************************************************/

struct sparse_clique_finder {
	const sparse_graph& g;
    std::vector<vertices_vec> cliques;
    std::vector<int> clique_sz;
    std::vector<vertices_vec> candidates; // Neighbors of the clique members
    std::vector<int> last_clique;
    int num_cliques;

    sparse_clique_finder(const sparse_graph& g);
	
	// clear previously cached results
	void clear();
	// initialize a new clique
	void new_clique();
    // initialize a new color
    void new_color();
    // insert v into the clq^th clique. assumes it fits
    void insert(int v, int clq);
    // insert v into the col^th color. assumes it fits. Puts vertices
    // added from candidates[i] into diff
    // void insert_color(int v, int col, bitset& diff);
    // heuristically find a set of cliques and return the size of the
    // largest

    template <class ordering> int find_cliques(ordering o, const int limit=0xfffffff)
    {
        clear();
        if (o.size() == 0)
            return 0;
        for (auto u : o) {
            bool found{false};
            for (int i = 0; i != num_cliques; ++i)
                if (std::find(candidates[i].begin(), candidates[i].end(), u) != candidates[i].end()) {
                    found = true;
                    insert(u, i);
                }
            if (!found && num_cliques < limit) {
                new_clique();
                insert(u, num_cliques - 1);
            }
        }

        for (auto u : o) {
            for (int i = last_clique[u] + 1; i < num_cliques; ++i)
                if (std::find(candidates[i].begin(), candidates[i].end(), u) != candidates[i].end()) {
                    insert(u, i);
                }
        }

        return *std::max_element(
            begin(clique_sz), begin(clique_sz) + num_cliques);
    }
	
	void display();
};

/****************************************************************/
//					BRONKERBOSCH CLIQUE FINDER					//
/****************************************************************/

// Finding all maximal cliques (recursive)
struct BronKerbosch {
	const sparse_graph& g;

	// Collected results
	std::vector<vertices_vec> maximal_cliques;
	std::vector<int> clique_sz;
	int num_cliques;
	
	// Algo variables (how to backtrack ?)
	vertices_vec actual_clique;
	vertices_vec candidates;
	vertices_vec banned;

	// Vertices ordering
	// Degree
	vertices_vec by_degree;
	std::vector<int> degree;
	void order_by_degree();

	// Degeneracy
	

	// Constructor
	BronKerbosch(const sparse_graph& g);
	
	void clear();
	
	void add_max_clique(vertices_vec& clique);
	
	void find_cliques_withoutPivot(vertices_vec clique, vertices_vec candidates, vertices_vec banned); 
	void find_cliques_withPivot(vertices_vec clique, vertices_vec candidates, vertices_vec banned);	

	void display();

	void bronkerbosch_calls_display(vertices_vec clique, vertices_vec candidates, vertices_vec banned);

	// Intersection, union
	vertices_vec intersect(vertices_vec const & v1, vertices_vec const & v2);
	vertices_vec unite_vector_element(vertices_vec & v,const int i);
};

} // namespace gc

#endif
