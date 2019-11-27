#ifndef __CG_CA_GRAPH_HH
#define __CG_CA_GRAPH_HH

#include "dsatur.hpp"
#include "intstack.hpp"
#include "statistics.hpp"

#include <vector>

namespace gc
{

// template <class graph_struct> struct cliquer {
//     const graph_struct& g;
//     std::vector<std::vector<int>> cliques;
//     std::vector<gc::bitset> candidates;
//     std::vector<int> last_clique;
//     int num_cliques;
//     int limit;
//
//
//     cliquer(const graph_struct& ig, const int c = 0xfffffff)
//         : g(ig)
//         , num_cliques(1)
//         , limit(c)
//     {
//         last_clique.resize(g.capacity());
//     }
//
//
//     // clear previously cached results
//     void clear()
//     {
//         num_cliques = 0;
//     }
//     // initialize a new clique
//     void new_clique();
//
//     // insert v into the clq^th clique. assumes it fits
//     void insert(int v, int clq);
//
//     // heuristically find a set of cliques and return the size of the
//     // largest
//     template <class ForwardIterator>
//     int find_cliques(
//         ForwardIterator beg_ordering, ForwardIterator end_ordering)
//     {
//         clear();
// 				for (auto it{beg_ordering}; it != end_ordering;
// ++it)
// {
// 	          bool found{false};
// 	          for (int i = 0; i != num_cliques; ++i)
// 	              if (candidates[i].fast_contain(*it)) {
// 	                  found = true;
// 	                  insert(*it, i);
// 	              }
// 	          if (!found && num_cliques < limit) {
// 	              new_clique();
// 	              insert(*it, num_cliques - 1);
// 	          }
// 				}
//
// 				for (auto it{beg_ordering}; it != end_ordering;
// ++it)
// {
// 	          for (int i = last_clique[*it] + 1; i < num_cliques; ++i)
// 	              if (candidates[i].fast_contain(*it)) {
// 	                  insert(*it, i);
// 	              }
// 				}
//
//         auto maxclique{*std::max_element(
//             begin(clique_sz), begin(clique_sz) + num_cliques)};
//
//         return maxclique;
//     }
//
//
//
// 	};
//
//
// 	template <class graph_struct>
// 	void cliquer<graph_struct>::new_clique()
// 	{
// 		if(cliques.size() == num_cliques) {
// 			cliques.resize(num_cliques+1);
// 			candidates.resize(num_cliques+1);
// 			candidates.back().initialise(0, g.capacity(),
// bitset::empt)
// 		}
//
// 	    cliques[num_cliques].clear();
// 	    clique_sz[num_cliques] = 0;
// 	    candidates[num_cliques].copy(g.nodeset);
// 	    ++num_cliques;
// 	}
// 	// initialize a new color
// 	template <class graph_struct>
// 	void cliquer<graph_struct>::new_color()
// 	{
// 	    assert(num_cliques < g.capacity());
// 	    cliques[num_cliques].clear();
// 	    clique_sz[num_cliques] = 0;
// 	    candidates[num_cliques].clear();
// 	    ++num_cliques;
// 	}
// 	// insert v into the clq^th clique. assumes it fits
// 	template <class graph_struct>
// 	void cliquer<graph_struct>::insert(int v, int clq)
// 	{
// 	    cliques[clq].fast_add(v);
// 	    ++clique_sz[clq];
// 	    candidates[clq].intersect_with(g.matrix[v]);
// 	    last_clique[v] = clq;
//
// 	    // if (update_nn) {
// 	    //     for (auto u : g.matrix[v]) {
// 	    //         if (g.nodeset.fast_contain(u))
// 	    //             ++num_neighbors_of[clq][u];
// 	    //     }
// 	    // }
// 	}

template <class InputIterator1, class InputIterator2, class OutputIterator,
    class DeltaIterator, class InterIterator>
OutputIterator set_union_delta(InputIterator1 first1, InputIterator1 last1,
    InputIterator2 first2, InputIterator2 last2, OutputIterator result_union,
    DeltaIterator delta, // 2 \ 1
    InterIterator inter)
{
    while (first1 != last1 and first2 != last2) {
        if (*first1 < *first2) {
            *result_union = *first1;
            ++result_union;
            ++first1;
        } else {
            *result_union = *first2;
            ++result_union;
            if (*first2 < *first1) {
                *delta = *first2;
                ++first2;
                ++delta;
            } else {
                *inter = *first1;
                ++inter;
                ++first1;
                ++first2;
            }
        }
    }

    std::copy(first2, last2, delta);
    std::copy(first2, last2, result_union);

    return std::copy(first1, last1, result_union);
}

template <class ForwardIt, class T>
ForwardIt find_in(ForwardIt first, ForwardIt last, const T& value)
{
    ForwardIt it;
    typename std::iterator_traits<ForwardIt>::difference_type count, step;
    count = std::distance(first, last);

    while (count > 0) {
        it = first;
        step = count / 2;
        std::advance(it, step);
        if (*it < value) {
            first = ++it;
            count -= step + 1;
        } else if (*it == value)
            return it;
        else
            count = step;
    }
    return first;
}

class ca_graph
{
private:
    std::vector<int> union_buffer;
    std::vector<int> inter_buffer;

public:
    bitset nodeset;
    intstack nodes;
    std::vector<std::vector<int>> matrix;

    unsigned int num_edges;

    /*
    - When backtracking on an edge, simply remove u (v) for v's (u's)
    neighborhood

    - When backtracking on a contraction, get the list of neighbors added to u
    and the removed vertex v, then:
            * remove those vertices from u's neighborhood
            * remove u from each of thier respective neighborhood
            * add u to all of its neighbors
    */

    // if >=0, list of neighbors preceded by the remove node, otherwise -(u+1),
    // -(v+1) where u,v was the added edge
    std::vector<int> trail;

    // struct for merge, merged nodes point to their parent, remaining nodes
    // point to themselves
    std::vector<int> parent;

    // struct for merge, rank of the tree whose root is v
    // struct for merge, size of the tree whose root is v

    ca_graph() {}
    explicit ca_graph(int nv)
        : matrix(nv)
        , num_edges(0)
        , parent(nv)
    {
        nodes.reserve(nv);
        nodes.fill();
        for (auto v : nodes)
            parent[v] = v;
    }

    size_t size() const { return nodes.size(); }
    size_t capacity() const { return nodes.capacity(); }
		
		bool contain(const int v) const {return nodeset[v];}

    void add_neighbor(const int u, const int v);

    void addition(const int u, const int v);

    void add_edge(const int u, const int v);

    void contract(const int u, const int v);

    void remove_neighbor(const int u, const int v);

    void swap_neighbor(const int u, const int v, const int w);

    void canonize();

    arc undo();

    arc backtrack(const int plevel);

    arc undo_addition();

    arc undo_contraction();

    arc any_non_edge();

    void search(statistics& stat, options& options);

    // void split(const int u, const int v);

    // find does not do the "flatten" operation to make backtracking easier
    int find(const int v) const;

    std::ostream& describe(std::ostream& os, const int verbosity) const;

    // void print() const;

    void check_consistency();
};

std::ostream& operator<<(std::ostream& os, const ca_graph& g);

} // namespace gc

#endif
