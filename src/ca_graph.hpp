#ifndef __CG_CA_GRAPH_HH
#define __CG_CA_GRAPH_HH

#include <limits>

#include "bi_graph.hpp"
#include "cliquer.hpp"
#include "dsatur.hpp"
#include "intstack.hpp"
#include "statistics.hpp"

#include <vector>

namespace gc
{

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
    std::vector<int> rank;

    // struct for merge, size of the tree whose root is v

    ca_graph() {}
    explicit ca_graph(int nv)
        : matrix(nv)
        , num_edges(0)
        , parent(nv)
        , rank(nv)
    {
        nodes.reserve(nv);
        nodes.fill();

        nodeset.initialise(0, nv - 1, gc::bitset::full);

        for (auto v : nodes) {
            parent[v] = v;
            rank[v] = 1;
        }
    }

    size_t size() const { return nodes.size(); }
    size_t capacity() const { return nodes.capacity(); }

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


    void get_subproblem(std::vector<int>& vertices, const int size_limit);

    // void split(const int u, const int v);

    // find does not do the "flatten" operation to make backtracking easier
    int find(const int v) const;

    std::ostream& describe(std::ostream& os, const int verbosity) const;

    // void print() const;

    void check_consistency(const char* msg);
};

std::ostream& operator<<(std::ostream& os, const ca_graph& g);
std::ostream& operator<<(std::ostream& os, const bi_graph& g);

} // namespace gc

#endif
