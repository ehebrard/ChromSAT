#ifndef __CG_BGRAPH_HH
#define __CG_BGRAPH_HH

#include "bitset.hpp"
#include "intstack.hpp"

#include "Heap.h"

#include <vector>
#include <algorithm>


namespace gc
{

using bitset = BitSet;

class graph
{
public:
    bitset nodeset;
    IntStack nodes;
    std::vector<bitset> matrix;

    graph() {}
    explicit graph(int nv)
        : nodeset(0, nv - 1, bitset::empt)
        , matrix(nv)
    {
        nodes.reserve(nv);
        nodes.fill();
        for (auto& bs : matrix) {
            bs.initialise(0, nv, bitset::empt);
        }
    }
    graph(graph&) = default;
    graph(graph&&) = default;
    graph& operator=(const graph& g)
    {
        nodes.clear();
        nodeset.clear();
        for (auto v : g.nodes) {
            add_node(v);
            matrix[v].copy(g.matrix[v]);
        }
        return *this;
    }
    graph& operator=(graph&&) = default;

    int capacity() const { return matrix.size(); }

    void add_edge(const int u, const int v)
    {
        matrix[u].add(v);
        matrix[v].add(u);
    }

    void add_edges(const int u, const bitset& N)
    {
        matrix[u].union_with(N);
        for (auto v : N) {
            matrix[v].add(u);
        }
    }

    void add_node(const int v)
    {
        nodes.add(v);
        nodeset.add(v);
    }

    void add_clique(const bitset& C)
    {
        for (auto v : C) {
            add_node(v);
            matrix[v].union_with(C);
            matrix[v].remove(v);
        }
    }

    void remove_edge(int u, int v)
    {
        matrix[u].remove(v);
        matrix[v].remove(u);
    }

    void remove_node(int v)
    {
        nodes.remove(v);
        nodeset.fast_remove(v);
    }

    void clear()
    {
        for (auto v : nodes) {
            matrix[v].clear();
        }
        nodeset.clear();
        nodes.clear();
    }
};


struct clique_finder {
    const graph& g;
    std::vector<bitset> cliques;
    std::vector<int> clique_sz;
    std::vector<bitset> candidates;
    std::vector<int> last_clique;
    int num_cliques;

    clique_finder(const graph& g);

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
    void insert_color(int v, int col, bitset& diff);
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
                if (candidates[i].fast_contain(u)) {
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
                if (candidates[i].fast_contain(u)) {
                    insert(u, i);
                }
        }

        return *std::max_element(
            begin(clique_sz), begin(clique_sz) + num_cliques);
    }
};

std::vector<int> brelaz_color(const graph& g);


} //namespace

#endif
