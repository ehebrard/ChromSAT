#ifndef __CG_GRAPH_HH
#define __CG_GRAPH_HH

#include "bitset.hpp"
#include "intstack.hpp"
#include "minicsp/core/solver.hpp"

#include <algorithm>
#include <list>
#include <vector>


// #define _DEBUG_MYCIEL
// #define _DEBUG_CLIQUE

namespace gc
{

using weight = int64_t;
using bitset = BitSet;
using edge = std::pair<int,int>;



class basic_graph
{
public:
    bitset nodeset;
    IntStack nodes;
    std::vector<bitset> matrix;

    basic_graph() {}
    explicit basic_graph(int nv)
        : nodeset(0, nv - 1, bitset::empt)
        , matrix(nv)
    {
        nodes.reserve(nv);
        nodes.fill();
        for (auto& bs : matrix) {
            bs.initialise(0, nv, bitset::empt);
        }
    }
    basic_graph(basic_graph&) = default;
    basic_graph(basic_graph&&) = default;
    basic_graph& operator=(const basic_graph& g)
    {
        nodes.clear();
        nodeset.clear();
        for (auto v : g.nodes) {
            add_node(v);
            matrix[v].copy(g.matrix[v]);
        }
        return *this;
    }
    basic_graph& operator=(basic_graph&&) = default;

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


class graph
{
public:
    bitset nodeset;
    IntStack nodes;
    std::vector<bitset> matrix;
    // we keep a copy of the original matrix because we modify matrix
    // when we do merge/separate
    std::vector<bitset> origmatrix;

    // checkpointing
    int cur_ckpt{0};
    std::vector<int> removed;



    // trail of removed edges, including a lim which delimits
    // checkpoints
    std::vector<std::pair<int, int>> edge_trail;
    std::vector<std::size_t> edge_lim;

    // the partitions generated by merging vertices. initially each
    // partition is a singleton, so rep_of[i]=i, partition[i] =
    // {i}. When we merge u and v, we keep one vertex (say u) as
    // representative of the partition, set partition[u] to be the
    // union of partition[u] and partition[v] (which were previously
    // disjoint) and set rep_of[i]=u for all i in partition[v],
    // including v. At the end g has u in its vertex set, but none of
    // the other vertices of partition[u].
    std::vector<int> rep_of;
    std::vector<std::vector<int>> partition;
    // information to backtrack the above: previous rep_of and
    // previous partition sizes
    std::vector<std::vector<size_t>> partition_size_trail;
    std::vector<std::vector<int>> rep_of_trail;

    // buffers
    bitset util_set, diff2;
    bitset partu, partv;

    //--------------------------------------------------
    // private, but out in the open

    // add an edge that will be backtracked later
    void add_dirty_edge(int u, int v);

public:
    graph() {}
    explicit graph(int nv)
        : nodeset(0, nv - 1, bitset::full)
        , matrix(nv)
        , origmatrix(nv)
        , rep_of(nv)
        , partition(nv)
        , util_set(0, nv - 1, bitset::empt)
        , diff2(0, nv - 1, bitset::empt)
        , partu(0, nv - 1, bitset::empt)
        , partv(0, nv - 1, bitset::empt)
    {
        nodes.reserve(nv);
        nodes.fill();
        for (auto& bs : matrix) {
            bs.initialise(0, nv, bitset::empt);
        }
        for (auto& bs : origmatrix) {
            bs.initialise(0, nv, bitset::empt);
        }
        for (auto v : nodes) {
            rep_of[v] = v;
            partition[v].push_back(v);
        }
    }
    graph(graph&) = default;
    graph(graph&&) = default;
    graph& operator=(const graph&) = default;
    graph& operator=(graph&&) = default;

    int capacity() const { return matrix.size(); }

    void add_edge(int u, int v)
    {
        matrix[u].add(v);
        matrix[v].add(u);
        origmatrix[u].add(v);
        origmatrix[v].add(u);
    }

    // merge vertices and return the id of the new vertex (one of u,
    // v)
    int merge(int u, int v);
    // separate u and v. Just adds an edge, but it is reversible through
    // checkpointing
    void separate(int u, int v);

		int contractPreprocess();


    int checkpoint();
    void restore(int ckpt);
    int current_checkpoint() const { return cur_ckpt; }

    void describe(std::ostream& os) const;

    // debugging
    void check_consistency() const;
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



struct mycielskan_subgraph_finder {

public:
    bool prune;

    const graph& g;
    const clique_finder& cf;

    basic_graph explanation_subgraph;
    int explanation_clique;

    mycielskan_subgraph_finder(const graph& g, const clique_finder& cf, const bool prune);

    // extend the subgraph G into a mycielski of subsequent order if possible,
    // the additional vertices go into "subgraph"
    int extends(const bitset& G);

		int improve_cliques_larger_than(const int size);
    int full_myciel(const int lb, const int ub, minicsp::Solver& s,
        const std::vector<std::vector<minicsp::Var>>& vars);

    int improve_cliques_larger_than(const int size, const int lb, const int ub,
        minicsp::Solver& s, const std::vector<std::vector<minicsp::Var>>& vars);

    int improve_greedy(const int size, const int lb, const int ub,
        minicsp::Solver& s, const std::vector<std::vector<minicsp::Var>>& vars);

private:
    // [tmp in "extends] subgraph that we try to build
    basic_graph subgraph;

    // [tmp in "extends] neighborhood of u
    bitset neighbors_w;

    // [tmp in "extends] neighborhood of S_v
    bitset neighbors_Sv;

    // [tmp in "extends] non neighborhood of v
    bitset non_neighbors;

    // [tmp in "extends] store the potential extra nodes of the mycielski
    std::vector<int> extra;

    // [tmp in extends] store the index of extra where S_extra[i] ends
    std::vector<int> endS;

    // [tmp in "extends] the set of "candidates" (intersection of the neighbors
    // of neighbors_Sv)
    bitset candidates;

    // [tmp in extends] edges to add
    std::vector<edge> new_edges;
    std::vector<int> u_layer;

    int ith_node;

    bitset pruning;
    bitset real_pruning;
    std::vector<edge> new_pruning;

    vec<minicsp::Lit> reason;

    // tries to find the possible u's starting from the ith v and returns the
    // rank for which it fails (or subgraph.size() if it succeeds)
    int another_myciel_layer(const int ith);

    // select the u's given a w and
    void select_middle_layer(const int w, const int beg_node,
        const int end_node, std::vector<int>& U, std::vector<edge>& edges);

    minicsp::Clause* do_prune(
        minicsp::Solver& s, const std::vector<std::vector<minicsp::Var>>& vars);
};


struct neighbors_wrapper {
    const graph& g;
    std::vector<IntStack> by_degree;
    std::vector<IntStack> neighbors;
    std::vector<int> degree;

    bitset buffer;

    neighbors_wrapper(const graph& g);

    // synchronize the neighbors structure with the current graph
    void synchronize();
    // heuristically find a set of cliques and return the size of the
    // largest
    void get_degeneracy_order(std::vector<int>& order);

    // debugging
    void check_consistency() const;
};

std::vector<int> brelaz_color(const graph& g);

} // namespace gc

#endif
