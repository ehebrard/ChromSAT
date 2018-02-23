#ifndef __CG_GRAPH_HH
#define __CG_GRAPH_HH

#include "bitset.hpp"
#include "intstack.hpp"

#include <vector>

namespace gc
{

using weight = int64_t;
using bitset = BitSet;

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
    std::vector<std::vector<bitset>> diffs;
    std::vector<bitset> dirty;
    std::vector<int> removed;

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
    bitset util_set;

    //--------------------------------------------------
    // private, but out in the open

    // add an edge that will be backtracked later
    void add_dirty_edge(int u, int v);
public:
    graph() {}
    explicit graph(int nv)
    {
        matrix.resize(nv);
        origmatrix.resize(nv);
        nodeset.initialise(0, nv, bitset::full);
        // nodes.initialise(0, nv, nv, true);
				nodes.reserve(nv);
				nodes.fill();
        for (auto& bs : matrix) {
            bs.initialise(0, nv, bitset::empt);
        }
        for (auto& bs : origmatrix) {
            bs.initialise(0, nv, bitset::empt);
        }
        rep_of.resize(capacity());
        partition.resize(capacity());
        for (auto v : nodes) {
            rep_of[v] = v;
            partition[v].push_back(v);
        }
        util_set.initialise(0, nv, bitset::empt);
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
    int num_cliques;

    clique_finder(const graph& g);

    // clear previously cached results
    void clear();
    // initialize a new clique
    void new_clique();
    // insert v into the clq^th clique. assumes it fits
    void insert(int v, int clq);
    // heuristically find a set of cliques and return the size of the
    // largest
    int find_cliques();
};


struct neighbors_wrapper {
    const graph& g;
    std::vector< IntStack > by_degree;
		std::vector< IntStack > neighbors;
		std::vector< int > degree;
		
		int size;
		
    neighbors_wrapper(const graph& g);

    // synchronize the neighbors structure with the current graph
    void synchronize();
    // heuristically find a set of cliques and return the size of the
    // largest
    void get_degeneracy_order( std::vector< int >& order  );
};

} // namespace gc

#endif
