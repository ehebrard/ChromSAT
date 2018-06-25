#ifndef __CG_GRAPH_HH
#define __CG_GRAPH_HH

#include "bitset.hpp"
#include "intstack.hpp"

#include <minicsp/mtl/Heap.h>

#include <algorithm>
#include <list>
#include <vector>


// #define _DEBUG_CLIQUE
// #define SUBCLASS

namespace gc
{

using weight = int64_t;
using bitset = BitSet;
using edge = std::pair<int,int>;




template< class adjacency_struct >
class basic_graph
{
public:
    bitset nodeset;
    IntStack nodes;
    std::vector<adjacency_struct> matrix;

    basic_graph() {}
    explicit basic_graph(int nv)
        : nodeset(0, nv - 1, bitset::full)
        // : nodeset(0, nv - 1, bitset::empt)
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
    int size() const { return nodes.size(); }

    void add_edge(const int u, const int v);

    void add_edges(const int u, const adjacency_struct& N);

    void add_node(const int v);

    void add_clique(const adjacency_struct& C);

    void remove_edge(int u, int v);

    void remove_node(int v);

    void clear();
};

template <class adjacency_struct>
class graph : public basic_graph<adjacency_struct>
{
public:
    using basic_graph<adjacency_struct>::matrix;
    using basic_graph<adjacency_struct>::nodes;
    using basic_graph<adjacency_struct>::nodeset;

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
    graph()
        : basic_graph<adjacency_struct>()
    {
    }
    explicit graph(int nv)
        : basic_graph<adjacency_struct>(nv)
        , origmatrix(nv)
        , rep_of(nv)
        , partition(nv)
        , util_set(0, nv - 1, bitset::empt)
        , diff2(0, nv - 1, bitset::empt)
        , partu(0, nv - 1, bitset::empt)
        , partv(0, nv - 1, bitset::empt)
    {
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

    void add_edge(int u, int v);

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




template< class adjacency_struct >
struct clique_finder {
    const graph<adjacency_struct>& g;
    std::vector<adjacency_struct> cliques;
    std::vector<int> clique_sz;
    std::vector<adjacency_struct> candidates;
    std::vector<int> last_clique;
    int num_cliques;

    clique_finder(const graph<adjacency_struct>& g)
        : g(g)
        , num_cliques(1)
    {
        last_clique.resize(g.capacity());
        cliques.resize(g.capacity());
        clique_sz.resize(g.capacity());
        candidates.resize(g.capacity());
        for (auto& b : cliques)
            b.initialise(0, g.capacity(), bitset::empt);
        for (auto& b : candidates)
            b.initialise(0, g.capacity(), bitset::full);
    }

    // clear previously cached results
    void clear() { num_cliques = 0; }
    // initialize a new clique
    void new_clique();

    // initialize a new color
    void new_color();

    // insert v into the clq^th clique. assumes it fits
    void insert(int v, int clq);

    // insert v into the col^th color. assumes it fits. Puts vertices
    // added from candidates[i] into diff
    void insert_color(int v, int clq, bitset& diff);

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


template< class adjacency_struct >
struct degeneracy_finder {

    const basic_graph<adjacency_struct>& g;
    int d;
    std::vector<int> order;
    std::vector<int> degrees;
    std::vector<std::list<int>::iterator> iterators;
    std::vector<bool> ordered;
    std::vector<std::list<int>> buckets;

    degeneracy_finder(const basic_graph<adjacency_struct>& g)
        : g(g)
        , degrees(g.size())
        , iterators(g.size())
        , ordered(g.size())
    {
    }

    void degeneracy_ordering();
};




/** IMPLEMENTATION **/

template< class adjacency_struct >
void basic_graph<adjacency_struct>::add_edge(const int u, const int v)
{
    matrix[u].add(v);
    matrix[v].add(u);
}

template< class adjacency_struct >
void basic_graph<adjacency_struct>::add_edges(const int u, const adjacency_struct& N)
{
    matrix[u].union_with(N);
    for (auto v : N) {
        matrix[v].add(u);
    }
}

template< class adjacency_struct >
void basic_graph<adjacency_struct>::add_node(const int v)
{
    nodes.add(v);
    nodeset.add(v);
}

template< class adjacency_struct >
void basic_graph<adjacency_struct>::add_clique(const adjacency_struct& C)
{
    for (auto v : C) {
        add_node(v);
        matrix[v].union_with(C);
        matrix[v].remove(v);
    }
}

template< class adjacency_struct >
void basic_graph<adjacency_struct>::remove_edge(int u, int v)
{
    matrix[u].remove(v);
    matrix[v].remove(u);
}

template< class adjacency_struct >
void basic_graph<adjacency_struct>::remove_node(int v)
{
    nodes.remove(v);
    nodeset.fast_remove(v);
}

template< class adjacency_struct >
void basic_graph<adjacency_struct>::clear()
{
    for (auto v : nodes) {
        matrix[v].clear();
    }
    nodeset.clear();
    nodes.clear();
}




template< class adjacency_struct >
void graph<adjacency_struct>::add_dirty_edge(int u, int v)
{
    if (matrix[u].fast_contain(v))
        return;
    for (auto vp : partition[v]) {
        matrix[u].fast_add(vp);
        edge_trail.push_back({u, vp});
    }
}

template< class adjacency_struct >
void graph<adjacency_struct>::add_edge(int u, int v)
{
    matrix[u].add(v);
    matrix[v].add(u);
    origmatrix[u].add(v);
    origmatrix[v].add(u);
}

template< class adjacency_struct >
int graph<adjacency_struct>::merge(int u, int v)
{
    // if( v < u ) {
    //      auto w = u;
    //      u = v;
    //      v = w;
    // }

    //   util_set.copy(matrix[v]);
    // util_set.intersect_with(nodeset);
    //
    //   diff2.copy(matrix[u]);
    // diff2.intersect_with(nodeset);
    //
    //      if( diff2.size() > util_set.size() ) {
    //              auto w = u;
    //              u = v;
    //              v = w;
    //      }

    // util_set.clear();
    util_set.copy(matrix[v]);
    util_set.setminus_with(matrix[u]);

    // diff2.clear();
    diff2.copy(matrix[u]);
    diff2.setminus_with(matrix[v]);

    matrix[u].union_with(matrix[v]);
    for (auto w : util_set) {
        add_dirty_edge(w, u);
        add_dirty_edge(u, w);
    }
    for (auto w : diff2) {
        for (auto vp : partition[v])
            add_dirty_edge(w, vp);
    }

    // update rep_of for the partition that was absorbed
    for (auto vp : partition[v])
        rep_of[vp] = u;
    std::copy(
        begin(partition[v]), end(partition[v]), back_inserter(partition[u]));

    nodes.remove(v);
    nodeset.fast_remove(v);
    removed.push_back(v);

    return u;
}

// separate u and v. Just adds an edge, but it is reversible through
// checkpointing
template< class adjacency_struct >
void graph<adjacency_struct>::separate(int u, int v)
{
    if (matrix[u].contain(v)) {
        assert(matrix[v].contain(u));
        return;
    }
    if (cur_ckpt > 0) {
        add_dirty_edge(u, v);
        add_dirty_edge(v, u);
    } else {
        matrix[u].fast_add(v);
        matrix[v].fast_add(u);
    }
}

template< class adjacency_struct >
int graph<adjacency_struct>::contractPreprocess() {
    int num_contractions = 0;
    bool some_propagation = true;
    while (some_propagation) {
        some_propagation = false;
        for(auto u : nodes) {
            for(auto v : nodes) {
                if(u != v && !matrix[u].fast_contain(v)) {
                    util_set.copy(matrix[v]);
                    util_set.setminus_with(matrix[u]);
                    if(!util_set.intersect(nodeset)) {
                        // N(v) <= N(U)s
                        some_propagation = true;
                        ++num_contractions;
                        merge(u,v);
                    }
                }
            }
        }
    }
    return num_contractions;
}

template< class adjacency_struct >
int graph<adjacency_struct>::checkpoint()
{
    ++cur_ckpt;

    if (static_cast<size_t>(cur_ckpt) >= rep_of_trail.size()) {
        // make space for the copying part
        rep_of_trail.resize(cur_ckpt);
        partition_size_trail.resize(cur_ckpt);
        partition_size_trail[cur_ckpt - 1].resize(capacity());
    }
    // trailing
    edge_lim.push_back(edge_trail.size());

    // copying
    rep_of_trail[cur_ckpt - 1] = rep_of;
    for (auto v : nodes)
        partition_size_trail[cur_ckpt - 1][v] = partition[v].size();

    return cur_ckpt;
}

template< class adjacency_struct >
void graph<adjacency_struct>::restore(int ckpt)
{
    cur_ckpt = ckpt;

    int trailnewsize = edge_lim[ckpt];
    for (int i = trailnewsize; static_cast<unsigned>(i) != edge_trail.size();
         ++i) {
        auto e = edge_trail[i];
        matrix[e.first].fast_remove(e.second);
        matrix[e.second].fast_remove(e.first);
    }
    edge_trail.resize(trailnewsize);
    edge_lim.resize(ckpt);

    rep_of = rep_of_trail[cur_ckpt];

    while (!removed.empty()) {
        auto v = removed.back();
        assert(!nodes.contain(v));
        if (rep_of[v] != v)
            break;
        nodes.add(v);
        nodeset.fast_add(v);
        removed.pop_back();
    }

    for (auto v : nodes)
        partition[v].resize(partition_size_trail[cur_ckpt][v]);
}

template< class adjacency_struct >
void graph<adjacency_struct>::describe(std::ostream& os) const
{
    os << "# vertices = " << capacity() << std::endl;
		// for( auto v : nodes ) {
		// 	os << matrix[v].size() << std::endl;
		// }
}

template< class adjacency_struct >
void graph<adjacency_struct>::check_consistency() const
{
    std::cout << "checking graph consistency" << std::endl;
    for (int i = 0; i != capacity(); ++i)
        assert((!nodes.contain(i) || rep_of[i] == i)
            && (rep_of[i] != i || nodes.contain(i)));

    bitset bs(0, capacity() - 1, bitset::full);
    for (auto v : removed) {
        assert(!nodes.contain(v));
        bs.fast_remove(v);
    }
    assert(bs == nodeset);

    for (auto v : nodes) {
        assert(!partition[v].empty());
        assert(partition[v][0] == v);
        assert(rep_of[v] == v);
        bs.clear();
        for (auto u : partition[v]) {
            assert(rep_of[u] == v);
            bs.union_with(origmatrix[u]);
        }
        if (!bs.included(matrix[v])) {
            std::cout << "    bs[" << v << "] = " << bs << "\n";
            std::cout << "matrix[" << v << "] = " << matrix[v] << "\n";
        }
        assert(bs.included(matrix[v]));
    }

    for (auto v : nodes) {
        for (auto u : nodes) {
            if (u < v)
                continue;
            if (matrix[v].fast_contain(u)) {
                for (auto vp : partition[v])
                    if (!matrix[u].fast_contain(vp)) {
                        std::cout << "matrix[v] = " << matrix[v] << std::endl;
                        std::cout << "matrix[u] = " << matrix[u] << std::endl;
                        assert(matrix[u].fast_contain(vp));
                    }
                for (auto up : partition[u]) {
                    if (!matrix[v].fast_contain(up)) {
                        std::cout << "matrix[v] = " << matrix[v] << std::endl;
                        std::cout << "matrix[u] = " << matrix[u] << std::endl;
                        assert(matrix[v].fast_contain(up));
                    }
                }
            }
        }
    }
}


template< class adjacency_struct >
void clique_finder<adjacency_struct>::new_clique()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].clear();
    clique_sz[num_cliques] = 0;
    candidates[num_cliques].copy(g.nodeset);
    ++num_cliques;
}
// initialize a new color
template< class adjacency_struct >
void clique_finder<adjacency_struct>::new_color()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].clear();
    clique_sz[num_cliques] = 0;
    candidates[num_cliques].clear();
    ++num_cliques;
}
// insert v into the clq^th clique. assumes it fits
template< class adjacency_struct >
void clique_finder<adjacency_struct>::insert(int v, int clq)
{
    cliques[clq].fast_add(v);
    ++clique_sz[clq];
    candidates[clq].intersect_with(g.matrix[v]);
    last_clique[v] = clq;
}
// insert v into the col^th color. assumes it fits. Puts vertices
// added from candidates[i] into diff
template< class adjacency_struct >
void clique_finder<adjacency_struct>::insert_color(int v, int clq, bitset& diff)
{
    cliques[clq].fast_add(v);
    ++clique_sz[clq];
    diff.copy(g.matrix[v]);
    diff.setminus_with(candidates[clq]);
    candidates[clq].union_with(g.matrix[v]);
}
// heuristically find a set of cliques and return the size of the
// largest


template< class adjacency_struct >
void degeneracy_finder<adjacency_struct>::degeneracy_ordering()
{
     for (auto v : g.nodes) {
         auto vd = g.matrix[v].size();
         if (vd >= buckets.size())
             buckets.resize(vd+1);
         buckets[vd].push_front(v);
         degrees[v] = vd;
         iterators[v] = buckets[vd].begin();
         ordered[v] = false;
     }

     while(true) {
         size_t i{0};
         for (; i != buckets.size(); ++i)
             if (!buckets[i].empty())
                 break;
         if (i == buckets.size())
             break;
         d = std::max(d,static_cast<int>(i));
         auto v = buckets[i].back();
         order.push_back(v);
         buckets[i].pop_back();
         ordered[v] = true;
         for (auto u : g.matrix[v]) {
             if (ordered[u])
                 continue;
             auto &ud = degrees[u];
             buckets[ud].erase(iterators[u]);
             --ud;
             buckets[ud].push_front(u);
             iterators[u] = buckets[ud].begin();
         }
     }

     std::cout << "graph has degeneracy " << d << std::endl;
}




namespace detail
{
    template <class adjacency_struct> struct brelaz_state {
        clique_finder<adjacency_struct> cf;

        // degrees are valid only if the corresponding bit is not set in dirty
        mutable std::vector<int> degrees;
        mutable bitset dirty;

        // current set of uncolored nodes
        bitset nodes;

        // saturation
        std::vector<int> saturation;

        // temp
        mutable bitset util_set;

        void remove(int v)
        {
            nodes.fast_remove(v);
            dirty.union_with(cf.g.matrix[v]);
        }
        int degree(int v) const
        {
            if (dirty.fast_contain(v)) {
                util_set.copy(cf.g.matrix[v]);
                util_set.intersect_with(nodes);
                degrees[v] = util_set.size();
                dirty.fast_remove(v);
            }
            return degrees[v];
        }

        brelaz_state(const graph<adjacency_struct>& g)
            : cf(g)
            , degrees(cf.g.capacity(), 0)
            , dirty(0, cf.g.capacity(), bitset::full)
            , nodes(cf.g.nodeset)
            , saturation(cf.g.capacity(), 0)
            , util_set(0, cf.g.capacity(), bitset::empt)
        {
        }
    };

    template <class adjacency_struct> struct saturation_gt {
        const brelaz_state<adjacency_struct>& bs;

        bool operator()(int u, int v)
        {
            int satlt = bs.saturation[u] - bs.saturation[v];
            if (satlt > 0)
                return true;
            if (satlt < 0)
                return false;
            return bs.degree(u) > bs.degree(v);
        }
    };
} // namespace detail

template< class adjacency_struct >
std::vector<int> brelaz_color(const graph<adjacency_struct>& g)
{
    detail::brelaz_state<adjacency_struct> state{g};
    auto& cf = state.cf;
    cf.clear();
    std::vector<int> solution(cf.g.capacity(), -1);
    Heap<detail::saturation_gt<adjacency_struct>> sheap(
        detail::saturation_gt<adjacency_struct>{state});

    for (auto v : cf.g.nodeset)
        sheap.insert(v);

    while (!sheap.empty()) {
        int v = sheap.removeMin();
        state.remove(v); // O(N/64) on bitsets, O(|N(v)|) on lists

        bool found{false};
        for (int i = 0; i != cf.num_cliques; ++i) {
            if (!cf.candidates[i].fast_contain(v)) {
                found = true;
                state.util_set.clear();
                cf.insert_color(v, i, state.util_set);
								
								std::cout << "color " << v << " with " << i << " (" << g.matrix[v].size() << "):";
								for( auto u : g.matrix[v] ) {
									if( state.nodes.contain(u) ) {
										std::cout << " " << u;
									}
								}
								std::cout << std::endl ;
								for (int j = 0; j != cf.num_cliques; ++j) {
									std::cout << "v[" << j << "] = " ;
									for( auto u : cf.candidates[j] )
										std::cout << " " << u;
									std::cout << std::endl;
								}
								std::cout << std::endl ;
								
								
                state.util_set.intersect_with(state.nodes);
                for (auto u : state.util_set) {
                    ++state.saturation[u];
                    sheap.update(u);
                }
                solution[v] = i;
                break;
            }
        }
        if (!found) {
            cf.new_color();
            cf.insert_color(v, cf.num_cliques - 1, state.util_set);
            solution[v] = cf.num_cliques - 1;
        }
    }

    return solution;
}

} // namespace gc

#endif
