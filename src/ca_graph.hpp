#ifndef __CG_CA_GRAPH_HH
#define __CG_CA_GRAPH_HH

#include <limits>

#include "dsatur.hpp"
#include "intstack.hpp"
#include "statistics.hpp"

#include <vector>

namespace gc
{

// template <class container> bool is_sorted(container& c)
// {
//     auto prev{c[0] - 1};
//     for (auto x : c) {
//         if (x < prev)
//             return false;
//         prev = x;
//     }
//
//     return true;
// }

template <class graph_struct> struct cliquer {
    const graph_struct& g;
    std::vector<std::vector<int>> cliques;
    std::vector<gc::bitset> candidates;
    std::vector<int> last_clique;
    size_t num_cliques;
    size_t maxcliquesize;

    cliquer(const graph_struct& ig)
        : g(ig)
        , num_cliques(0)
        , maxcliquesize(0)
    {
        last_clique.resize(g.capacity());
    }

    // clear previously cached results
    void clear()
    {
        num_cliques = 0;
        maxcliquesize = 0;
    }
    // initialize a new clique
    void new_clique();

    // insert v into the clq^th clique. assumes it fits
    void insert(int v, int clq);

    // heuristically find a set of cliques and return the size of the
    // largest
    template <class ForwardIterator>
    int find_cliques(ForwardIterator beg_ordering, ForwardIterator end_ordering, const int limit=0xfffffff)
    {
        clear();
        for (auto it{beg_ordering}; it != end_ordering; ++it) {

            // std::cout << *it ;

            bool found{false};
            for (int i = 0; i != num_cliques; ++i)
                if (candidates[i].fast_contain(*it)) {
                    found = true;
                    insert(*it, i);
                    // std::cout << ":" << i << " ";
                }
            if (!found && num_cliques < limit) {
                new_clique();
                insert(*it, num_cliques - 1);
                // std::cout << ":" << (num_cliques - 1) << "* ";
            }
        }

        // std::cout << "| ";

        for (auto it{beg_ordering}; it != end_ordering; ++it) {
            for (int i = last_clique[*it] + 1; i < num_cliques; ++i)
                if (candidates[i].fast_contain(*it)) {
                    insert(*it, i);
                    // std::cout << *it << ":" << (num_cliques - 1) << " ";
                }
        }

        // std::cout << " -> " << maxcliquesize << std::endl;

        return maxcliquesize;
    }

    // heuristically find a set of cliques and return the size of the
    // largest
    template <typename ForwardIterator, typename set>
    int find_cliques_in(
        ForwardIterator beg_ordering, ForwardIterator end_ordering, set nodes)
    {
        clear();
        for (auto it{beg_ordering}; it != end_ordering; ++it) {
            // std::cout << *it ;

            if (!nodes.contain(*it))
                continue;

            bool found{false};
            for (int i = 0; i != num_cliques; ++i)
                if (candidates[i].fast_contain(*it)) {
                    found = true;
                    insert(*it, i);
                    // std::cout << ":" << i << " ";
                }
            if (!found) {
                new_clique();
                insert(*it, num_cliques - 1);
                // std::cout << ":" << (num_cliques - 1) << "* ";
            }
        }

        // std::cout << "| ";

        for (auto it{beg_ordering}; it != end_ordering; ++it) {

            if (!nodes.contain(*it))
                continue;

            for (int i = last_clique[*it] + 1; i < num_cliques; ++i)
                if (candidates[i].fast_contain(*it)) {
                    insert(*it, i);
                    // std::cout << *it << ":" << (num_cliques - 1) << " ";
                }
        }

        // std::cout << " -> " << maxcliquesize << std::endl;

        return maxcliquesize;
    }
};

template <class graph_struct> void cliquer<graph_struct>::new_clique()
{
    if (cliques.size() == num_cliques) {
        cliques.resize(num_cliques + 1);
        candidates.resize(num_cliques + 1);
        candidates.back().initialise(0, g.capacity(), bitset::empt);
    }

    cliques[num_cliques].clear();
    candidates[num_cliques].copy(g.nodeset);
    ++num_cliques;
}

// insert v into the clq^th clique. assumes it fits
template <class graph_struct> void cliquer<graph_struct>::insert(int v, int clq)
{
    cliques[clq].push_back(v);
    candidates[clq].intersect_with(begin(g.matrix[v]), end(g.matrix[v]));
    last_clique[v] = clq;
    maxcliquesize = std::max(maxcliquesize, cliques[clq].size());
}

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
        } else if (*it == value) {
            return it;
        } else
            count = step;
    }

    return first;
}

class bi_graph
{

public:
    const int INFTY = std::numeric_limits<int>::max();

private:
    // static const int T = -1;

    std::vector<int> vmap;
    bitset first;

    std::vector<int> dist;

    std::vector<int> Q;
    int q;

public:
    int N;
    int T;
    int I;

    std::vector<int> original;
    std::vector<int> matching;

    // intstack nodes;

    std::vector<std::vector<int>> matrix;

    bi_graph() {}

    template <class graph_struct, class clique_struct>
    void get_from_cliques(
        graph_struct& g, clique_struct& c1, clique_struct& c2);

    bool dfs(const int u);
    bool bfs();
    int hopcroftKarp();

    template <class graph_struct, class clique_struct>
    int get_bound(graph_struct& g, clique_struct& c1, clique_struct& c2);

    std::ostream& describe(std::ostream& os) const;
};

template <class graph_struct, class clique_struct>
int bi_graph::get_bound(graph_struct& g, clique_struct& c1, clique_struct& c2)
{
    get_from_cliques(g, c1, c2);
    return c1.size() + c2.size() - I - hopcroftKarp();
}

template <class graph_struct, class clique_struct>
void bi_graph::get_from_cliques(
    graph_struct& g, clique_struct& c1, clique_struct& c2)
{
    // std::cout << g << std::endl;

    Q.clear();
    q = 0;

    vmap.clear();
    vmap.resize(g.capacity(), -1);

    // nodes.clear();
    // nodes.reserve(c1.size() + c2.size());

    matrix.resize(g.capacity());

    original.clear();

    first.reinitialise(0, g.capacity() - 1, bitset::empt);

    std::sort(begin(c1), end(c1));
    std::sort(begin(c2), end(c2));

    // for (auto u : c1) {
    //     std::cout << " " << u;
    // }
    // std::cout << std::endl;
    // for (auto u : c2) {
    //     std::cout << " " << u;
    // }
    // std::cout << std::endl;
    //
    // std::cout << c1.size() << " " << c2.size() << std::endl;

    N = 0;
    I = 0;
    int i{0}, j{0};
    while (i < c1.size() or j < c2.size()) {

        // std::cout << i << " " << j;
        if (i < c1.size() and j < c2.size()) {
            if (c1[i] < c2[j]) {
                // std::cout << " -> " << c1[i];
                original.push_back(c1[i]);
                first.add(c1[i++]);
                ++N;
            } else if (c1[i] > c2[j]) {
                // std::cout << " -> " << c2[j];
                original.push_back(c2[j++]);
            } else {
                ++I;
                ++i;
                ++j;
            }
        } else if (i == c1.size()) {
            while (j < c2.size()) {
                // std::cout << " -> " << c2[j];
                original.push_back(c2[j++]);
            }
        } else {
            while (i < c1.size()) {
                // std::cout << " -> " << c1[i];
                original.push_back(c1[i]);
                first.add(c1[i++]);
                ++N;
            }
        }

        // std::cout << std::endl;
    }
    // for (auto u : original) {
    //     std::cout << " " << u;
    // }
    // std::cout << std::endl;

    std::sort(begin(original), end(original), [&](int x, int y) {
        return first.contain(x) > first.contain(y)
            or (first.contain(x) == first.contain(y) and x < y);
    });

    i = 0;
    for (auto u : original) {
        // std::cout << " " << u;
        vmap[u] = i++;

        // if (i == N)
        //     std::cout << " |";
    }
    // std::cout << std::endl;

    T = original.size();
    // nodes.reserve(original.size());
    // nodes.fill();

    // T = nodes.size();

    matching.clear();
    matching.resize(T + 1, T);

    dist.clear();
    dist.resize(T + 1, INFTY);

    // for (auto u : c1) {
    //     std::cout << " " << u;
    //     if (vmap[u] == -1) {
    //         original.push_back(u);
    //         vmap[u] = nodes.size();
    //         nodes.push(vmap[u]);
    //         matrix[vmap[u]].clear();
    //     }
    // }
    // std::cout << std::endl;
    //
    // int N{static_cast<int>(nodes.size())};
    //
    // for (auto u : c2) {
    //     std::cout << " " << u;
    //     if (vmap[u] == -1) {
    //         original.push_back(u);
    //         vmap[u] = nodes.size();
    //         nodes.push(vmap[u]);
    //         matrix[vmap[u]].clear();
    //     }
    // }
    // std::cout << std::endl;

    for (i = 0; i < T; ++i) {
        matrix[i].clear();
    }

    for (i = 0; i < N; ++i) {
        auto u{original[i]};
        auto prev{N};
        auto p{-1};
        for (auto v : g.matrix[u]) {
            assert(p < v);
            p = v;
            if (vmap[v] >= N) {
                auto next{vmap[v]};
                for (int w = prev; w < next; ++w) {
                    matrix[i].push_back(w);
                    matrix[w].push_back(i);
                }
                prev = next + 1;
            }
        }
    }
    // for (auto u : ) {
    //     int prev{N};
    //     for (auto v : g.matrix[u]) {
    //         if (u < v and vmap[v]) {
    //             int next{vmap[v]};
    //             for (int w = prev; w < next; ++w) {
    //                 matrix[vmap[u]].push_back(w);
    // 										matrix[w].push_back(vmap[u]);
    // 								}
    //             prev = next + 1;
    //         }
    //     }
    // }

    // this->describe(std::cout);

    // exit(1);
    // std::cout << *this << std::endl;
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

    void search(statistics& stat, options& options);

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
