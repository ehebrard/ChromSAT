#ifndef __CG_BI_GRAPH_HH
#define __CG_BI_GRAPH_HH

#include <vector>

#include "bitset.hpp"


namespace gc
{


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

    Q.clear();
    q = 0;

    vmap.clear();
    vmap.resize(g.capacity(), -1);

    matrix.resize(g.capacity());

    original.clear();

    first.reinitialise(0, g.capacity() - 1, bitset::empt);

    std::sort(begin(c1), end(c1));
    std::sort(begin(c2), end(c2));

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
								
								// std::cout << " [" << original[prev] << ".." << v << "[";
								
                for (int w = prev; w < next; ++w) {
                    matrix[i].push_back(w);
                    matrix[w].push_back(i);
                }
                prev = next + 1;
            }
        }
				
        for (int w = prev; w < original.size(); ++w) {
            matrix[i].push_back(w);
            matrix[w].push_back(i);
        }
				
				// std::cout << std::endl;
				
    }
 
}

} // namespace gc

#endif
