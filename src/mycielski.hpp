#ifndef __CG_MYCIELSKI_HH
#define __CG_MYCIELSKI_HH

#include "bitset.hpp"
#include "intstack.hpp"
#include "graph.hpp"
#include "minicsp/core/solver.hpp"

#include <algorithm>
#include <list>
#include <vector>


// #define _DEBUG_MYCIEL

namespace gc
{

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

}

#endif
