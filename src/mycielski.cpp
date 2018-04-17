#include "mycielski.hpp"
#include "utils.hpp"

#include <minicsp/mtl/Heap.h>

#include <algorithm>

#include <iomanip>

using namespace minicsp;


namespace gc
{

mycielskan_subgraph_finder::mycielskan_subgraph_finder(
    const graph& g, const clique_finder& cf, const bool prune)
    : prune(prune)
    , g(g)
    , cf(cf)
    , explanation_subgraph(g.capacity())
    , subgraph(g.capacity())
    , ith_node(0)
{
    non_neighbors.initialise(0, g.capacity() + 1, bitset::empt);
    neighbors_Sv.initialise(0, g.capacity() + 1, bitset::empt);
    neighbors_w.initialise(0, g.capacity() + 1, bitset::empt);
    candidates.initialise(0, g.capacity() + 1, bitset::empt);
    pruning.initialise(0, g.capacity() + 1, bitset::empt);
    real_pruning.initialise(0, g.capacity() + 1, bitset::empt);

    explanation_clique = -1;
}

// here we assume that the current upper bound is k+1 and we just failed extending the current subgraph to order k+1
// ith_node is the node after which the subgraph at which candidates would become empty, we then go through the nodes from ith_node to the last
Clause* mycielskan_subgraph_finder::do_prune(Solver& s, const std::vector<std::vector<Var>>& vars)
{
    pruning.copy(candidates);
    candidates.fill();

    if (another_myciel_layer(ith_node)
        == static_cast<int>(subgraph.nodes.size())) {
        // first compute the potential pruning
        new_pruning.clear();
        for (auto u : pruning) {
            real_pruning.copy(candidates);
            real_pruning.setminus_with(g.origmatrix[u]);
            for (auto v : real_pruning) {
                if (s.value(Lit(vars[u][v])) == l_Undef)
                    new_pruning.push_back(edge{u, v});
            }
        }

#ifdef _DEBUG_MYCIEL
        std::cout << " prune " << pruning << " x " << candidates << " ("
                  << new_pruning.size() << ")\n";
#endif

        if (new_pruning.size() > 0) {
            // base graph same for every pruning
            reason.clear();
            for (auto u : subgraph.nodes) {
                for (auto v : subgraph.matrix[u]) {
                    if (!g.origmatrix[u].fast_contain(v)) {
                        reason.push(Lit(vars[u][v]));
                    }
                }
            }
            auto size_reason = reason.size();

            for (auto e : new_pruning) {
                // the edge gives us w1 and w2
                auto w1{e.first};
                auto w2{e.second};

                new_edges.clear();
                u_layer.clear();

#ifdef _DEBUG_MYCIEL
                std::cout << " select U layer for " << w1 << " in [" << 0
                          << ".." << ith_node << "]\n";
                int j = 0;
                for (int i = 0; i < ith_node; ++i) {
                    std::cout << " | " << subgraph.nodes[i] << ":";
                    while (j < endS[i]) {
                        std::cout << " " << extra[j++];
                    }
                }
                std::cout << std::endl;
#endif

                select_middle_layer(w1, 0, ith_node, u_layer, new_edges);

                for (auto ee : new_edges) {
                    if (!g.origmatrix[ee.first].fast_contain(ee.second)) {
                        reason.push(Lit(vars[ee.first][ee.second]));
                    }
                }

#ifdef _DEBUG_MYCIEL
                std::cout << " select U layer for " << w2 << " in [" << ith_node
                          << ".." << subgraph.nodes.size() << "]\n";
                for (int i = ith_node; i < subgraph.nodes.size(); ++i) {
                    std::cout << "| " << subgraph.nodes[i] << ":";
                    while (j < endS[i]) {
                        std::cout << " " << extra[j++];
                    }
                }
                std::cout << std::endl;
#endif

                new_edges.clear();
                select_middle_layer(
                    w2, ith_node, subgraph.nodes.size(), u_layer, new_edges);

                for (auto ee : new_edges) {
                    if (!g.origmatrix[ee.first].fast_contain(ee.second)) {
                        reason.push(Lit(vars[ee.first][ee.second]));
                    }
                }

                DO_OR_RETURN(s.enqueueFill(~Lit(vars[w1][w2]), reason));

                // std::cout << "shrink from " <<
                // reason.size() << " to " <<
                // size_reason << std::endl;
                reason.shrink_(reason.size() - size_reason);
            }
        }
    }
#ifdef _DEBUG_MYCIEL
    else
        std::cout << "\n";
#endif

    return minicsp::NO_REASON;
}

int mycielskan_subgraph_finder::another_myciel_layer(const int start)
{
    int ith = start;
    while (ith < static_cast<int>(subgraph.nodes.size())) {
        auto v{subgraph.nodes[ith++]};

        extra.push_back(v);
        neighbors_Sv.copy(g.matrix[v]);

#ifdef _DEBUG_MYCIEL
        std::cout << " - " << v << ":";
#endif
				
        // non_neighbors.copy(g.nodeset);
        // non_neighbors.setminus_with(g.matrix[v]);
        // non_neighbors.fast_remove(v);
        // for (auto u : non_neighbors) {
				for(auto u : g.nodes) {
						if (u==v)
								continue;

            if (g.rep_of[u] != u)
                continue; // use only representatives ?

            if (g.matrix[u].includes(subgraph.matrix[v])) {
                extra.push_back(u);
                neighbors_Sv.union_with(g.matrix[u]);

#ifdef _DEBUG_MYCIEL
                std::cout << " " << u;
#endif
            }
        }

        if (!candidates.intersect(neighbors_Sv)) {
#ifdef _DEBUG_MYCIEL
            std::cout << " -> {} at node " << (ith) << "/"
                      << subgraph.nodes.size() << "\n";
#endif
            extra.resize(endS.back());
            return ith - 1;
        }

        // stop early when there is no candidate for w
        candidates.intersect_with(neighbors_Sv);

#ifdef _DEBUG_MYCIEL
        std::cout << " -> " << candidates << std::endl;
#endif

        endS.push_back(extra.size());
    }

    return ith;
}


void mycielskan_subgraph_finder::select_middle_layer(const int w, const int beg_node, const int end_node, std::vector<int>& nodes, std::vector<edge>& edges) {

    // put all potential extra nodes into a bitset so that we can intersect with
    // N(u)
    neighbors_w.clear();
    for (auto v : extra) {
        neighbors_w.fast_add(v);
    }
    neighbors_w.intersect_with(g.matrix[w]);

    // now for every v, select any element of Sv that is also a neighbor of u
    auto j{(beg_node ? endS[beg_node - 1] : 0)}; // index in extra

    for (auto i = beg_node; i < end_node; ++i) {

#ifdef _DEBUG_MYCIEL
        std::cout << " " << subgraph.nodes[i];
#endif

        do {
            assert(j < endS[i]);
            auto u{extra[j]};

            if (neighbors_w.fast_contain(u)) {
                edges.push_back(edge{u, w});
                if (!subgraph.nodeset.fast_contain(u)) {
                    nodes.push_back(u);
                    for (auto v : subgraph.matrix[subgraph.nodes[i]]) {
                        edges.push_back(edge{u, v});
                    }

#ifdef _DEBUG_MYCIEL
                    std::cout << "--" << u << "--" << w << std::endl;
#endif
                }
#ifdef _DEBUG_MYCIEL
                else
                    std::cout << "--" << w << std::endl;
#endif
                break;
            }
            ++j;
        } while (true);

        j = endS[i];
    }
}


int mycielskan_subgraph_finder::extends( const bitset& G )
{
    int iter = 0;

    subgraph.clear();
    subgraph.add_clique(G);

    while (subgraph.nodes.size() < g.nodes.size()) {
        // new ietration, clear the struct for the extra nodes and the potential
        // w's
        extra.clear(); // extra contains the union of the Sv's
        endS.clear(); // where do the Sv's en in the vector extra
        candidates.fill(); // potential candidates for w

#ifdef _DEBUG_MYCIEL
        std::cout << "\nextends " << subgraph.nodeset << std::endl;
#endif

        ith_node = another_myciel_layer(0);

        assert(ith_node == endS.size());

        if (ith_node < subgraph.nodes.size())
            return iter;

        // select any (?) w
        auto w = g.rep_of[candidates.min()];

        new_edges.clear();
        u_layer.clear();

        select_middle_layer(w, 0, subgraph.nodes.size(), u_layer, new_edges);

        for (auto u : u_layer) {
            subgraph.add_node(u);
        }
        subgraph.add_node(w);
        for (auto e : new_edges) {
            subgraph.add_edge(e.first, e.second);
        }

        ++iter;
    }

    return iter;
}

int mycielskan_subgraph_finder::improve_cliques_larger_than(const int size)
{
    explanation_clique = -1;
    auto lb{0};
    for (auto cl = 0; cl < cf.num_cliques; ++cl) {
                                if (cf.clique_sz[cl] >= size) {
                        auto mycielski_lb = cf.clique_sz[cl] + extends(cf.cliques[cl]);
                        if (mycielski_lb > lb) {
                            lb = mycielski_lb;
                        }
                                }
    }
    return lb;
}

int mycielskan_subgraph_finder::full_myciel(const int curlb, const int ub,
    minicsp::Solver& s, const std::vector<std::vector<minicsp::Var>>& vars)
{
    explanation_clique = -1;
    auto lb{curlb};
    for (auto cl = 0; cl < cf.num_cliques; ++cl) {
        auto mycielski_lb = cf.clique_sz[cl] + extends(cf.cliques[cl]);
        if (mycielski_lb > lb) {
            lb = mycielski_lb;
            explanation_subgraph = subgraph;
            explanation_clique = cl;

            if (prune && ub - lb == 1) {
                do_prune(s, vars);
            }
        }
    }
    return lb;
}

int mycielskan_subgraph_finder::improve_cliques_larger_than(const int size,
    const int curlb, const int ub, minicsp::Solver& s,
    const std::vector<std::vector<minicsp::Var>>& vars)
{
    explanation_clique = -1;
    auto lb{curlb};
    for (auto cl = 0; cl < cf.num_cliques; ++cl) {
        if (cf.clique_sz[cl] >= size) {
            auto niters{extends(cf.cliques[cl])};
            auto mycielski_lb = cf.clique_sz[cl] + niters;
            if (mycielski_lb > lb) {
                lb = mycielski_lb;
                explanation_subgraph = subgraph;
                explanation_clique = cl;

                if (prune && ub - lb == 1) {
                    do_prune(s, vars);
                }
            }
        }
    }

    return lb;
}

int mycielskan_subgraph_finder::improve_greedy(const int size, const int curlb,
    const int ub, minicsp::Solver& s,
    const std::vector<std::vector<minicsp::Var>>& vars)
{
    explanation_clique = -1;
    auto lb{curlb};
    for (auto cl = 0; cl < cf.num_cliques; ++cl) {
        if (cf.clique_sz[cl] >= lb) {
            auto mycielski_lb = cf.clique_sz[cl] + extends(cf.cliques[cl]);
            if (mycielski_lb > lb) {
                lb = mycielski_lb;
                explanation_subgraph = subgraph;
                explanation_clique = cl;

                if (prune && ub - lb == 1) {
                    do_prune(s, vars);
                }
            }
        }
    }
    return lb;
}


}
