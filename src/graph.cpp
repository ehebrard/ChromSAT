#include "graph.hpp"
#include "utils.hpp"

#include <minicsp/mtl/Heap.h>

#include <algorithm>

#include <iomanip>

using namespace minicsp;


namespace gc
{

void graph::add_dirty_edge(int u, int v)
{
    if (matrix[u].fast_contain(v))
        return;
    for (auto vp : partition[v]) {
        matrix[u].fast_add(vp);
        edge_trail.push_back({u, vp});
    }
}

int graph::merge(int u, int v)
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

void graph::separate(int u, int v)
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

void contractPreprocess() {
		int num_contractions = 0; 
		bool some_propagation = true;
		while(some_propagation) {
				some_propagation = false;
				for(auto u : nodes) {
						for(auto v : nodes) {
								if(u != v && !matrix[u].fast_contain(v)) {
										neighborhood.copy(matrix[v]);
										neighborhood.setminus_with(matrix[u]);
										if(!neighborhood.intersect(nodeset)) {
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

int graph::checkpoint()
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

void graph::restore(int ckpt)
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

void graph::describe(std::ostream& os) const
{
    os << "# vertices = " << capacity() << std::endl;
}

void graph::check_consistency() const
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

clique_finder::clique_finder(const graph& g)
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

void clique_finder::clear() { num_cliques = 0; }

void clique_finder::new_clique()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].clear();
    clique_sz[num_cliques] = 0;
    candidates[num_cliques].copy(g.nodeset);
    ++num_cliques;
}

void clique_finder::new_color()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].clear();
    clique_sz[num_cliques] = 0;
    candidates[num_cliques].clear();
    ++num_cliques;
}

void clique_finder::insert(int v, int clq)
{
    cliques[clq].fast_add(v);
    ++clique_sz[clq];
    candidates[clq].intersect_with(g.matrix[v]);
    last_clique[v] = clq;
}

void clique_finder::insert_color(int v, int clq, bitset& diff)
{
    cliques[clq].fast_add(v);
    ++clique_sz[clq];
    diff.copy(g.matrix[v]);
    diff.setminus_with(candidates[clq]);
    candidates[clq].union_with(g.matrix[v]);
}

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

        non_neighbors.copy(g.nodeset);
        non_neighbors.setminus_with(g.matrix[v]);
        non_neighbors.fast_remove(v);

        // Sv contains v
        extra.push_back(v);
        neighbors_Sv.copy(g.matrix[v]);

#ifdef _DEBUG_MYCIEL
        std::cout << " - " << v << ":";
#endif

        for (auto u : non_neighbors) {

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


neighbors_wrapper::neighbors_wrapper(const graph& g)
    : g(g)
{
    by_degree.resize(g.capacity());
    neighbors.resize(g.capacity());
    degree.resize(g.capacity());

    for (auto v : g.nodes) {
        by_degree[v].reserve(g.capacity());
        neighbors[v].reserve(g.capacity());
    }

    for (auto v : g.nodes) {
        by_degree[v].reserve(g.capacity());
        neighbors[v].reserve(g.capacity());
    }

    buffer.initialise(0, g.capacity(), bitset::empt);
}

void neighbors_wrapper::synchronize()
{
    // // the nodes betweem g.nodes.size() and size have been removed
    // while( size > g.nodes.size() ) {
    //              auto v = g.nodes[--size];
    //              for( auto u : neighbors[v] ) {
    //                              neighbors[u].remove( v );
    //              }
    // }
    // // the nodes betweem size and g.nodes.size() have been added
    // while( size < g.nodes.size() ) {
    //              auto v = g.nodes[size++];
    //              for( auto u : neighbors[v] ) {
    //                              neighbors[u].push( v );
    //              }
    // }

    for (auto v : g.nodes) {
        neighbors[v].clear();
        if (!g.matrix[v].empty()) {
            buffer.copy(g.matrix[v]);
            buffer.intersect_with(g.nodeset);
            for (auto u : buffer)
                neighbors[v].push(u);
        }
    }

    check_consistency();
}

void neighbors_wrapper::get_degeneracy_order(std::vector<int>& order)
{
    synchronize();
    for (auto v : g.nodes) {
        degree[v] = neighbors[v].size();
        by_degree[degree[v]].push(v);
    }

    for (auto i = g.nodes.size(); i-- > 0;) {
        for (auto& vertices : by_degree) {
            if (!vertices.empty()) {
                auto v = vertices.back();
                vertices.pop_back();

                for (auto ni = 0; ni < degree[v]; ++ni) {
                    auto u = neighbors[v][ni];
                    by_degree[degree[u]].remove(u);
                    neighbors[u].move_up(v, --degree[u]);
                    by_degree[degree[u]].push(u);
                }

                order.push_back(v);
                break;
            }
        }
    }
}

void neighbors_wrapper::check_consistency() const
{

    bitset buffer(0, g.capacity(), bitset::empt);

    // assert( size == g.nodes.size() );

    for (auto v : g.nodes) {
        buffer.copy(g.matrix[v]);
        buffer.intersect_with(g.nodeset);

        assert(buffer.size() == neighbors[v].size());

        for (auto u : neighbors[v]) {
            assert(buffer.contain(u));
        }
    }
}

// degeneracy_finder::degeneracy_finder(const graph& g)
//     : g(g)
// {
//
//      degrees.resize(g.capacity());
//      iterators.resize(g.capacity());
//      ordered.initialise(0, g.capacity()-1, bitset::empt);
// }
//
//
//
// void degeneracy_finder::get_degeneracy_order( std::vector< int >& order  ) {
//              for (auto v : g.nodes) {
//                  // auto vd = g.neighbor(v).size();
//                  if (vd >= buckets.size())
//                      buckets.resize(vd + 1);
//                  buckets[vd].push_front(v);
//                  degrees[v] = vd;
//                  iterators[v] = buckets[vd].begin();
//              }
//
//              ordered.clear();
//     while (true) {
//         size_t i{0};
//         for (; i != buckets.size(); ++i)
//             if (!buckets[i].empty())
//                 break;
//         if (i == buckets.size())
//             break;
//         auto v = buckets[i].back();
//         order.push_back(v);
//         buckets[i].pop_back();
//         ordered.fast_add(v);
//         for (auto u : neighbor(v)) {
//             if (ordered.fast_contain(u))
//                 continue;
//             auto& ud = degrees[u];
//             buckets[ud].erase(iterators[u]);
//             --ud;
//             buckets[ud].push_front(u);
//             iterators[u] = buckets[ud].begin();
//         }
//     }
// }

namespace detail
{
    struct brelaz_state {
        clique_finder cf;

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

        brelaz_state(const graph& g)
            : cf(g)
            , degrees(cf.g.capacity(), 0)
            , dirty(0, cf.g.capacity(), bitset::full)
            , nodes(cf.g.nodeset)
            , saturation(cf.g.capacity(), 0)
            , util_set(0, cf.g.capacity(), bitset::empt)
        {
        }
    };

    struct saturation_gt {
        const brelaz_state& bs;

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

std::vector<int> brelaz_color(const graph& g)
{
    detail::brelaz_state state{g};
    auto& cf = state.cf;
    cf.clear();
    std::vector<int> solution(cf.g.capacity(), -1);
    auto sheap = Heap<detail::saturation_gt>(detail::saturation_gt{state});

    for (auto v : cf.g.nodeset)
        sheap.insert(v);

    while (!sheap.empty()) {
        int v = sheap.removeMin();
        state.remove(v);

        bool found{false};
        for (int i = 0; i != cf.num_cliques; ++i) {
            if (!cf.candidates[i].fast_contain(v)) {
                found = true;
                state.util_set.clear();
                cf.insert_color(v, i, state.util_set);
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
