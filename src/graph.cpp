#include "graph.hpp"
#include "utils.hpp"

#include <algorithm>
#include <iomanip>

namespace gc
{

void graph::add_dirty_edge(int u, int v)
{
    if (matrix[u].fast_contain(v))
        return;
    if (cur_ckpt > 0) {
        if (!dirty[cur_ckpt].fast_contain(u)) {
            dirty[cur_ckpt].fast_add(u);
            diffs[cur_ckpt][u].clear();
        }
        diffs[cur_ckpt][u].fast_add(v);
    }
    matrix[u].fast_add(v);
}

int graph::merge(int u, int v)
{
    util_set.clear();
    util_set.copy(matrix[v]);
    util_set.setminus_with(matrix[u]);
    matrix[u].union_with(matrix[v]);

    // update rep_of for the partition that was absorbed
    for (auto vp : partition[v])
        rep_of[vp] = u;
    std::copy(
        begin(partition[v]), end(partition[v]), back_inserter(partition[u]));

    nodes.remove(v);
    nodeset.fast_remove(v);
    removed.push_back(v);

    if (cur_ckpt > 0) {
        if (!dirty[cur_ckpt].fast_contain(u)) {
            diffs[cur_ckpt][u].copy(util_set);
            dirty[cur_ckpt].fast_add(u);
        } else {
            diffs[cur_ckpt][u].union_with(util_set);
        }
    }
    for(auto w : util_set) {
        add_dirty_edge(w, u);
        add_dirty_edge(u, w);
    }

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

int graph::checkpoint()
{
    ++cur_ckpt;
    if (static_cast<size_t>(cur_ckpt) >= diffs.size()) {
        // trailing part
        diffs.resize(cur_ckpt + 1);
        dirty.resize(cur_ckpt + 1);
        diffs[cur_ckpt].resize(capacity());
        for (auto& bs : diffs[cur_ckpt])
            bs.initialise(0, capacity() + 1, bitset::empt);
        dirty[cur_ckpt].initialise(0, capacity() + 1, bitset::empt);

        // copying part
        rep_of_trail.resize(cur_ckpt);
        rep_of_trail[cur_ckpt-1] = rep_of;
        partition_size_trail.resize(cur_ckpt);
        partition_size_trail[cur_ckpt-1].resize(capacity());
        for (auto v : nodes)
            partition_size_trail[cur_ckpt-1][v] = partition[v].size();
    } else {
        // trailing
        for (auto& bs : diffs[cur_ckpt])
            bs.clear();
        dirty[cur_ckpt].clear();

        // copying
        rep_of_trail[cur_ckpt-1] = rep_of;
        for (auto v : nodes)
            partition_size_trail[cur_ckpt-1][v] = partition[v].size();
    }
    return cur_ckpt;
}

void graph::restore(int ckpt)
{
    for (int i = cur_ckpt; i > ckpt; --i) {
        for (auto v : dirty[i]) {
            matrix[v].setminus_with(diffs[i][v]);
        }
    }
    cur_ckpt = ckpt;
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

    bitset bs(0, capacity(), bitset::full);
    for (auto v : removed) {
        assert(!nodes.contain(v));
        bs.fast_remove(v);
    }
    assert(bs == nodeset);

    for (auto v : nodes) {
        assert(!partition[v].empty());
        assert(partition[v][0] == v);
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
}

clique_finder::clique_finder(const graph& g)
    : g(g)
    , num_cliques(1)
{
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

void clique_finder::insert(int v, int clq)
{
    cliques[clq].fast_add(v);
    ++clique_sz[clq];
    candidates[clq].intersect_with(g.matrix[v]);
}

int clique_finder::find_cliques()
{
    clear();
    if (g.nodes.size() == 0)
        return 0;
    for (auto u : g.nodes) {
        bool found{false};
        for (int i = 0; i != num_cliques; ++i)
            if (candidates[i].fast_contain(u)) {
                found = true;
                insert(u, i);
            }
        if (!found) {
            new_clique();
            insert(u, num_cliques - 1);
        }
    }
    return *std::max_element(begin(clique_sz), begin(clique_sz) + num_cliques);
}


neighbors_wrapper::neighbors_wrapper(const graph& g)
    : g(g)
		, size(g.nodes.size())
{
    by_degree.resize(g.capacity());
		neighbors.resize(g.capacity());
    degree.resize(g.capacity());
		
		for( auto v : g.nodes ) {
				by_degree[v].reserve(g.capacity());
				neighbors[v].reserve(g.capacity());
				if( !g.matrix[v].empty() ) {
						auto u{0}, unext{g.matrix[v].min()};
						do {
								u = unext;
								neighbors[v].push(u);
								unext = g.matrix[v].next(u);
						} while( u != unext );
				}
		}
}

void neighbors_wrapper::synchronize() { 
		// the nodes betweem g.nodes.size() and size have been removed
		while( size > g.nodes.size() ) {
				auto v = g.nodes[--size];
				for( auto u : neighbors[v] ) {
						neighbors[v].remove( u );
				}
		}	
		// the nodes betweem size and g.nodes.size() have been added
		while( size < g.nodes.size() ) {
				auto v = g.nodes[size++];
				for( auto u : neighbors[v] ) {
						neighbors[v].push( u );
				}
		}
}

void neighbors_wrapper::get_degeneracy_order( std::vector< int >& order ) {
 	
		synchronize();
		for( auto v : g.nodes ) {
				degree[v] = neighbors[v].size();
				by_degree[degree[v]].push( v );
		}
		
		for( auto i = g.nodes.size(); i --> 0; ) {
				for( auto& vertices : by_degree ) {
						if( !vertices.empty() ) {
								auto v = vertices.back();
								vertices.pop_back();	
							
								for( auto ni = 0; ni < degree[v]; ++ni ) {									
										auto u = neighbors[v][ni];
										by_degree[degree[u]].remove( u );
										neighbors[u].move_up(v, --degree[u]);
										by_degree[degree[u]].push( u );
								}
							
								order.push_back(v);
								break;
						}
				}
		}
}




} // namespace gc
