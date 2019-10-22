

#include <algorithm>
#include <cassert>
#include <iomanip>

#include "ca_graph.hpp"
// #include "graph.hpp"
// #include "dsatur.hpp"

namespace gc
{

// let's just do it greedily because I'm not sure about the objective function
// anyway
void ca_graph::get_subproblem(std::vector<int>& vertices, const int size_limit)
{
    bi_graph B;
    degeneracy_finder df(*this);
    df.degeneracy_ordering();
    cliquer cq(*this);

    gc::bitset visited(0, nodes.capacity() - 1, gc::bitset::empt);

    auto clique_sz = cq.find_cliques(rbegin(df.order), rend(df.order));

    std::vector<double> unmatched_score(cq.num_cliques, 0);
		

    for (auto i{0}; i < cq.num_cliques; ++i) {
        std::cout << std::setw(2) << cq.cliques[i].size() << " " << std::setw(3)
                  << i << " ";

        for (auto j{0}; j < cq.num_cliques; ++j) 
					if(i!=j)
				{

            B.get_from_cliques(*this, cq.cliques[i], cq.cliques[j]);
            auto mm{B.hopcroftKarp()};

            auto mm_bound
                = (cq.cliques[i].size() + cq.cliques[j].size() - B.I - mm);

            if (mm_bound - std::max(cq.cliques[i].size(), cq.cliques[j].size())
                > 0)
                std::cout //<< " " << std::setw(2)
                    << (mm_bound - std::max(cq.cliques[i].size(),
                                       cq.cliques[j].size()));
            else
                std::cout << " ";

        } else
            std::cout << " ";

        std::cout << std::endl;
    }

    // contains the cliques for which we computed the matching to every other
    // cliques
    intstack explored_cliques;
    explored_cliques.reserve(cq.num_cliques);

    //
    intstack chosen_cliques;
    chosen_cliques.reserve(cq.num_cliques);

    std::vector<std::vector<int>> matching_debt;

    std::vector<int> largest_cliques;
    for (auto i{0}; i < cq.num_cliques; ++i) {
        if (cq.cliques[i].size() == clique_sz) {
            largest_cliques.push_back(i);
        }
    }

    // while (vertices.size() < size_limit) {
    for (auto i : largest_cliques) {

        // std::cout << std::setw(2) << cq.cliques[i].size() << " " <<
        // std::setw(2)
        //           << i;

        matching_debt.resize(matching_debt.size() + 1);
        for (auto j{0}; j < cq.num_cliques; ++j) {
            if (i != j) {
                if (explored_cliques.contain(j)) {

                    matching_debt.back().push_back(
                        matching_debt[explored_cliques.index(j)][i]);
										
										// std::cout << " .";
										
                } else {
                    B.get_from_cliques(*this, cq.cliques[i], cq.cliques[j]);
                    int mm{B.hopcroftKarp()};

                    int mm_bound{static_cast<int>(cq.cliques[i].size()
                                     + cq.cliques[j].size())
                        - B.I - mm};

                    int cq_bound{static_cast<int>(
                        std::max(cq.cliques[i].size(), cq.cliques[j].size()))};

                    assert(mm_bound >= cq_bound);
                    matching_debt.back().push_back(mm_bound - cq_bound);

                    unmatched_score[i] += std::pow(
                        (double)(cq.num_cliques), (mm_bound - cq_bound));
										
										// std::cout << " *";
                }

            } else
                matching_debt.back().push_back(0);

            // std::cout << " " << std::setw(2) << matching_debt.back()[j];
        }

        explored_cliques.add(i);

        // std::cout << " -> " << unmatched_score[i] << std::endl;
    }

    // choose the best next clique
    auto choice{*std::max_element(begin(largest_cliques), end(largest_cliques),
        [&](int x, int y) { return unmatched_score[x] < unmatched_score[y]; })};

    chosen_cliques.add(choice);
    for (auto u : cq.cliques[choice]) {
        visited.add(u);
        vertices.push_back(u);
    }

    // std::cout << " -> CHOOSE " << choice << ": " << visited << std::endl;
    // for (int i{0}; i < cq.num_cliques; ++i)
    //     std::cout << std::setw(5) << i;
    // std::cout << std::endl;

    while (vertices.size() < size_limit) {

        for (auto cl{chosen_cliques.begin_not_in()};
             cl != chosen_cliques.end_not_in(); ++cl) {
            auto i{*cl};
            unmatched_score[i] = 0;

            for (auto j : chosen_cliques)
                unmatched_score[i]
                    += (double)(matching_debt[explored_cliques.index(j)][i]);
        }

        // for (int i{0}; i < cq.num_cliques; ++i)
        //     if (!chosen_cliques.contain(i))
        //         std::cout << std::setw(5) << " ";
        //     else
        //         std::cout << std::setw(5) << (int)(unmatched_score[i]);

        choice = *std::max_element(chosen_cliques.begin_not_in(),
            chosen_cliques.end_not_in(), [&](int x, int y) {
                return unmatched_score[x] < unmatched_score[y];
            });

        // std::cout << std::endl;

        // compute the matching with the chosen clique

        matching_debt.resize(matching_debt.size() + 1);
        for (int i{0}; i < cq.num_cliques; ++i)
            if (i != choice and !chosen_cliques.contain(i)) {
                if (explored_cliques.contain(i)) {

                    matching_debt.back().push_back(
                        matching_debt[explored_cliques.index(i)][choice]);
                } else {

                    B.get_from_cliques(*this, cq.cliques[choice], cq.cliques[i]);
                    int mm{B.hopcroftKarp()};

                    int mm_bound{static_cast<int>(cq.cliques[choice].size()
                                     + cq.cliques[i].size())
                        - B.I - mm};

                    int cq_bound{static_cast<int>(
                        std::max(cq.cliques[choice].size(), cq.cliques[i].size()))};

                    assert(mm_bound >= cq_bound);
                    matching_debt.back().push_back(mm_bound - cq_bound);
                }

            } else
                matching_debt.back().push_back(0);

        chosen_cliques.add(choice);
        if (!explored_cliques.contain(choice))
            explored_cliques.add(choice);
        for (auto u : cq.cliques[choice])
            if (!visited.contain(u)) {
                visited.add(u);
                vertices.push_back(u);
            }

        // std::cout << " -> CHOOSE " << choice << ": " << visited << std::endl;
    }

    // exit(1);
}

int ca_graph::find(const int v) const
{
    auto x{v};
    while (parent[x] != x)
        x = parent[x];
    return x;
}

void ca_graph::add_neighbor(const int u, const int v)
{
    matrix[u].push_back(v);
    auto vptr = rbegin(matrix[u]);
    auto stop = rend(matrix[u]);
    while (++vptr != stop) {
        if (*(vptr - 1) > *vptr)
            break;
        std::swap(*vptr, *(vptr - 1));
    }
}

void ca_graph::addition(const int u, const int v)
{
    assert(parent[u] == u);
    assert(parent[v] == v);
    add_neighbor(u, v);
    add_neighbor(v, u);

    trail.push_back(-u - 1);
    trail.push_back(-v - 1);

    ++num_edges;
}

void ca_graph::add_edge(const int u, const int v)
{
    matrix[u].push_back(v);
    matrix[v].push_back(u);
    ++num_edges;
}

void ca_graph::canonize()
{
    num_edges = 0;
    for (auto v : nodes) {
        std::sort(begin(matrix[v]), end(matrix[v]));
        matrix[v].erase(
            std::unique(begin(matrix[v]), end(matrix[v])), end(matrix[v]));
        num_edges += matrix[v].size();
    }
    num_edges /= 2;
}

void ca_graph::remove_neighbor(const int u, const int v)
{
    auto it{end(matrix[u])};
    do {
        --it;
    } while (*it != v);
    matrix[u].erase(it);
}

void ca_graph::swap_neighbor(const int u, const int v, const int w)
{
    auto r_v{gc::find_in(begin(matrix[u]), end(matrix[u]), v)};

    *r_v = w;
    if (w > v) {
        while (++r_v != end(matrix[u]) and *r_v < *(r_v - 1)) {
            std::swap(*r_v, *(r_v - 1));
        }
    } else {
        while (r_v != begin(matrix[u]) and *r_v < *(r_v - 1)) {
            std::swap(*r_v, *(r_v - 1));
            --r_v;
        }
    }
}

void ca_graph::contract(const int u, const int v)
{
    // check_consistency("before contract");

    parent[v] = u;
    rank[u] += rank[v];
    nodes.remove(v);
    nodeset.remove(v);

    union_buffer.clear();
    inter_buffer.clear();

    // keep track of the current trail's tail
    auto sz{trail.size()};

    // union the neighborhood and put N(v) \ N(u) in the trail
    gc::set_union_delta(begin(matrix[u]), end(matrix[u]), begin(matrix[v]),
        end(matrix[v]), std::back_inserter(union_buffer),
        std::back_inserter(trail), std::back_inserter(inter_buffer));

    num_edges -= inter_buffer.size();

    std::swap(matrix[u], union_buffer);

    // vertices in N(v) \ N(u)
    for (auto w{begin(trail) + sz}; w != end(trail); ++w)
        swap_neighbor(*w, v, u);

    trail.push_back(sz);
    sz = trail.size();

    // vertices in N(u) \cap N(v)
    for (auto w : inter_buffer) {
        remove_neighbor(w, v);
        trail.push_back(w);
    }

    trail.push_back(sz);
    trail.push_back(v);

    // check_consistency("after contract");
}

arc ca_graph::undo()
{
    arc pdecision;
    if (trail.back() < 0) {
        undo_addition();
    } else {
        pdecision = undo_contraction();
    }

    return pdecision;
}

arc ca_graph::backtrack(const int plevel)
{
    arc pdecision;
    while (trail.back() != plevel) {
        pdecision = undo();
    }
    trail.pop_back();

    return pdecision;
}

arc ca_graph::undo_addition()
{
    int v{-trail.back() - 1};
    trail.pop_back();
    int u{-trail.back() - 1};
    trail.pop_back();

    remove_neighbor(u, v);
    remove_neighbor(v, u);

    --num_edges;

    return arc{u, v};
}

arc ca_graph::undo_contraction()
{
    int v{trail.back()};
    trail.pop_back();
    int u{parent[v]};
    parent[v] = v;
    rank[u] -= rank[v];
    nodes.add(v);
    nodeset.add(v);

    int sz = trail.back();
    trail.pop_back();
    for (auto w{begin(trail) + sz}; w != end(trail); ++w)
        add_neighbor(*w, v);

    num_edges += (trail.size() - sz);

    trail.resize(sz);
    sz = trail.back();
    trail.pop_back();

    for (auto w{begin(trail) + sz}; w != end(trail); ++w)
        swap_neighbor(*w, u, v);

    union_buffer.clear();
    std::set_difference(begin(matrix[u]), end(matrix[u]), begin(trail) + sz,
        end(trail), std::back_inserter(union_buffer));
    std::swap(matrix[u], union_buffer);
    trail.resize(sz);

    return arc{u, v};
}

arc ca_graph::any_non_edge()
{
    for (auto u : nodes) {
        if (matrix[u].size() < nodes.size() - 1) {
            int w{0};
            while (!nodes.contain(w) or u == w)
                ++w;
            for (auto v : matrix[u]) {
                if (v > w) {
                    return arc{u, w};
                } else {
                    w = v + 1;
                    while (!nodes.contain(w) or u == w)
                        ++w;
                }
            }
        }
    }
    return arc{0, 0};
}


std::ostream& ca_graph::describe(std::ostream& os, const int verbosity) const
{
    // void gc::ca_graph::describe(std::ostream& os, const int verbosity) const
    // {
    switch (verbosity) {
    case 0:
        os << "n=" << nodes.size() << " m=" << num_edges << std::endl;
        break;

    case 1:
        os << "V=" << nodes << " m=" << num_edges;
        break;

    case 2:
        for (auto v : nodes) {
            os << v << ":";
            for (auto u : matrix[v]) {
                os << " " << u;
            }
            os << std::endl;
        }
        break;

    case 3:
        for (auto v : nodes) {
            os << v << ":";
            for (auto u : matrix[v]) {
                os << " " << u;
            }
            os << std::endl;
        }
        for (auto t : trail)
            std::cout << " " << t;
        std::cout << std::endl;
    }

    return os;
}

void ca_graph::check_consistency(const char* msg)
{
    auto count{2 * num_edges};
    std::vector<int> f(nodes.capacity(), 0);
    std::vector<int> b(nodes.capacity(), 0);

    if (nodeset.size() != nodes.size()) {
        std::cout << msg << ": |nodes|=" << nodes.size()
                  << " and |nodeset|=" << nodeset.size() << std::endl;
        exit(1);
    }
    // assert(nodeset.size() == nodes.size());

    int total_rank{0};
    for (auto u : nodes) {
        if (!nodeset.contain(u)) {
            std::cout << msg << ": " << u << " in nodes but not in nodeset\n";
            exit(1);
        }
        // assert(nodeset.contain(u));
        total_rank += rank[u];

        auto prev{-1};
        for (auto v : matrix[u]) {

            if (!nodeset.contain(u)) {
                std::cout << msg << ": " << v << " in N(" << u
                          << ") but not in nodes\n"
                          << std::endl;
                exit(1);
            }

            if (v <= prev) {
                std::cout << msg << ": nodes not sorted N(" << u << ") =";
                for (auto x : matrix[u]) {
                    std::cout << " " << x;
                    if (x == v) {
                        std::cout << "...\n";
                        break;
                    }
                }
                exit(1);
            }

            // assert(nodes.contain(v));
            ++f[u];
            ++b[v];
            --count;
        }
    }

    assert(total_rank == nodes.capacity());

    if (count)
        std::cout << count << std::endl;

    assert(count == 0);

    for (auto u : nodes) {

        if (f[u] != b[u]) {
            std::cout << msg << ": " << u << " has " << f[u]
                      << " neighbors but appears in " << b[u]
                      << " neighborhoods\n"
                      << std::endl;
            exit(1);
        }
        // assert(f[u] == b[u]);
    }
}

// Returns size of maximum matching
int bi_graph::hopcroftKarp()
{
    // Initialize result
    int result = 0;

    // Keep updating the result while there is an
    // augmenting path.
    while (bfs()) {

        // std::cout << "\nmatching (" << result << "):\n";
        //
        // for (int u = 0; u < N; ++u) {
        //     if (matching[u] != T)
        //         std::cout << u << " -- " << matching[u] << std::endl;
        // }
        //
        // std::cout << "BFS\n";

        // Find a free vertex
        for (int u = N; u < T; ++u)
            // If current vertex is free and there is
            // an augmenting path from current vertex
            if (matching[u] == T && dfs(u))
                result++;

        // for (int u = 0; u < T; ++u) {
        //     std::cout << std::setw(2) << u << " " << std::setw(3) <<
        //     original[u]
        //               << " ";
        //     if (dist[u] != INFTY)
        //         std::cout << dist[u] << std::endl;
        //     else
        //         std::cout << "inf" << std::endl;
        // }
        // if (dist[T] != INFTY)
        //     std::cout << std::setw(2) << T << "     " << dist[T] <<
        //     std::endl;
        // else
        //     std::cout << std::setw(2) << T << "     inf\n";
        // std::cout << std::endl;
    }
    return result;
}

// Returns true if there is an augmenting path, else returns
// false
bool bi_graph::bfs()
{

    // First layer of vertices (set distance as 0)
    for (int u = N; u < T; ++u) {
        // If this is a free vertex, add it to queue
        if (matching[u] == T) {
            // u is not matched
            dist[u] = 0;
            Q.push_back(u);
            // std::cout << " " << u;
        } else {
            dist[u] = INFTY;
            // std::cout << " -" << u;
        }
    }

    // std::cout << std::endl;

    // Initialize distance to T as infinite
    dist[T] = INFTY;

    // Q is going to contain vertices of left side only.
    while (q < Q.size()) {
        // Dequeue a vertex
        int u = Q[q++];

        // std::cout << u << ":";

        // If this node is not T and can provide a shorter path to T
        if (dist[u] < dist[T]) {
            // Get all adjacent vertices of the dequeued vertex u
            for (auto v : matrix[u]) {
                // If pair of v is not considered so far
                // (v, pairV[V]) is not yet explored edge.
                if (dist[matching[v]] == INFTY) {
                    // Consider the pair and add it to queue
                    dist[matching[v]] = dist[u] + 1;
                    Q.push_back(matching[v]);

                    // std::cout << " " << matching[v];
                }
            }
        }

        // std::cout << std::endl;
    }

    // If we could come back to T using alternating path of distinct
    // vertices then there is an augmenting path
    return (dist[T] != INFTY);
}

// Returns true if there is an augmenting path beginning with free vertex u
bool bi_graph::dfs(const int u)
{
    if (u != T) {
        for (auto v : matrix[u]) {
            if (dist[matching[v]] == dist[u] + 1) {
                if (dfs(matching[v])) {
                    matching[u] = v;
                    matching[v] = u;
                    return true;
                }
            }
        }

        dist[u] = INFTY;
        return false;
    }

    return true;
}

std::ostream& bi_graph::describe(std::ostream& os) const
{

    for (int i = 0; i < T; ++i) {
        if (i == N)
            std::cout << std::endl;

        os << original[i] << ":";
        for (auto j : matrix[i]) {
            os << " " << original[j];
        }
        os << std::endl;
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const gc::ca_graph& g)
{
    return g.describe(os, 3);
}

std::ostream& operator<<(std::ostream& os, const gc::bi_graph& g)
{
    return g.describe(os);
}
}