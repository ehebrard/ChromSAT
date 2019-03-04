

#include <algorithm>
#include <cassert>
#include <iomanip>

#include "ca_graph.hpp"
// #include "graph.hpp"
// #include "dsatur.hpp"

namespace gc
{

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
    for (auto v : nodes)
        std::sort(begin(matrix[v]), end(matrix[v]));
}

void ca_graph::remove_neighbor(const int u, const int v)
{
    auto it{end(matrix[u])};
    do
        --it;
    while (*it != v);
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
            std::swap(*r_v, *(--r_v));
        }
    }
}

void ca_graph::contract(const int u, const int v)
{
    parent[v] = u;
    nodes.remove(v);

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
    nodes.add(v);

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

void ca_graph::search(gc::statistics& stats, gc::options& options)
{
    int limit{static_cast<int>(nodes.capacity())};

    check_consistency();

    int depth{0};

    degeneracy_finder df(*this);
    dsatur brelaz;
    brelaz.use_recolor = false;
    df.degeneracy_ordering();

    // assert(df.degeneracy == *std::max_element(begin(df.degrees),
    // end(df.degrees)));

    int lb{1 + (num_edges > 0)};
    stats.notify_lb(lb);

    int ub{df.degeneracy + 1};

    // std::cout << "ub = " << ub << " / " << limit << std::endl;

    std::vector<int> coloring(limit, 0);
    // int pedge[2] = {0, 0};
    bool contraction{true};
    arc pedge{0, 0};
    // decision prev{0,0};
    int nub{ub};

    int cur_lb{lb};

    int period = 1000;

    while (lb < ub) {

        if (options.verbosity > 1 && ++stats.total_iteration % period == 0)
            stats.force_display(std::cout);

        // check_consistency();

        assert(pedge[0] >= 0);

        if (2 * num_edges == (size() * (size() - 1))) {
            // cur_lb = size();
            // if( ub > size() )
            // 	ub = size();

            std::cout << "this happens\n";

            exit(1);
        }

        // std::cout << "dsatur (" << size() << ")";

        if (pedge[0] == pedge[1])
            nub = brelaz.brelaz_color(*this, ub - 1, size(), 12345);
        else if (contraction) {
            coloring[pedge[0]] = 0;
            nub = brelaz.brelaz_color_guided(
                *this, ub - 1, pedge.v, pedge.v + 1, coloring, size(), 12345);
        } else {
            coloring[pedge[0]] = 0;
            coloring[pedge[1]] = 1;
            nub = brelaz.brelaz_color_guided(
                *this, ub - 1, pedge.v, pedge.v + 2, coloring, size(), 12345);
        }

        // std::cout << "nub = " << nub << "/" << size() << std::endl;

        if (nub == size())
            exit(1);

        // assert(brelaz.order.size() == size());
        // for (auto v : brelaz.order) {
        // 		std::cout << " " << v;
        // }
        // std::cout << std::endl;

        // cur_lb = 0;
        int clique_sz{0};
        for (auto v : brelaz.order) {
            // std::cout << v << " " << brelaz.color[v] << std::endl;
            if (brelaz.color[v] < clique_sz) {
                assert(brelaz.color[brelaz.order[brelaz.color[v]]]
                    == brelaz.color[v]);

                // std::cout << "--> branch on " << v << "," <<
                // brelaz.order[brelaz.color[v]] << std::endl;
                pedge = arc{v, brelaz.order[brelaz.color[v]]};
                break;
            }
            ++clique_sz;
        }

        if (clique_sz > cur_lb)
            cur_lb = clique_sz;

        if (depth == 0 and cur_lb > lb)
            lb = cur_lb;

        if (ub > nub) {
            ub = nub;
            stats.notify_ub(ub);
        }

        brelaz.clear();

        // if (size() * (size() - 1) == 2 * num_edges)
        //     cur_lb = size();
        // else
        //     cur_lb = 0;

        if (options.verbosity > 4) {
            if (options.verbosity > 5)
                std::cout << std::endl << *this;

            std::cout //<< std::endl << *this
                << nodes.size() << "/" << num_edges << ": [" << cur_lb << ".."
                << ub << "]" << std::endl;
        }

        //

        if (ub > cur_lb) {
            // pedge = any_non_edge();

            ++depth;

            if (options.verbosity > 3) {
                std::cout << std::setw(3) << std::left << depth << std::right;
                for (int i = 0; i < depth; i++)
                    std::cout << " ";
                std::cout << "/(" << pedge[0] << "," << pedge[1] << ")\n";
            }

            trail.push_back(limit);
            contract(pedge[0], pedge[1]);
            contraction = true;
            // std::cout << *this << std::endl;
        } else if (depth > 0) {

            ++stats.total_conflicts;

            --depth;
            pedge = backtrack(limit);
            addition(pedge[0], pedge[1]);
            contraction = false;

            if (options.verbosity > 3) {
                std::cout << std::setw(3) << std::left << depth << std::right;
                for (int i = 0; i < depth; i++)
                    std::cout << " ";
                std::cout << "+(" << pedge[0] << "," << pedge[1] << ")\n";
            }

            cur_lb = ub - 1;

        } else {
            break;
        }

        // if (stats.total_iteration == 10000)
        //     break;
    }
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

void ca_graph::check_consistency()
{
    auto count{2 * num_edges};
    std::vector<int> f(nodes.capacity(), 0);
    std::vector<int> b(nodes.capacity(), 0);

    for (auto u : nodes) {
        for (auto v : matrix[u]) {
            assert(nodes.contain(v));
            ++f[u];
            ++b[v];
            --count;
        }
    }

    if (count)
        std::cout << count << std::endl;

    assert(count == 0);

    for (auto u : nodes) {
        assert(f[u] == b[u]);
    }
}

std::ostream& operator<<(std::ostream& os, const gc::ca_graph& g)
{
    return g.describe(os, 3);
}
}