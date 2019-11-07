#include <iomanip>
#include <iomanip>
#include <iostream>

#include "graph.hpp"
#include "heap.hpp"
#include "intstack.hpp"
#include "options.hpp"
#include "partition.hpp"
#include "statistics.hpp"

#include <boost/dynamic_bitset.hpp>

#include <minicsp/core/utils.hpp>

#ifndef __COLORING_HEURISTIC_HPP
#define __COLORING_HEURISTIC_HPP

// #define _DEBUG_DSATUR
// #define CHECK

namespace gc
{

template <typename T> class coldomain
{

public:
    coldomain();
    coldomain(const int ub);

    // empty domain of size ub
    void initialise(const int ub);
    void resize(const size_t s);

    inline size_t size() const;
    inline T count() const;
    int get_first_allowed(const int prev = -1) const;
    bool contain(const int elt) const;
    T weight(const int elt) const;

    // add w neighborhing color
    bool add(const int c, const T w = 1);
    bool remove(const int c, const T w = 1);
    void clear();

    std::ostream& display(std::ostream& os) const;

private:
    std::vector<T> b;
    T count_;
};

//..............{D+1}{-}(d+1)(-)x....(d)......(d-1)(d-2)(d-3).....(d-4)...(-)
//[-----------colored----------][---------------uncolored-------------------]

class coloring_heuristic
{

public:
    const std::vector<int>& get_coloring() const;

    // initialise the structures w. r. t. graph g and upper bound ub
    template <class graph_struct, class compare>
    void initialise(graph_struct& g, const int ub, compare comp);

    template <class graph_struct>
    void assign_color(graph_struct& g, const int x, const int c);

    template <class graph_struct>
    void unassign_color(graph_struct& g, const int x);

    // satur[y] was d, and is now d+1
    void move_up(const int y, const int d);

    // satur[y] was d+1, and is now d
    void move_down(const int y, const int d);

    template <class graph_struct, typename tiebreaker>
    int brelaz_greedy(graph_struct& g, const int ub,
        std::vector<int>::iterator start, const int limit,
        tiebreaker criterion);

    template <class graph_struct>
    int dsatur(
        graph_struct& g, const int ub, const int limit = 1, const int seed = 1);

    template <class graph_struct>
    int sdatur(
        graph_struct& g, const int ub, const int limit = 1, const int seed = 1);

    // update the saturation degrees; reorder by increasing SD; fix the
    // last_vertex pointers
    template <class graph_struct> void close(graph_struct& g);

    template <class graph_struct>
    int degeneracy(graph_struct& g, const bool use_recolor = false);

    // try to find a color for y within x's neighborhood
    template <class graph_struct>
    bool recolor(graph_struct& g, const int y, const int x);

    // try to find a color for x within y's neighborhood
    template <class graph_struct> bool reduce(graph_struct& g, const int x);

    void clear();

    template <class graph_struct> void print(graph_struct& g);

private:
    // // store the best full coloring
    // std::vector<int> best;

    // number of distinct colors in 'best_coloring'
    size_t best_numcolors;

    // current (partial)  coloring
    std::vector<int> color;
    // number of distinct colors in 'color'
    size_t numcolors;

    // store info about the colors of the niehgbors in 'color'
    std::vector<coldomain<int>> neighbor_colors;

    // degree of the vertices in the graph
    std::vector<int> degree;

    // ordered list of the vertices of the graph
    std::vector<int> order;

    // rank of every vertex in 'order'
    std::vector<std::vector<int>::iterator> rank;

    // 'last_vertex[d]' points just after the last unassigned vertex of
    // saturation degree d
    std::vector<std::vector<int>::iterator> last_vertex;

    std::vector<int>::iterator beg_update;
    std::vector<int>::iterator end_update;

    std::mt19937 random_generator;

    template <class graph_struct> void check_consistency(graph_struct& g);
};

class independent_set_heuristic
{

public:
    std::vector<int> order;
    boost::dynamic_bitset<> candidate;

    template <class graph_struct, class RandomIt, class Compare>
    void initialise(graph_struct& g, RandomIt beg, RandomIt end, Compare c);

    template <class graph_struct, class container>
    void find_is(graph_struct& g, container is);

    void seed(const int s) { random_generator.seed(s); }

private:
    std::mt19937 random_generator;
};

template <typename T> class multi_coloring_heuristic
{

public:
    template <class graph_struct, class map_struct>
    void greedy_color(graph_struct& g, map_struct& w);

    template <class graph_struct, class map_struct>
    void greedy_uncolor(graph_struct& g, map_struct& w);

    T numcolors{0};
    T degeneracy{0};

    // the current set of ISs
    std::vector<std::vector<int>> independent_set;
    // the converse structure
    std::vector<std::vector<int>> colors;

    // for each IS, its contribution to the dual cost
    std::vector<T> ncolors;

    //
    std::vector<coldomain<T>> neighbor_colors;

    void seed(const int s) { ish.seed(s); }

private:
    // the residual weight after some vertices have been colored
    std::vector<T> residual_weight;

    independent_set_heuristic ish;
    intstack nodes;

    // // std::vector<int> order;
    // heap::Heap<int> order;
};

template <typename T>
template <class graph_struct, class map_struct>
void multi_coloring_heuristic<T>::greedy_color(graph_struct& g, map_struct& w)
{

    nodes.reserve(g.size());
    nodes.fill();

    residual_weight.resize(g.size());
    for (auto v : g.nodes) {
        residual_weight[v] = w[v];
        if (residual_weight[v] == 0)
            nodes.remove(v);
    }

    while (!nodes.empty()) {
        independent_set.resize(independent_set.size() + 1);

        ish.initialise(
            g, nodes.begin(), nodes.end(), [=](const int x, const int y) {
                // return g.matrix[x].size() * residual_weight[y] <
                // g.matrix[y].size() * residual_weight[x];
                return residual_weight[y] < residual_weight[x];
            });
        ish.find_is(g, std::back_inserter(independent_set.back()));

        T minw{residual_weight[independent_set.back().back()]};
        for (auto v : independent_set.back())
            if (residual_weight[v] < minw)
                minw = residual_weight[v];

        for (auto v : independent_set.back()) {
            // std::cout << " " << v;
            residual_weight[v] -= minw;
            if (residual_weight[v] == 0)
                nodes.remove(v);
        }

        // std::cout << " -> " << minw << std::endl;

        ncolors.push_back(minw);
        numcolors += minw;
    }
}

template <typename T>
template <class graph_struct, class map_struct>
void multi_coloring_heuristic<T>::greedy_uncolor(graph_struct& g, map_struct& w)
{
    colors.clear();
    colors.resize(g.size());
    for (auto i{0}; i < independent_set.size(); ++i) {
        // std::cout << i << " (" << ncolors[i] << "):";
        for (auto v : independent_set[i]) {
            colors[v].push_back(i);
            // std::cout << " " << v;
        }
        // std::cout << std::endl;
    }

    neighbor_colors.resize(g.size());
    for (auto& n : neighbor_colors) {
        n.clear();
        n.initialise(colors.size());
    }

    std::vector<int> order;
    // std::vector<int> reverse;
    for (auto v : g.nodes) {
        order.push_back(v);
        // reverse.push_back(v);
        for (auto u : g.matrix[v])
            for (auto c : colors[v])
                neighbor_colors[u].add(c, ncolors[c]);
    }

#ifdef CHECK
    for (auto v : g.nodes) {
        T total_v{0};
        for (auto c : colors[v])
            total_v += ncolors[c];
        assert(w[v] == total_v);

        T total_n{0};
        std::vector<int> nc;
        for (auto u : g.matrix[v])
            for (auto c : colors[u])
                nc.push_back(c);

        std::sort(nc.begin(), nc.end());
        nc.erase(std::unique(nc.begin(), nc.end()), nc.end());

        for (auto c : nc)
            total_n += ncolors[c];

        assert(neighbor_colors[v].count() == total_n);
        // {}

        // std::cout << v << ": " << neighbor_colors[v] << std::endl;
    }
#endif

    auto comp{[=](const int x, const int y) {
        return (w[x] + neighbor_colors[x].count())
            < (w[y] + neighbor_colors[y].count());
    }};
    heap::heapify(order.begin(), order.end(), comp);

    std::vector<int> rank(g.size());
    auto r{0};
    for (auto v : order) {
        // std::cout << neighbor_colors[v] << std::endl;
        rank[v] = r++;
    }

    // std::cout << g << std::endl;

    // T cdegeneracy{0};
    while (not order.empty()) {

        order.pop_back();
        std::swap(*begin(order), *end(order));
        auto i{heap::percolate_down(begin(order), end(order), 0, comp)};

        // std::cout << *end(order) << ": 0 -> " << i << std::endl;

        // for()
        while (true) {

            // std::cout << "rank[" << order[i] << "] = " << i << std::endl;

            rank[order[i]] = i;
            if (i == 0)
                break;
            i = --i / 2;
        }

        auto v{*end(order)};

        auto d{w[v] + neighbor_colors[v].count()};

        degeneracy = std::max(degeneracy, d);

        // std::cout << "uncolor " << v << ": (" << (w[v] +
        // neighbor_colors[v].count()) << ")\n " ;//<< neighbor_colors[v] <<
        // std::endl;

        for (auto u : g.matrix[v]) {
            if (rank[u] < order.size()) {
                auto update{false};

                // auto old{neighbor_colors[u].count()};

                for (auto c : colors[v])
                    update |= neighbor_colors[u].remove(c, ncolors[c]);

                if (update) {

                    // std::cout << "update " << u << " was " << old << " now "
                    //           << neighbor_colors[u].count() << " (+" << w[u]
                    //           << ")" << std::endl;

                    auto j{heap::percolate_up(
                        begin(order), end(order), rank[u], comp)};
                    for (auto k{rank[u]};; k = --k / 2) {

                        // std::cout << "rank[" << order[k] << "] = " << k <<
                        // "("
                        //           << neighbor_colors[order[k]].count() << " +
                        //           "
                        //           << w[order[k]] << ")" << std::endl;

                        rank[order[k]] = k;
                        if (k == j)
                            break;
                    }
                }
            }
        }

// heap::rank(begin(order), 0, i, rank);

#ifdef CHECK
        auto r{0};
        for (auto v : order) {
            // std::cout << v << std::endl;
            // std::cout << neighbor_colors[v] << std::endl;
            // assert(rank[v] == r++);
            if (rank[v] != r++) {
                std::cout << "rank[" << v << "] != " << (r - 1) << std::endl;
                exit(1);
            }

            auto u{*begin(order)};

            if ((w[v] + neighbor_colors[v].count())
                < (w[u] + neighbor_colors[u].count())) {

                std::cout << u << " (" << w[u] << " + "
                          << neighbor_colors[u].count()
                          << ") is first in the heap, however " << v
                          << " is s.t. (" << u << " (" << w[v] << " + "
                          << neighbor_colors[v].count() << ")\n";

                exit(1);
            }
        }
#endif

        // exit(1);
    }

    // order.resize(g.size());

    // for(auto v : reverse)
    // 	std::cout << neighbor_colors[v] << std::endl;

    // return cdegeneracy;
}

/************* IMPlEMENTATION **************/

template <typename T> coldomain<T>::coldomain() {}

template <typename T>
coldomain<T>::coldomain(const int ub)
    : b(ub, 0)
{
    count_ = 0;
}

template <typename T> void coldomain<T>::resize(const size_t s)
{
    b.resize(s, 0);
}

template <typename T> void coldomain<T>::initialise(const int ub)
{
    b.clear();
    count_ = 0;
    b.resize(ub, 0);
}

template <typename T> T coldomain<T>::count() const { return count_; }

template <typename T> size_t coldomain<T>::size() const { return b.size(); }

template <typename T> int coldomain<T>::get_first_allowed(const int prev) const
{
    for (auto i{begin(b) + prev + 1}; i != end(b); ++i)
        if (*i == 0)
            return i - begin(b);
    return b.size();
}

template <typename T> bool coldomain<T>::contain(const int elt) const
{
    return b[elt] > 0;
}

template <typename T> T coldomain<T>::weight(const int elt) const
{
    return b[elt];
}

template <typename T> bool coldomain<T>::add(const int x, const T w)
{
    b[x] += w;
    if (b[x] == w) {
        count_ += w;
        return true;
    }
    return false;
}

template <typename T> bool coldomain<T>::remove(const int x, const T w)
{
    b[x] -= w;
    if (b[x] == 0) {
        count_ -= w;
        return true;
    }
    return false;
}

template <typename T> void coldomain<T>::clear()
{
    for (auto i{begin(b)}; i != end(b); ++i)
        *i = 0;
    count_ = 0;
}

template <typename T>
std::ostream& coldomain<T>::display(std::ostream& os) const
{
    os << std::setw(5) << count() << " [";
    for (auto i{begin(b)}; i != end(b); ++i)
        if (*i != 0)
            os << " " << (i - begin(b)) << "(" << *i << ")";
    os << " ]";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const coldomain<T>& x)
{
    x.display(os);
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const coldomain<T>* x)
{
    x->display(os);
    return os;
}

const std::vector<int>& coloring_heuristic::get_coloring() const
{
    return color;
}

// initialise the structures w. r. t. graph g and upper bound ub
template <class graph_struct, class compare>
void coloring_heuristic::initialise(graph_struct& g, const int ub,
    compare comp) //=[&](const int, const int) { return 1; }
{
    numcolors = 0;
    best_numcolors = 0;
    rank.resize(g.capacity());
    color.resize(g.capacity(), -1);
    degree.resize(g.capacity());

    for (auto v : g.nodes) {
        order.push_back(v);
        degree[v] = g.matrix[v].size();
    }
    std::shuffle(begin(order), end(order), random_generator);
    std::sort(begin(order), end(order), comp);
    for (auto vptr{begin(order)}; vptr != end(order); ++vptr) {
        rank[*vptr] = vptr;
    }

    // reuse the available memory
    for (auto it{begin(neighbor_colors)}; it != end(neighbor_colors); ++it) {
        it->initialise(ub);
    }
    // allocate as much as needed
    neighbor_colors.resize(g.capacity(), coldomain<int>(ub));

    // the saturation degree is null for every vertex
    last_vertex.resize(ub + 2, begin(order));
    // hence last_vertex[0] oints at the end of order
    *begin(last_vertex) = end(order);

    // status = 1;
}

template <class graph_struct>
void coloring_heuristic::assign_color(graph_struct& g, const int x, const int c)
{
    color[x] = c;

    // update the saturation degree of x's neighbors
    for (auto y : g.matrix[x])
        if (rank[y] >= beg_update and rank[y] < end_update)
            if (neighbor_colors[y].add(c))
                move_up(y, neighbor_colors[y].count());
}

template <class graph_struct>
void coloring_heuristic::unassign_color(graph_struct& g, const int x)
{
    auto c{color[x]};
    color[x] = -1;

    // update the saturation degree of x's neighbors
    for (auto y : g.matrix[x])
        if (rank[y] >= beg_update and rank[y] < end_update)
            if (neighbor_colors[y].remove(c))
                move_down(y, neighbor_colors[y].count() + 1);
}

// satur[y] was d, and is now d+1
void coloring_heuristic::move_up(const int y, const int d)
{
    // swap y with *last_vertex[d]
    auto l{*last_vertex[d]};

    rank[l] = rank[y];
    rank[y] = last_vertex[d];

    *rank[y] = y;
    *rank[l] = l;

    ++last_vertex[d];
}

// satur[y] was d+1, and is now d
void coloring_heuristic::move_down(const int y, const int d)
{

    assert(d >= 0);
    assert(d < last_vertex.size());
    assert(last_vertex[d] != begin(order));

    // swap y with *last_vertex[d]-1
    auto l{*(--last_vertex[d])};

    rank[l] = rank[y];
    rank[y] = last_vertex[d];

    *rank[y] = y;
    *rank[l] = l;
}

template <class graph_struct, typename tiebreaker>
int coloring_heuristic::brelaz_greedy(graph_struct& g, const int ub,
    std::vector<int>::iterator start, const int limit, tiebreaker criterion)
{
    beg_update = start;
    end_update = end(order);

    int potential_colors = begin(neighbor_colors)->size();

    int c, d;

    std::vector<int>::iterator candidate{start};

    while (candidate != end(order)) {

#ifdef _DEBUG_DSATUR
        check_consistency(g);
#endif

        // get the highest saturation degree
        d = neighbor_colors[*candidate].count();

        if (limit > 1) {
            auto best{std::min_element(candidate,
                std::min(last_vertex[d], candidate + limit), criterion)};

            std::swap(rank[*best], rank[*candidate]); // not sure this is useful
            std::swap(*best, *candidate);
        }

        // use the first possible color for x
        c = neighbor_colors[*candidate].get_first_allowed();

        if (c == numcolors) {
            assert(potential_colors > c);
            ++numcolors;
        }

        if (numcolors > ub) {
            color[*candidate] = c;
            return g.size();
        }

        // move all the pointers >= d
        while (++d < last_vertex.size())
            ++last_vertex[d];

        ++beg_update;
        assign_color(g, *candidate, c);

        // update degrees
        if (limit > 1)
            for (auto y : g.matrix[*candidate])
                --degree[y];

        ++candidate;
    }

    return numcolors;
}

template <class graph_struct>
int coloring_heuristic::dsatur(
    graph_struct& g, const int ub, const int limit, const int seed)
{
    if (g.nodes.empty())
        return 0;

    random_generator.seed(seed);

    initialise(g, ub,
        [&](const int x_, const int y_) { return (degree[x_] > degree[y_]); });

    return brelaz_greedy(g, ub, begin(order), limit,
        [&](int x, int y) { return degree[x] > degree[y]; });
}

template <class graph_struct>
int coloring_heuristic::sdatur(
    graph_struct& g, const int ub, const int limit, const int seed)
{
    if (g.nodes.empty())
        return 0;

    random_generator.seed(seed);

    initialise(g, ub,
        [&](const int x_, const int y_) { return (degree[x_] < degree[y_]); });

    return brelaz_greedy(g, ub, begin(order), limit,
        [&](int x, int y) { return degree[x] > degree[y]; });
}

// update the saturation degrees; reorder by increasing SD; fix the
// last_vertex pointers
template <class graph_struct> void coloring_heuristic::close(graph_struct& g)
{

    // update the color neighborhood ()
    for (auto v : order) {
        for (auto u : g.matrix[v])
            if (rank[u] < rank[v])
                neighbor_colors[u].add(color[v]);
        degree[v] = g.matrix[v].size();
    }

    std::sort(begin(order), end(order), [&](const int x_, const int y_) {
        return (neighbor_colors[x_].count() > neighbor_colors[y_].count()
            or (neighbor_colors[x_].count() == neighbor_colors[y_].count()
                   and degree[x_] > degree[y_]));
    });

    for (auto vptr{begin(order)}; vptr != end(order); ++vptr) {
        rank[*vptr] = vptr;
    }

    last_vertex.resize(numcolors + 1);
    auto d{0};
    for (auto it{end(order)}; it != begin(order); --it) {
        auto v{*(it - 1)};
        while (d < numcolors and neighbor_colors[v].count() >= d) {
            last_vertex[d++] = it;
        }
    }
    last_vertex.back() = begin(order);
}

template <class graph_struct>
int coloring_heuristic::degeneracy(graph_struct& g, const bool use_recolor)
{
    beg_update = begin(order);
    end_update = end(order);

    int dgn{0};
    for (auto vptr{rbegin(order)}; vptr != rend(order); ++vptr) {
        int d = neighbor_colors[*vptr].count();

        if (use_recolor) {
            for (auto it{rank[*vptr]}; it >= last_vertex[d + 1]; --it)
                if (reduce(g, *it))
                    break;
            d = neighbor_colors[*vptr].count();
        }

        dgn = std::max(dgn, d);

        // move all the pointers >= d
        while (d >= 0)
            --last_vertex[d--];

        --end_update;
        unassign_color(g, *vptr);
    }

    return dgn;
}

// try to find a color for y within x's neighborhood
template <class graph_struct>
bool coloring_heuristic::recolor(graph_struct& g, const int y, const int x)
{
    for (auto c{0}; c < numcolors; ++c) {
        if (color[y] != c and !neighbor_colors[y].contain(c)
            and neighbor_colors[x].contain(c)) {
            unassign_color(g, y);
            assign_color(g, y, c);
            return true;
        }
    }
    return false;
}

// try to find a color for x within y's neighborhood
template <class graph_struct>
bool coloring_heuristic::reduce(graph_struct& g, const int x)
{
    int reduction{0};
    for (auto y : g.matrix[x])
        if (color[y] >= 0 and neighbor_colors[x].weight(color[y]) == 1)
            // neighbor_colors[x].b[color[y]] == 1)
            reduction += recolor(g, y, x);
    return reduction > 0;
}

void coloring_heuristic::clear()
{
    last_vertex.clear();
    color.clear();
    if (neighbor_colors.size() > 0)
        for (auto v : order)
            neighbor_colors[v].clear();
    order.clear();
}

template <class graph_struct> void coloring_heuristic::print(graph_struct& g)
{

    // std::cout << std::endl;
    //
    // // dd = neighbor_colors[*begin(order)].count()+1;
    // //
    // //  for (auto v{begin(order)}; v != end(order); ++v) {
    // //
    // //      if (v == last_vertex[dd]) {
    // //          std::cout << "-------(" << dd << ")" << std::endl;
    // //  								--dd;
    // //      }
    // //      std::cout << std::setw(3) << *v << " " << std::setw(3)
    // //                << neighbor_colors[*v].count() << " " <<
    // std::setw(3)
    // //                << degree[*v] << std::endl;
    // //  }

    int d = last_vertex.size() - 1;
    for (auto r{begin(order)}; r != end(order); ++r) {

        bool lim{false};
        while (d and last_vertex[d] == r) {
            lim = true;
            std::cout << "start[" << d - 1 << "] ";
            --d;
        }
        if (lim) {
            std::cout << std::endl;
        }
        auto v{*r};

        std::cout << std::setw(3) << v << ": ";

        std::cout << "(" << neighbor_colors[v].count();

        std::cout << "|" << degree[v] << ") " << neighbor_colors[v];

        if (color[v] >= 0) {
            std::cout << " ** " << color[v] << " **";
        }

        std::cout << std::endl;
    }
}

template <class graph_struct>
void coloring_heuristic::check_consistency(graph_struct& g)
{

    // assert(!full);

    for (size_t d{last_vertex.size() - 1}; d > 0; --d) {
        assert(last_vertex[d] <= last_vertex[d - 1]);
        for (auto r{last_vertex[d]}; r != last_vertex[d - 1]; ++r) {

            if (color[*r] < 0 and neighbor_colors[*r].count() != (d - 1)) {

                std::cout << *r << " has satur degree "
                          << neighbor_colors[*r].count() << " but is in bucket "
                          << (d - 1) << std::endl;
            }

            assert(color[*r] >= 0 or neighbor_colors[*r].count() == (d - 1));
        }
    }

    std::vector<int> colv(numcolors);
    for (auto r{begin(order)}; r != end(order); ++r) {
        auto v{*r};
        auto d{neighbor_colors[v].count()};

        if (color[v] < -1) {
            for (auto c{0}; c < numcolors; ++c) {
                colv[c] = neighbor_colors[v].weight(c); // b[c];
            }
            for (auto u : g.matrix[v]) {
                if (color[u] >= 0 and rank[u] < rank[v])
                    --colv[color[u]];
            }

            for (auto c{0}; c < numcolors; ++c) {

                if (colv[c] != 0) {

                    std::cout << "problem in coldomain of " << v << " (" << c
                              << ") ";
                    neighbor_colors[v].display(std::cout);
                    std::cout << std::endl;

                    for (auto b{0}; b < numcolors; ++b) {
                        std::cout << " "
                                  << neighbor_colors[v].weight(c); //.b[b];
                    }
                    // std::cout << std::endl << "--------\n";

                    for (auto u : g.matrix[v]) {
                        if (color[u] >= 0 and rank[u] < rank[v])
                            std::cout << " " << u << " " << color[u]
                                      << std::endl;
                    }

                    exit(1);
                }

                assert(colv[c] == 0);
            }
        }

        if (color[v] >= 0) {
            for (auto u : g.matrix[v]) {

                if (color[u] == color[v]) {
                    std::cout << "N(" << v << ") ="; //<< g.matrix[v]
                    for (auto w : g.matrix[v])
                        std::cout << " " << w;
                    std::cout << std::endl
                              << "ERROR: " << u << ":=" << color[u] << " and "
                              << v << ":=" << color[v] << std::endl;
                }

                assert(color[u] != color[v]);

                if (color[u] < 0) {
                    if (!neighbor_colors[u].contain(color[v]))
                        std::cout << "ERROR: NC(" << u << ") = ";
                    neighbor_colors[u].display(std::cout);
                    //<< neighbor_colors[u]
                    std::cout << " - c[" << v << "] = " << color[v]
                              << std::endl;

                    assert(neighbor_colors[u].contain(color[v]));
                }
            }
        } else {
            assert(last_vertex[d] > r);
            assert(d + 1 == last_vertex.size() or last_vertex[d + 1] <= r);
        }
    }
}

template <class graph_struct, class RandomIt, class Compare>
void independent_set_heuristic::initialise(
    graph_struct& g, RandomIt b, RandomIt e, Compare c)
{
    order.clear();
    candidate.reset();

    order.reserve(e - b);
    candidate.resize(g.size(), false);
    for (auto vp{b}; vp != e; ++vp) {
        order.push_back(*vp);
        candidate.set(*vp);
    }

    std::shuffle(begin(order), end(order), random_generator);
    std::sort(begin(order), end(order), c);
}

template <class graph_struct, class container>
void independent_set_heuristic::find_is(graph_struct& g, container is)
{
    // initialise(g, beg, end, [=](const int x) {return g.matrix[x].size();});

    for (auto v : order) {
        if (candidate[v]) {
            is = v;
            ++is;
            candidate.reset(v);
            for (auto u : g.matrix[v])
                candidate.reset(u);
        }
    }
}

// class DualProblem
// {
// public:
//     // the graph
//     Graph* G;
//     // number of virtual/real cliques in the current clique cover(s)
//     weight_type dual_cost{0};
//     // cumulative size of each "layer" of clique covers
//     std::vector<int> layer_csize;
//     // // weight of each "layer" of clique covers
//     // std::vector<int> layer_weight;
//     // whether the clique cover is sorted
//     std::vector<int> sorted;
//     // the current clique cover
//     std::vector<BitSet> partition;
//     // also, maintains its size because each call to size() is costly
//     std::vector<int> partition_sz;
//     // the sets of vertices that could possibly be added to the cliques
//     std::vector<BitSet> candidate;
//     // also, maintains its size because each call to size() is costly
//     std::vector<int> candidate_sz;
//     // for each partition, its contribution to the dual cost
//     std::vector<weight_type> partition_weight;
//     // for each node, the partitions it belongs to and how much does it
//     // contribute to it
//     std::vector<std::vector<int>> partition_of;
//     std::vector<std::vector<weight_type>> weight_contribution;
//     // helper to compute the residual weights / partition_of efficiently
//     std::vector<std::vector<int>> partition_list;
//     std::vector<int> singleton_partitions;
//     // std::vector<int> partition_of;
//     // bool is_sorted;
//
//     DualProblem(Graph* g);
//     virtual ~DualProblem();
//     // void initialise(Graph* g);
//
//     void resize(const int n);
//
//     // insert x in the first available clique and return its index
//     int insert(const int x);
//     // insert x in the first available clique, return its index and the
//     removed
//     // candidate
//     int insert(const int x, BitSet& delta);
//     // compute all candidate size
//     void close(std::vector<int>&, std::vector<weight_type>&);
//     // compute all candidate sizes, but consume weight less eagerly
//     void close_weak(std::vector<int>&, std::vector<weight_type>&);
//     // add a new clique with unique element x
//     int push(const int x, const int nx);
//     // add a new clique G
//     void push(const Graph& G);
//
//     void repair(BitSet& util_set);
//     // try to remove the ith clique by moving its vertices to subsequent
//     // cliques -- explore the cliques in the order given by "sorted"
//     void remove_partition(std::vector<int>::iterator cl, BitSet& util_set);
//     // update the clique cover in response to "removed_nodes", adds the
//     cliques
//     // that have changes in "changed_partitions"
//     void update(BitSet& removed_nodes, BitSet& util_set);
//     void clear();
//     void to_size(const int k);
//     // void clear(const int actual_size);
//     // set the order of the clique indices in "sorted" by non-increasing
//     // candidate-size (ties broken by non-decreasing partition-size)
//     void sort();
//     // move nodes to lower cliques greedily
//     void percolate(BitSet& delta);
//     void check_conflict_set(std::vector<int>& cset, std::vector<BitSet>&
//     part,
//         bool cons, bool verbose);
//
//     int get_partition(const int v) const;
//     int num_partition() const;
//     void new_layer();
//     // void new_clique();
//
//     void check(const char* msg, const bool covering = false);
//
// private:
//     int _insert(const int x);
// };

} // namespace gc

#endif // __COLORING_HEURISTIC_HPP


