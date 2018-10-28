#include <iostream>

#include "graph.hpp"
#include "options.hpp"
#include "sparseset.hpp"

#ifndef __DSATUR_HPP
#define __DSATUR_HPP

#define _DEBUG_DSATUR

namespace gc
{

// enum class core_type : uint8_t {
//     all,
//     witness,
//     lower,
// };

struct colvector {

    std::vector<int> b;
    size_t size_;

    inline size_t size() { return size_; }

    colvector() {}
    colvector(const int ub)
        : b(ub, 0)
    {
        size_ = 0;
    }

    bool add(const int x)
    {
        ++b[x];
        if (b[x] == 1) {
            ++size_;
            return true;
        }
        return false;
    }
    bool remove(const int x)
    {
        --b[x];
        if (b[x] == 0) {
            --size_;
            return true;
        }
        return false;
    }
    // int get_first_allowed() const
    // {
    //     for (auto i{begin(b)}; i != end(b); ++i)
    //         if (*i == 0)
    //             return i - begin(b);
    //     return b.size();
    // }
    int get_first_allowed(const int prev = -1) const
    {
        for (auto i{begin(b) + prev + 1}; i != end(b); ++i)
            if (*i == 0)
                return i - begin(b);
        return b.size();
    }
    bool contain(const int elt) const { return b[elt] > 0; }
    bool num_neighbors_of_color(const int elt) const { return b[elt]; }
    void clear()
    {
        for (auto i{begin(b)}; i != end(b); ++i)
            *i = 0;
        size_ = 0;
    }

    void initialise(const int ub)
    {
        b.clear();
        size_ = 0;
        b.resize(ub, 0);
    }

    std::ostream& display(std::ostream& os) const
    {
        os << "[";
        for (auto i{begin(b)}; i != end(b); ++i)
            if (*i != 0)
                os << " " << (i - begin(b));
        os << " ]";
        return os;
    }
};

std::ostream& operator<<(std::ostream& os, const colvector& x)
{
    return x.display(os);
}

std::ostream& operator<<(std::ostream& os, const colvector* x)
{
    return x->display(os);
}

struct dsatur {

    std::vector<int> color;
    std::vector<int> degree;
    std::vector<int> order;
    std::vector<colvector> neighbor_colors;

    std::vector<std::vector<int>::iterator> rank;
    std::vector<std::vector<int>::iterator> last_vertex;

    std::random_device rd;
    int limit;

    std::vector<int> single;
    colvector colorbag;

    std::vector<int> ncolor;
    std::vector<int>::iterator frontier;

    std::vector<int> core;

    int numcolors;

    bool full{false};
    bool use_recolor{true};

    std::vector<size_t> bag_idx;
    std::vector<sparseset> color_bag;

    void swap_colors(const int a, const int b)
    {
        std::swap(color_bag[a].list_, color_bag[b].list_);
        std::swap(color_bag[a].size_, color_bag[b].size_);
    }

    std::vector<int> first_of_color;
    std::vector<int> col_bag;

    // std::vector<int> branch_color;
    gc::bitset visited;
    std::vector<int> prev;
    std::vector<int> stack;
    std::vector<int> trail;

    // return true is recoloring was succesful (the next color to
    // use is <
    // numcolors) and false otherwise. numcolors is set to the color
    // to be used,
    // and the first vertex in the order list might change
    template <class graph_struct>
    bool recolor(graph_struct& g, const int x, int& col)
    {
        // std::cout << "try to avoid " << x << " @" << (rank[x] - begin(order))
        // << " <- " << col << std::endl;

        auto rx{rank[x]};

        // count the number of x's neighbors in each color bag
        single.clear();
        single.resize(col, -1);
        for (auto y : g.matrix[x]) {
            if (color[y] >= 0 and color[y] < col - 1) {
                if (single[color[y]] == -1)
                    single[color[y]] = y;
                else
                    single[color[y]] = -2;
            }
        }

        // find a color bag where x has only one neighbor
        for (int b{0}; b < col - 1; ++b) {
            if (single[b] >= 0) {
                auto w{single[b]};
                // w is the only neighbor of x colored with c

                colorbag.initialise(col);

                // make sur that we do not recolor vertices to
                // smaller colors,
                // as it could entail infinite loops
                for (int c = 0; c <= b; ++c)
                    colorbag.add(c);

                for (auto y : g.matrix[w])
                    if (color[y] >= 0) {

                        // std::cout << color[y] << " in N(" << w <<
                        // ")\n";

                        colorbag.add(color[y]);
                        if (colorbag.size() == col)
                            break;
                    }

                if (colorbag.size() < col) {
                    auto a{colorbag.get_first_allowed()};
                    // color w with a instead of b and x with b

                    unassign_color(g, w, b);
                    assign_color(g, w, a);
//                     // perhaps we should swap x and w so that
// std::swap(*(rx-1), *rank[w]);
// std::swap(rank[*(rx-1)], rank[w]);

// std::cout << " recolor " << w << " <- " << a << " to free up " << b <<
// std::endl;
// for(auto z{begin(order)}; z!=rx; ++z) {
// 	std::cout << std::std::setw(3) << *z << " " << std::std::setw(3) <<
// color[*z] << "
// [" ;
// 	for(auto colz{0}; colz<=col; ++colz) {
// 		if(neighbor_colors[*z].contain(colz)) {
// 			std::cout << " " << colz;
// 		}
// 	}
// 	std::cout << "]\n";
// }

#ifdef _DEBUG_DSATUR
                    check_consistency(g);
#endif

                    auto y{*rx};
                    if (y != x) {
                        // recoloring the neighbors of x decreased its
                        // sat-degree, but increased that of y (or it was
                        // already high)

                        if (neighbor_colors[y].size() == col) {
                            // y's sat-degree is not satisfying, let's try to
                            // recolor it as well
                            return recolor(g, y, col);
                        } else {
                            // ok, but we switch to assigning y before x
                            col = neighbor_colors[y].get_first_allowed();
                        }

                    } else {
                        // everything went as planned
                        col = b;
                    }

                    return true;
                }
            }
        }

        return false;
    }

    template <class graph_struct, class RandomIt>
    int brelaz_color_guided(graph_struct& g, const int ub, RandomIt beg,
        RandomIt stop, std::vector<int>& coloring, const int limit = 1,
        const int seed = 1)
    {

        if (beg == stop)
            return 0;

        brelaz_init(g, ub, limit, seed);

        auto first{begin(order)};

        for (auto vptr{begin(order)}; vptr != end(order); ++vptr)
            rank[*vptr] = vptr;

        std::vector<int> color_map(ub, -1);

        for (auto it{beg}; it != stop; ++it) {

            auto v{*it};
            std::swap(*first, *(rank[v]));
            rank[*(rank[v])] = rank[v];
            rank[v] = first;

            auto c{coloring[v]};

            // std::cout << " " << v << ":" << c;

            // std::cout << c << " " << color_map.size() << " " << ub <<
            // std::endl;
            assert(c < color_map.size());

            if (color_map[c] < 0) {
                color_map[c] = numcolors++;
            }
            c = color_map[c];

            ncolor.push_back(numcolors);

            color[v] = c;

            // update the saturation degree of x's neighbors
            for (auto u : g.matrix[v]) {
                if (color[u] < 0)
                    neighbor_colors[u].add(c);
                --degree[u];
            }

            ++first;

        }

        // std::cout << "GUIDE";
        // for (auto it{beg}; it != stop; ++it) {
        //     std::cout << " " << std::setw(3) << (*it);
        // }
        // std::cout << "\nCOLOR";
        // for (auto it{beg}; it != stop; ++it) {
        //     std::cout << " " << std::setw(3) << color_map[coloring[*it]];
        // }
        // std::cout << "\n";

        std::sort(first, end(order), [&](const int x_, const int y_) {
            return (neighbor_colors[x_].size() > neighbor_colors[y_].size()
                or (neighbor_colors[x_].size() == neighbor_colors[y_].size()
                       and degree[x_] > degree[y_]));
        });
        for (auto vptr{begin(order)}; vptr != end(order); ++vptr)
            rank[*vptr] = vptr;

        last_vertex.clear();
        last_vertex.resize(ub + 1, first);
        *begin(last_vertex) = end(order);

        int d{1};
        for (auto it{end(order)}; it-- != first;) {
            auto v{*it};
            while (neighbor_colors[v].size() >= d)
                last_vertex[d++] = it + 1;
        }

        auto ncol{brelaz_greedy(g, ub, first, limit)};

        for (auto v : g.nodes) {
            coloring[v] = color[v];
        }

        return ncol;
    }

    //
    template <class graph_struct>
    void brelaz_init(
        graph_struct& g, const int ub, const int limit, const int seed)
    {
        numcolors = 0;
        rank.resize(g.capacity());
        color.resize(g.capacity(), -1);
        degree.resize(g.capacity());
        // ncolor.reserve(g.size());

        // nodes are stored by non-increasing saturation
        // degree
        for (auto v : g.nodes) {
            order.push_back(v);
            degree[v] = g.matrix[v].size();
        }
        if (seed > 0) {
            // std::mt19937 s(rd());
            std::mt19937 s((seed * 256) + 75);
            std::shuffle(begin(order), end(order), s);
        }

        neighbor_colors.resize(g.capacity(), colvector(ub));

        last_vertex.resize(ub + 1, begin(order));
        *begin(last_vertex) = end(order);
    }

    template <class graph_struct>
    int brelaz_greedy(graph_struct& g, const int ub,
        std::vector<int>::iterator start, const int limit)
    {
        std::cout << "DSATUR\n";

        int c, d;

        for (auto vptr{start}; vptr != end(order); ++vptr)
            rank[*vptr] = vptr;

        std::vector<int>::iterator candidate{start};

        while (candidate != end(order)) {

#ifdef _DEBUG_DSATUR
            check_consistency(g);
#endif
						
            // get the highest saturation degree
            d = neighbor_colors[*candidate].size();

            if (limit > 1) {
                auto best{std::max_element(candidate,
                    std::min(last_vertex[d], candidate + limit),
                    [&](const int x_, const int y_) {
                        return (degree[x_] < degree[y_]);
                    })};
                std::swap(
                    rank[*best], rank[*candidate]); // not sure this is useful
                std::swap(*best, *candidate);
            }

            // use the first possible color for x
            c = neighbor_colors[*candidate].get_first_allowed();

            if (c == numcolors) {
                if (!use_recolor or (last_vertex[d] - last_vertex[d + 1]) > 2
                    or !recolor(g, *candidate, c)) {
                    ++numcolors;
                    frontier = candidate;
                } else {
                    --d;
                }
            }

            ncolor.push_back(numcolors);

            if (numcolors > ub) {
                color[*candidate] = c;
                return g.size();
            }

            // move all the pointers >= d
            while (++d < last_vertex.size())
                ++last_vertex[d];
						
            assign_color(g, *candidate, c);
            // std::cout << "vertex " << *candidate << " <- " << c << std::endl;

            // update degrees
            if (limit > 1)
                for (auto y : g.matrix[*candidate])
                    --degree[y];

            ++candidate;
        }

        //        for (auto it{rbegin(order)}; it != rend(order); ++it) {
        //             auto v{*it};
        //
        //
        // std::cout << "update " << v << std::endl;
        //
        //             // for (auto u : g.matrix[v]) {
        //             //     if (rank[u] < rank[v]) {
        //             //         neighbor_colors[u].add(color[v]);
        //             //     }
        //             // }
        //         }

        return numcolors;
    }

    //     template <class graph_struct>
    //     int brelaz_greedy(graph_struct& g, const int ub,
    //         std::vector<int>::iterator start, const int limit)
    //     {
    //
    //         int c, d;
    //
    //         for (auto vptr{start}; vptr != end(order); ++vptr)
    //             rank[*vptr] = vptr;
    //
    //         std::vector<int>::iterator candidate{start};
    //
    //         while (candidate != end(order)) {
    //
    // 						check_full_consistency(g);
    // #ifdef _DEBUG_DSATUR
    //             check_consistency(g);
    // #endif
    //
    // 						assert(color[*candidate] < 0);
    //
    //             // get the highest saturation degree
    //             d = neighbor_colors[*candidate].size();
    //
    //             if (limit > 1) {
    //                 auto best{std::max_element(candidate,
    //                     std::min(last_vertex[d], candidate + limit),
    //                     [&](const int x_, const int y_) {
    //                         return (degree[x_] < degree[y_]);
    //                     })};
    //                 std::swap(
    //                     rank[*best], rank[*candidate]); // not sure this is
    //                     useful
    //                 std::swap(*best, *candidate);
    //             }
    //
    //             // use the first possible color for x
    //             c = neighbor_colors[*candidate].get_first_allowed();
    //
    //             if (c == numcolors) {
    //                 if (!use_recolor or (last_vertex[d] - last_vertex[d + 1])
    //                 > 2
    //                     or !recolor(g, *candidate, c)) {
    //                     ++numcolors;
    //                     frontier = candidate;
    //                 } else {
    //                     --d;
    //                 }
    //             }
    //
    //             ncolor.push_back(numcolors);
    //
    //             if (numcolors > ub) {
    //                 color[*candidate] = c;
    //                 return g.size();
    //             }
    //
    //             // // move all the pointers >= d
    //             // while (++d < last_vertex.size())
    //             //     ++last_vertex[d];
    //
    //
    // 						assert(color[*candidate] < 0);
    //
    //             assign_color(g, *candidate, c);
    // 						std::cout << "vertex " <<
    // *candidate
    // <<
    // "
    // <-
    // "
    // <<
    // c
    // <<
    // std::endl;
    //
    //             // update degrees
    //             if (limit > 1)
    //                 for (auto y : g.matrix[*candidate])
    //                     --degree[y];
    //
    //             ++candidate;
    //         }
    //
    //         return numcolors;
    //     }

    template <class graph_struct>
    int brelaz_color(
        graph_struct& g, const int ub, const int limit = 1, const int seed = 1)
    {
        if (g.nodes.empty())
            return 0;

        brelaz_init(g, ub, limit, seed);

        std::sort(begin(order), end(order), [&](const int x_, const int y_) {
            return (degree[x_] > degree[y_]);
        });

        return brelaz_greedy(g, ub, begin(order), limit);
    }

    template <class graph_struct>
    void get_core(graph_struct& g, const gc::options::core_type t, const int lb,
        const int ub)
    {

        gc::bitset visited_vertex(0, g.capacity() - 1, gc::bitset::empt);
        // gc::bitset visited_color(0, numcolors - 1, gc::bitset::empt);
        std::vector<int> color_witness(numcolors, -1);
        // std::vector<int> core;
        core.clear();
        if (t == gc::options::core_type::ALL) {
            copy(order.begin(), frontier + 1, back_inserter(core));
        } else if (t == gc::options::core_type::LB) {
            for (int i = 0; i < order.size(); ++i) {
                core.push_back(order[i]);
                if (ncolor[i] > lb)
                    break;
            }
        } else {

            core.reserve(g.size());
            core.push_back(*frontier);

            for (auto vi{begin(core)}; vi != end(core); ++vi) {
                auto v{*vi};
                auto r{rank[v]};
                auto c{ncolor[r - begin(order)]};

                assert(color[v] < c);
								assert(color[v] >= 0);

                // visited_color.clear();
                color_witness.clear();
                color_witness.resize(numcolors, -1);

                // std::cout << " reason for " << v << "
                // = " << color[v] << " (" <<
                // c
                //           << ") @" << (r -
                //           begin(order)) << "\n";

                // get a witness for each color that u
                // cannot take
                // prefer 1/ a witness already in the
                // core, 2/ a witness with lower
                // rank
                auto maxc{(t == gc::options::core_type::WITNESS ? c : color[v])};
                for (auto u : g.matrix[v])
                    if (rank[u] < r and color[u] < maxc) {
                        if (visited_vertex.fast_contain(u))
                            color_witness[color[u]] = u;
                        else if (color_witness[color[u]] < 0
                            or rank[u] < rank[color_witness[color[u]]]) {
                            color_witness[color[u]] = u;
                        }
                    }
                for (auto k{0}; k < c; ++k) {
                    auto u{color_witness[k]};
                    if (color_witness[k] > 0) {
                        if (!visited_vertex.fast_contain(u)) {
                            visited_vertex.fast_add(u);
                            core.push_back(u);
                            // std::cout << " -> " << u
                            // << std::endl;
                        }
                    }
                }
            }
        }

        // std::cout << "core size = " << core.size() << " / "
        //           << (frontier - begin(order) + 1) << std::endl;
    }

    void select()
    {
        auto first{begin(order)};
        for (auto it{begin(core)}; it != end(core); ++it) {
            auto v{*it};
            auto u{*first};
            std::swap(*(rank[v]), *first);
            rank[u] = rank[v];
            rank[v] = first++;
        }

        int x{0};
        for (auto it{begin(order)}; it != end(order); ++it) {
            assert(it == rank[*it]);
            assert(x >= core.size() or core[x] == *it);

            ++x;
        }
    }

    template <class graph_struct>
    void assign_color(graph_struct& g, const int x, const int c)
    {
        color[x] = c;

        // #ifdef _DEBUG_DSATUR
        //         std::cout << "assigns " << c << " to " << x << "\n";
        // #endif

        // update the saturation degree of x's neighbors
        for (auto y : g.matrix[x])
            if (color[y] < 0 or full) {
                if (neighbor_colors[y].add(c)) {
                    auto d{neighbor_colors[y].size()};
                    // move y one partition up in the saturation degree
                    // list
                    move_up(y, d);
                }
            }
    }

    template <class graph_struct>
    void unassign_color(graph_struct& g, const int x, const int c)
    {

        // #ifdef _DEBUG_DSATUR
        //         std::cout << "unassigns " << x << " (" << c << ")\n";
        // #endif

        color[x] = -1;

        // update the saturation degree of x's neighbors
        for (auto y : g.matrix[x])
            if (color[y] < 0 or full) {
                if (neighbor_colors[y].remove(c)) {
										// move y one partition down in the saturation degree
                    // list
                    move_down(y, neighbor_colors[y].size() + 1);
                }
            }
        // else
        //                 assert(color[y] != c);
    }

    template <class graph_struct>
    void re_assign(graph_struct& g, const int x, const int c)
    {
        color_bag[color[x]].remove(x);
        color_bag[c].add(x);

        unassign_color(g, x, color[x]);
        assign_color(g, x, c);
    }

    // satur[y] was d, and is now d+1
    void move_up(const int y, const int d)
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
    void move_down(const int y, const int d)
    {
        // swap y with *last_vertex[d]-1
        auto l{*(--last_vertex[d])};

        rank[l] = rank[y];
        rank[y] = last_vertex[d];

        *rank[y] = y;
        *rank[l] = l;
    }

    void clear()
    {
        last_vertex.clear();
        color.clear();
        for (auto v : order)
            neighbor_colors[v].clear();
        order.clear();
        ncolor.clear();

        std::cout << "CLEAR\n";
    }

    template <class graph_struct> void init_local_search(graph_struct& g)
    {
        full = true;

        if (use_recolor)
            for (auto v : order)
                neighbor_colors[v].clear();

        // update the color neighborhood ()
        for (auto it{rbegin(order)}; it != rend(order); ++it) {
            auto v{*it};
            for (auto u : g.matrix[v]) {
                if (rank[u] < rank[v] or use_recolor) {
                    neighbor_colors[u].add(color[v]);
                }
            }
            degree[v] = g.matrix[v].size();
        }

        std::sort(begin(order), end(order), [&](const int x_, const int y_) {
            return (neighbor_colors[x_].size() > neighbor_colors[y_].size()
                or (neighbor_colors[x_].size() == neighbor_colors[y_].size()
                       and degree[x_] > degree[y_]));
        });

        for (auto vptr{begin(order)}; vptr != end(order); ++vptr)
            rank[*vptr] = vptr;

        auto d{last_vertex.size()};
        for (auto it{begin(order)}; it != end(order); ++it) {
            auto v{*it};
            while (d and neighbor_colors[v].size() < d) {
                last_vertex[d--] = it;
            }
        }
        for (auto i{0}; i <= d; ++i) {
            last_vertex[i] = end(order);
        }

        visited.initialise(0, g.capacity() - 1, gc::bitset::empt);
        prev.resize(g.capacity(), -1);

        color_bag.resize(numcolors);
        // std::vector<size_t> idx(order.size(), sparseset::NOVAL);
        bag_idx.resize(order.size(), NOINDEX);
        // std::vector<sparseset> bags(numcolors);
        for (auto it{begin(color_bag)}; it != end(color_bag); ++it)
            it->binds(&bag_idx);

        for (auto it{begin(order)}; it != end(order); ++it) {
            bags[color[*it]].add(*it);
        }
    }

    template <class graph_struct> void percolate(graph_struct& g)
    {
        compute_color_bags();
        for (auto r{begin(col_bag)}; r != end(col_bag); ++r) {
            auto v{*r};
            auto c{color[v]};
            auto b{neighbor_colors[v].get_first_allowed()};
            if (b < c) {
                unassign_color(g, v, c);
                assign_color(g, v, b);
            }
        }
        compute_color_bags();
    }

    // explore randomly a path from x
    template <class graph_struct> bool randpath(graph_struct& g, const int x)
    {

        // check_full_consistency(g);

        single.clear();
        single.resize(numcolors, -1);
        for (auto y : g.matrix[x]) {
            if (single[color[y]] == -1)
                single[color[y]] = y;
            else
                single[color[y]] = -2;
        }
        single[color[x]] = -2;

        auto success{-1};
        stack.clear();
        for (int b{0}; b < numcolors and success < 0; ++b)
            // if (b == numcolors - 1 and single[b] == -1)
            if (single[b] >= 0)
                stack.push_back(single[b]);
            else if (single[b] == -1 and b < numcolors - 1)
                success = b;

        if (stack.size() == 0) {
            // auto y{x};

            // std::cout << " -> :(\n";

            while (trail.size() > 0) {

                auto y{trail.back()};
                trail.pop_back();

                // std::cout << "unroll " << y << " <- " << trail.back()
                //           << std::endl;

                re_assign(g, y, trail.back());
                trail.pop_back();
                // y = prev[y];
            }

            return false;

        } else if (success >= 0) {
            re_assign(g, x, success);
            // std::cout << " -> " << success << "!\n";
            trail.clear();
            return true;

        } else {

            if (!visited.fast_contain(x)) {
                trail.push_back(color[x]);
                trail.push_back(x);
                visited.fast_add(x);
            }

            auto r{rand() % stack.size()};

            auto y{stack[r]};
            auto c{color[y]};
            re_assign(g, x, c);
            // std::cout << " -> (" << r << "|" << stack.size() << ") " << y <<
            // ":"
            //           << c;

            prev[y] = x;

            return randpath(g, y);
        }
    }

    // find a path from the vertex in stack to a free vertex
    template <class graph_struct>
    bool findpath_out(graph_struct& g, const int col)
    {
        trail.clear();
        visited.clear();
        for (auto x : stack) {
            visited.fast_add(x);
            prev[x] = x;
        }

        while (stack.size() > 0) {
            auto x = stack.back();
            stack.pop_back();

            if (neighbor_colors[x].size() - neighbor_colors[x].contain(col)
                < numcolors - 2) {
                auto a{neighbor_colors[x].get_first_allowed()};
                while (a == col or a == color[x])
                    a = neighbor_colors[x].get_first_allowed(a);
                assert(a != col);
                assert(a != color[x]);

                if (x != prev[x])
                    re_assign(g, prev[x], color[x]);
                re_assign(g, x, a);

                return true;
            }

            // unroll that branch
            while (trail.size() > 0 and trail.back() != prev[prev[x]]) {
                auto p{trail.back()};
                trail.pop_back();
                auto c{trail.back()};
                trail.pop_back();

                // std::cout << "b";

                // std::cout << p << " <-b- " << c << std::endl;

                re_assign(g, p, c);
            }

            // std::cout << "(" << trail.size() << ")[";
            // for(int i=0; i<trail.size(); i+=2) {
            // 	std::cout << " " << trail[i] << ":" << trail[i+1] ;
            // }
            // std::cout << " ]\n";

            // px = x;

            // for (int i = 0; i < trail.size(); ++i)
            //     std::cout << " ";
            // std::cout << x << " (" << color[x] << "):"; //<< std::end;

            single.clear();
            single.resize(numcolors, -1);
            for (auto y : g.matrix[x]) {
                if (single[(color)[y]] == -1)
                    single[(color)[y]] = y;
                else
                    single[(color)[y]] = -2;
            }

            // find a color bag where x has only one neighbor
            auto backtrack{true};
            for (int b{0}; b < numcolors; ++b) {
                if (single[b] >= 0) {
                    auto w{single[b]};
                    // w is the only neighbor of x colored with c

                    // auto c{color[x]};
                    if (!visited.fast_contain(w)) {
                        visited.fast_add(w);
                        backtrack = false;
                        prev[w] = x;

                        // std::cout << " " << w << "(" << color[w] << ")";

                        assert(color[w] == b);

                        // if (neighbor_colors[w].b[b] > 0) {
                        //     std::cout << std::endl
                        //               << w << " " << neighbor_colors[w]
                        //               << std::endl;
                        // }

                        assert(neighbor_colors[w].b[b] == 0);

                        stack.push_back(w);
                    }
                }
            }
            // std::cout << std::endl;

            if (!backtrack and prev[x] != x) {
                // std::cout << prev[x] << " <-f- " << color[x] << std::endl;
                // 	<< " store " << prev[x] << " = " << color[prev[x]] <<
                // "(" << (trail.size() == 0) << " or " << (trail.size() ?
                // trail.back() != prev[x] : 0) << ")\n";

                if (trail.size() == 0
                    or trail.back()
                        != prev[x]) { // store only at the first change
                    trail.push_back(color[prev[x]]);
                    trail.push_back(prev[x]);
                }
                re_assign(g, prev[x], color[x]);
            }
            // std::cout << "trail.size() = " << trail.size() << std::endl;
        }

        while (trail.size() > 0) {
            auto p{trail.back()};
            trail.pop_back();
            auto c{trail.back()};
            trail.pop_back();

            // std::cout << "b";

            // std::cout << p << " <-b- " << c << std::endl;

            re_assign(g, p, c);
        }

        return false;
    }

    // find a path from the vertex in stack to a free vertex
    template <class graph_struct> bool findpath(graph_struct& g)
    {
        // branch_color.clear();
        // std::copy(color.begin(),color.end(),back_inserter(branch_color));

        // std::cout << std::endl;

        trail.clear();
        visited.clear();
        for (auto x : stack) {
            visited.fast_add(x);
            prev[x] = x;
        }

        // auto px{*begin(stack)};
        while (stack.size() > 0) {
            auto x = stack.back();
            stack.pop_back();

            // std::cout << " " << x << " " << prev[x] << " " << prev[prev[x]]
            // << " " << (trail.size() ? trail.back() : -1) << std::endl;

            // unroll that branch
            while (trail.size() > 0 and trail.back() != prev[prev[x]]) {
                auto p{trail.back()};
                trail.pop_back();
                auto c{trail.back()};
                trail.pop_back();

                // std::cout << "b";

                // std::cout << p << " <-b- " << c << std::endl;

                re_assign(g, p, c);
            }

            // std::cout << "(" << trail.size() << ")[";
            // for(int i=0; i<trail.size(); i+=2) {
            // 	std::cout << " " << trail[i] << ":" << trail[i+1] ;
            // }
            // std::cout << " ]\n";

            // px = x;

            // for (int i = 0; i < trail.size(); ++i)
            //     std::cout << " ";
            // std::cout << x << " (" << color[x] << "):"; //<< std::end;

            single.clear();
            single.resize(numcolors, -1);
            for (auto y : g.matrix[x]) {
                if (single[(color)[y]] == -1)
                    single[(color)[y]] = y;
                else
                    single[(color)[y]] = -2;
            }

            // find a color bag where x has only one neighbor
            auto backtrack{true};
            for (int b{0}; b < numcolors; ++b) {
                if (single[b] >= 0) {
                    auto w{single[b]};
                    // w is the only neighbor of x colored with c

                    // auto c{color[x]};
                    if (!visited.fast_contain(w)) {
                        visited.fast_add(w);
                        backtrack = false;
                        prev[w] = x;

                        // std::cout << " " << w << "(" << color[w] << ")";

                        assert(color[w] == b);

                        // if (neighbor_colors[w].b[b] > 0) {
                        //     std::cout << std::endl
                        //               << w << " " << neighbor_colors[w]
                        //               << std::endl;
                        // }

                        assert(neighbor_colors[w].b[b] == 0);

                        auto Nsize{neighbor_colors[w].size()};

                        // if (Nsize == numcolors - 2
                        //     and neighbor_colors[w].get_first_allowed(b)
                        //         != numcolors - 1) {
                        //     std::cout << std::endl
                        //               << w << " " << neighbor_colors[w]
                        //               << std::endl;
                        // }
                        //
                        // assert((Nsize == numcolors - 2) <=
                        // (neighbor_colors[w].get_first_allowed(b) == numcolors
                        // - 1));

                        if (Nsize > numcolors - 3
                                + neighbor_colors[w].contain(numcolors - 1)) {
                            stack.push_back(w);
                        } else {
                            auto a{neighbor_colors[w].get_first_allowed(b)};

                            assert(a != b);

                            assert(a < numcolors - 1);

                            // std::cout << " -> " << a << " whoohoo!\n";
                            //
                            // if (prev[w] != prev[prev[w]])
                            //     std::cout << prev[prev[w]] << " <-f- "
                            //               << color[prev[w]] << std::endl;
                            // std::cout << prev[w] << " <-f- " << b <<
                            // std::endl;
                            // std::cout << w << " <-f- " << a << std::endl;

                            if (prev[w] != prev[prev[w]])
                                re_assign(g, prev[prev[w]], color[prev[w]]);

                            re_assign(g, prev[w], b);

                            // whoohoo!
                            unassign_color(g, w, b);
                            assign_color(g, w, a);

                            return true;
                        }
                    }
                }
            }
            // std::cout << std::endl;

            if (!backtrack and prev[x] != x) {
                // std::cout << prev[x] << " <-f- " << color[x] << std::endl;
                // 	<< " store " << prev[x] << " = " << color[prev[x]] <<
                // "(" << (trail.size() == 0) << " or " << (trail.size() ?
                // trail.back() != prev[x] : 0) << ")\n";

                if (trail.size() == 0
                    or trail.back()
                        != prev[x]) { // store only at the first change
                    trail.push_back(color[prev[x]]);
                    trail.push_back(prev[x]);
                }
                re_assign(g, prev[x], color[x]);
            }
            // std::cout << "trail.size() = " << trail.size() << std::endl;
        }

        while (trail.size() > 0) {
            auto p{trail.back()};
            trail.pop_back();
            auto c{trail.back()};
            trail.pop_back();

            // std::cout << "b";

            // std::cout << p << " <-b- " << c << std::endl;

            re_assign(g, p, c);
        }

        return false;
    }

    template <class graph_struct> bool descent(graph_struct& g)
    {
        percolate(g);

        int nmax{
            static_cast<int>(col_bag.size()) - first_of_color[numcolors - 1]};
        int ith = first_of_color[numcolors - 1];

        int improvements{0};
        while (ith < col_bag.size()) {
            if (color[col_bag[ith++]] < numcolors - 1)
                continue;

            stack.clear();
            stack.push_back(col_bag[ith - 1]);
            improvements += findpath(g);
        }

        trail.clear();

        int rnmax = 0;
        for (auto it{begin(order)}; it != end(order); ++it) {
            rnmax += (color[*it] == numcolors - 1);
        }
        assert(rnmax == nmax - improvements);

        if (improvements == nmax) {
            --numcolors;
            return true;
        }

        return false;
    }

    template <class graph_struct> void descents(graph_struct& g)
    {
        auto improvement{true};
        while (improvement) {
            improvement = false;
            int c{numcolors - 1};
            while (c >= 0) {
                while (!color_bag[c].empty()) {
                    stack.clear();
                    stack.push_back(color_bag[c].back());
                    if (!findpath_out(g, c))
                        break;
                }
            }
            if (c < 0) {
                swap_colors(c, --numcolors);
                improvement = true;
            }
        }
    }

    template <class graph_struct>
    bool randomwalk(graph_struct& g, const int limit)
    {
        percolate(g);

        int nmax{
            static_cast<int>(col_bag.size()) - first_of_color[numcolors - 1]};

        // print_col();

        int rnmax = 0;
        for (auto it{begin(order)}; it != end(order); ++it) {
            rnmax += (color[*it] == numcolors - 1);
            // if (color[*it] == numcolors - 1)
            // std::cout << " " << *it ;
        }

        assert(nmax == rnmax);

        int iter{0};
        int improvements{0};
        int successes{0};
        while (iter < limit and improvements < nmax) {
            auto x{rand() % order.size()};

            int improve{color[x] == numcolors - 1};
            // assert(visited.empty());
            assert(trail.size() == 0);

            visited.clear();
            int success{randpath(g, x)};

            successes += success;
            improvements += (improve & success);
            ++iter;
        }

        rnmax = 0;
        for (auto it{begin(order)}; it != end(order); ++it) {
            rnmax += (color[*it] == numcolors - 1);
        }
        assert(rnmax == nmax - improvements);

        // std::cout << nmax << " - " << improvements << " = " << rnmax <<
        // std::endl;

        // print_col();

        if (improvements == nmax) {
            --numcolors;
            return true;
        }

        return false;
    }

    template <class graph_struct> void local_search(graph_struct& g)
    {

        // std::vector<size_t> idx(order.size(), sparseset::NOVAL);
        // std::vector<sparseset> bags(numcolors);
        // for(auto it{begin(bags)}; it!=end(bags); ++it)
        // 	it->binds(&idx);

        // std::vector<sparseset> bags;
        // sparseset s(idx, order.size());
        // bags.push_back(s);
        // for(int i=1; i<numcolors; ++i) {
        // 	sparseset t(s);
        // 	bags.push_back(t);
        // }

        // sparseset* bags = new sparseset[numcolors];
        // bags[i] = new sparseset(idx, order.size());

        // for(auto it{begin(order)}; it!=end(order); ++it) {
        // 	bags[color[*it]].add(*it);
        // }
        //
        // for(int i=0; i<numcolors; ++i) {
        // 	std::cout << bags[i] << std::endl;
        // }
        //
        //
        // exit(1);

        init_local_search(g);

        print_col();

        descents(g);

        exit(1);

        // int limit{100000};
        // while (limit--) {
        //
        //     // std::cout << "descent " << numcolors << std::endl;
        //     //
        //     // while (descent(g)) {
        //     //     std::cout << " -- " << numcolors << std::endl;
        //     //     print_col();
        //     // }
        //
        //     // std::cout << "randwalk " << numcolors << std::endl;
        //
        //     if (randomwalk(g, 1000)) {
        //         std::cout << " -- " << numcolors << std::endl;
        //         print_col();
        //     }
        //     // visited.clear();
        // }

        //
        //
        // percolate(g);
        //
        // check_full_consistency(g);
        //
        // print_col();
        //
        // percolate(g);
        //
        // // int prev=-1;
        // int ith = first_of_color[numcolors - 1];
        // // int x = 0;
        // // auto x{col_bag[first_of_color[numcolors - 1]]};
        //
        // do {
        //     check_full_consistency(g);
        //
        //     while (color[col_bag[ith]] < numcolors - 1)
        //         ++ith;
        //     if (ith >= col_bag.size())
        //         break;
        //
        //     stack.clear();
        //     stack.push_back(col_bag[ith]);
        // } while (findpath(g));
        //
        // check_full_consistency(g);
        //
        // percolate(g);
        //
        // print_col();
        //
        // ith = first_of_color[numcolors - 1];
        // // int x = 0;
        // // auto x{col_bag[first_of_color[numcolors - 1]]};
        //
        // do {
        //     check_full_consistency(g);
        //
        //     while (color[col_bag[ith]] < numcolors - 1)
        //         ++ith;
        //     if (ith >= col_bag.size())
        //         break;
        //
        //     stack.clear();
        //     stack.push_back(col_bag[ith]);
        // } while (findpath(g));
        //
        // check_full_consistency(g);
        //
        // percolate(g);
        //
        // print_col();
        //
        // ith = first_of_color[numcolors - 2];
        // // int x = 0;
        // // auto x{col_bag[first_of_color[numcolors - 1]]};
        //
        // do {
        //     check_full_consistency(g);
        //
        //     while (color[col_bag[ith]] < numcolors - 2)
        //         ++ith;
        //     if (ith >= first_of_color[numcolors - 1])
        //         break;
        //
        //     stack.clear();
        //     stack.push_back(col_bag[ith]);
        // } while (findpath(g));
        //
        // check_full_consistency(g);
        //
        // percolate(g);
        //
        // print_col();
        //
        // // std::mt19937 s(1234);
        //
        // int count{0};
        // int improvements{0};
        // int successes{0};
        // while (count < 1000) {
        //     auto x{rand() % order.size()};
        //
        //     int improve{color[x] == numcolors - 1};
        //
        //     std::cout << x << ":" << color[x];
        //
        //     prev[x] = x;
        //
        //     // assert(visited.empty());
        //     assert(trail.size() == 0);
        //
        //     visited.clear();
        //     int success{randpath(g, x)};
        //
        //     successes += success;
        //     improvements += (improve & success);
        //
        //     check_full_consistency(g);
        //
        //     ++count;
        // }
        //
        // std::cout << improvements << " " << successes << " " << count
        //           << std::endl;
        //
        // percolate(g);
        //
        // print_col();
        //
        // --numcolors;
        // ith = first_of_color[numcolors - 1];
        // // int x = 0;
        // // auto x{col_bag[first_of_color[numcolors - 1]]};
        //
        // do {
        //     check_full_consistency(g);
        //
        //     while (color[col_bag[ith]] < numcolors - 1)
        //         ++ith;
        //     if (ith >= col_bag.size())
        //         break;
        //
        //     stack.clear();
        //     stack.push_back(col_bag[ith]);
        // } while (findpath(g));
        //
        // check_full_consistency(g);
        //
        // percolate(g);
        //
        // print_col();
        // // auto d{rand() % (n * n)};
        // //
        // //
        // //         auto x{col_bag[col_bag.size() / 2]};
        // //
        // //         std::cout << x;
        // //
        // //         assert(trail.size() == 0);
        // //         randpath(g, x);
        // //
        // //         print_col();
        // //
        // //         x = col_bag[col_bag.size() / 3];
        // //
        // //         std::cout << x;
        // //         assert(trail.size() == 0);
        // //         randpath(g, x);
        // //
        // //         print_col();
        //
        // // findpath(g);
        //
        // //         int iter{100};
        // //         while (iter--) {
        // //             check_full_consistency(g);
        // //
        // //             std::cout << std::endl;
        // //
        // //             print_col();
        // //
        // //             std::cout << "pick the " << ith << "-th element of the
        // //             last bag: ";
        // //
        // //             // pick a vertex of the last color
        // //             while (
        // //                 ith < col_bag.size() and color[col_bag[ith]] <
        // //                 numcolors - 1) {
        // //                 std::cout << "(" << col_bag[ith] << ") ";
        // //                 ++ith;
        // //             }
        // //             if (ith < col_bag.size())
        // //                 x = col_bag[ith++];
        // //             else { // if there is none, decrease the number of
        // //             colors,
        // //                 // percolate
        // //                 // and continue
        // //
        // //                 std::cout << " improvement!n\n";
        // //                 print_col();
        // //                 std::cout << "pick: ";
        // //
        // //                 --numcolors;
        // //                 percolate(g);
        // //                 ith = first_of_color[numcolors - 1];
        // //                 x = col_bag[ith++];
        // //             }
        // //
        // //             std::cout << x << "\n";
        // //
        // //             // count the number of x's neighbors in each color bag
        // //             single.clear();
        // //             single.resize(numcolors, -1);
        // //             for (auto y : g.matrix[x]) {
        // //                 if (single[color[y]] == -1)
        // //                     single[color[y]] = y;
        // //                 else
        // //                     single[color[y]] = -2;
        // //             }
        // //
        // //             // find a color bag where x has only one neighbor
        // //             for (int b{0}; b < numcolors; ++b) {
        // //                 if (single[b] >= 0) {
        // //                     auto w{single[b]};
        // //                     // w is the only neighbor of x colored with c
        // //
        // //                     std::cout << x << " has a single neighbor in
        // bag
        // //                     " << b
        // //                               << ": " << w; //<< std::endl;
        // //
        // //                     if (neighbor_colors[w].size() < numcolors - 1)
        // {
        // //                         auto a{colorbag.get_first_allowed()};
        // //
        // //                         std::cout << " -> " << a << "\n";
        // //
        // //                         unassign_color(g, w, b);
        // //                         assign_color(g, w, a);
        // //
        // // #ifdef _DEBUG_DSATUR
        // //                         check_full_consistency(g);
        // // #endif
        // //
        // //                         unassign_color(g, x, numcolors - 1);
        // //                         assign_color(g, x, b);
        // //
        // // #ifdef _DEBUG_DSATUR
        // //                         check_full_consistency(g);
        // // #endif
        // //
        // //                         // NEED TO UPDATE THE COLOR BAGS
        // //                         // w goes from b to a, and x from
        // numcolors-1
        // //                         to
        // //                         // b
        // //
        // //                         break;
        // //
        // //                     }
        // //
        // //                     else {
        // //                         std::cout << " (too constrained: " <<
        // //                         neighbor_colors[w]
        // //                                   << ")\n";
        // //                     }
        // //                 }
        // //             }
        // //         }
    }

    void compute_color_bags()
    {
        first_of_color.clear();
        first_of_color.resize(last_vertex.size(), 0);
        col_bag.resize(order.size());

        for (auto r{begin(order)}; r != end(order); ++r) {
            auto v{*r};
            ++first_of_color[color[v]];
        }

        for (auto d{begin(first_of_color) + 1}; d != end(first_of_color); ++d) {
            *d += *(d - 1);
        }

        for (auto r{begin(order)}; r != end(order); ++r) {
            auto v{*r};
            col_bag[--first_of_color[color[v]]] = v;
        }

        // for(int c=0; c<numcolors; ++c) {
        // 	std::cout << " " << first_of_color[c];
        // }
        // std::cout << " " << col_bag.size() << std::endl;
    }

    void print_col()
    {
        for (auto it{begin(color_bag)}; it != end(color_bag); ++it)
            std::cout << *it << std::endl;

        // gc::bitset is(0, order.size() - 1, gc::bitset::empt);
        //
        // compute_color_bags();
        //
        // auto it{begin(col_bag)};
        // for (auto c{color[*begin(col_bag)]}; c <= color[*rbegin(col_bag)];
        //      ++c) {
        //     std::cout << c << ": ";
        //
        //     while (it != end(col_bag) and color[*it] == c) {
        //         is.add(*it);
        //         ++it;
        //     }
        //
        //     std::cout << is << std::endl;
        //     is.clear();
        // }
    }

    template <class graph_struct> void print(graph_struct& g)
    {

        std::cout << std::endl;
        int d = last_vertex.size() - 1;
        for (auto r{begin(order)}; r != end(order); ++r) {

            bool lim{false};
            while (last_vertex[d] == r) {
                lim = true;
                std::cout << "start[" << d - 1 << "] ";
                --d;
            }
            if (lim) {
                std::cout << std::endl;
            }
            auto v{*r};

            std::cout << std::setw(3) << v << ": ";

            std::cout << "(" << neighbor_colors[v].size() << ") "
                      << neighbor_colors[v];

            if (color[v] >= 0) {
                std::cout << " ** " << color[v] << " **";
            }

            std::cout << std::endl;

            // if (color[v] >= 0) {
            //      std::cout << color[v] << " | N(" << v << ")";
            //      for (auto u : g.matrix[v]) {
            //          if (color[u] < 0)
            //              std::cout << " " << u;
            //      }
            //      std::cout << std::endl;
            //  } else {
            //      std::cout << " (" << neighbor_colors[v].size() <<
            //      ")"
            //                << neighbor_colors[v] << std::endl;
            //  }
        }
    }

    template <class graph_struct> void check_consistency(graph_struct& g)
    {

        assert(!full);

        // print(g);
        for (auto r{begin(order)}; r != end(order); ++r) {
            assert(rank[*r] == r);
        }

        for (size_t d{last_vertex.size() - 1}; d > 0; --d) {
            assert(last_vertex[d] <= last_vertex[d - 1]);
            for (auto r{last_vertex[d]}; r != last_vertex[d - 1]; ++r) {

                if (color[*r] < 0 and neighbor_colors[*r].size() != (d - 1)) {

                    std::cout << *r << " has satur degree "
                              << neighbor_colors[*r].size()
                              << " but is in bucket " << (d - 1) << std::endl;
                }

                assert(color[*r] >= 0 or neighbor_colors[*r].size() == (d - 1));
            }
        }

        std::vector<int> colv(numcolors);
        for (auto r{begin(order)}; r != end(order); ++r) {
            auto v{*r};
            auto d{neighbor_colors[v].size()};

            if (color[v] < -1) {
                for (auto c{0}; c < numcolors; ++c) {
                    colv[c] = neighbor_colors[v].b[c];
                }
                for (auto u : g.matrix[v]) {
                    if (color[u] >= 0 and rank[u] < rank[v])
                        --colv[color[u]];
                }
                // for(auto c{0}; c<numcolors; ++c) {
                // std::cout << " " << colv[c];
                // }
                // std::cout << std::endl;

                for (auto c{0}; c < numcolors; ++c) {

                    if (colv[c] != 0) {

                        std::cout << "problem in colvector of " << v << " ("
                                  << c << ") " << neighbor_colors[v]
                                  << std::endl;

                        for (auto b{0}; b < numcolors; ++b) {
                            std::cout << " " << neighbor_colors[v].b[b];
                        }
                        std::cout << std::endl;

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
                        std::cout << "N(" << v << ") = " << g.matrix[v]
                                  << std::endl;
                        std::cout << "ERROR: " << u << ":=" << color[u]
                                  << " and " << v << ":=" << color[v]
                                  << std::endl;
                    }

                    assert(color[u] != color[v]);

                    if (color[u] < 0) {
                        if (!neighbor_colors[u].contain(color[v]))
                            std::cout << "ERROR: NC(" << u
                                      << ") = " << neighbor_colors[u] << " - c["
                                      << v << "] = " << color[v] << std::endl;

                        assert(neighbor_colors[u].contain(color[v]));
                    }
                }
            } else {
                assert(last_vertex[d] > r);
                assert(last_vertex[d + 1] <= r);
            }
        }
    }

    template <class graph_struct> void check_full_consistency(graph_struct& g)
    {
        // print(g);
        for (auto r{begin(order)}; r != end(order); ++r) {
            assert(rank[*r] == r);
        }

        for (size_t d{last_vertex.size() - 1}; d > 0; --d) {
            assert(last_vertex[d] <= last_vertex[d - 1]);
            for (auto r{last_vertex[d]}; r != last_vertex[d - 1]; ++r) {

                if (neighbor_colors[*r].size() != (d - 1)) {

                    std::cout << *r << " has satur degree "
                              << neighbor_colors[*r].size()
                              << " but is in bucket " << (d - 1) << std::endl;
                }

                assert(neighbor_colors[*r].size() == (d - 1));
            }
        }

        std::vector<int> colv(numcolors);
        for (auto r{begin(order)}; r != end(order); ++r) {
            auto v{*r};
            auto d{neighbor_colors[v].size()};

            for (auto c{0}; c < numcolors; ++c) {
                colv[c] = neighbor_colors[v].b[c];
            }
            for (auto u : g.matrix[v]) {
                --colv[color[u]];
            }
            // for (auto c{0}; c < numcolors; ++c) {
            //     std::cout << " " << colv[c];
            // }
            // std::cout << std::endl;

            for (auto c{0}; c < numcolors; ++c) {

                if (colv[c] != 0) {

                    std::cout << "problem in colvector of " << v << " (" << c
                              << ") " << neighbor_colors[v] << std::endl;

                    for (auto b{0}; b < numcolors; ++b) {
                        std::cout << " " << neighbor_colors[v].b[b];
                    }
                    std::cout << std::endl;

                    for (auto b{0}; b < numcolors; ++b) {
                        std::cout << " " << colv[b];
                    }
                    std::cout << std::endl;

                    for (auto u : g.matrix[v]) {
                        std::cout << " " << u << " " << color[u] << std::endl;
                    }

                    exit(1);
                }

                assert(colv[c] == 0);
            }

            if (color[v] >= 0) {
                for (auto u : g.matrix[v]) {

                    if (color[u] == color[v]) {
                        std::cout << "N(" << v << ") = " << g.matrix[v]
                                  << std::endl;
                        std::cout << "ERROR: " << u << ":=" << color[u]
                                  << " and " << v << ":=" << color[v]
                                  << std::endl;
                    }

                    assert(color[u] != color[v]);

                    if (color[u] < 0) {
                        if (!neighbor_colors[u].contain(color[v]))
                            std::cout << "ERROR: NC(" << u
                                      << ") = " << neighbor_colors[u] << " - c["
                                      << v << "] = " << color[v] << std::endl;

                        assert(neighbor_colors[u].contain(color[v]));
                    }
                }
            } else {

                if (last_vertex[d] <= r or last_vertex[d + 1] > r)
                    std::cout << v << " @" << (rank[v] - begin(order))
                              << " d=" << d << " start[" << d
                              << "]=" << (last_vertex[d + 1] - begin(order))
                              << " start[" << ((int)d - 1)
                              << "]=" << (last_vertex[d] - begin(order))
                              << std::endl;

                assert(last_vertex[d] > r);
                assert(last_vertex[d + 1] <= r);
            }
        }
    }
};

} // namespace gc

#endif // __DSATUR_HPP


