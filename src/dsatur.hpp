#include <iostream>

#include "graph.hpp"
#include "options.hpp"
#include "partition.hpp"

#ifndef __DSATUR_HPP
#define __DSATUR_HPP

// #define _DEBUG_DSATUR

namespace gc
{

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

    void resize(const size_t s) { b.resize(s, 0); }

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
    int num_neighbors_of_color(const int elt) const { return b[elt]; }
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
                os << " " << (i - begin(b)) << "(" << *i << ")";
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

    // std::vector<size_t> bag_idx;
    // std::vector<sparseset> color_bag;
    partition color_bag;

    // void swap_colors(const int a, const int b)
    // {
    //     std::swap(color_bag[a].list_, color_bag[b].list_);
    //     std::swap(color_bag[a].size_, color_bag[b].size_);
    // }

    std::vector<int> first_of_color;
    std::vector<int> col_bag;

    // std::vector<int> branch_color;
    gc::bitset visited;
    std::vector<int> prev;
    std::vector<int> stack;
    std::vector<int> trail;

    //     std::vector<int> last_update;
    //     std::vector<double> avg_dsatur;
    // double decay{.99};

    long long int num_reassign{0};

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

        // std::cout << "INIT (" << ub << ")\n";

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

        // std::cout << (ub + 1) << std::endl;
        last_vertex.resize(ub + 1, begin(order));
        *begin(last_vertex) = end(order);
    }

    template <class graph_struct>
    int brelaz_greedy(graph_struct& g, const int ub,
        std::vector<int>::iterator start, const int limit)
    {

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

                // std::cout << "CANNOT IMPROVE ON " << (ub + 1) << std::endl;

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

        return numcolors;
    }


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
    }

    template <class graph_struct>
    void re_assign(graph_struct& g, const int x, const int c)
    {
        ++num_reassign;
        auto o{color[x]};

        color_bag.move(x, o, c);

        // if(color_bag[o])

        unassign_color(g, x, o);
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

        // if (avg_dsatur.size() > 0) {
        //     avg_dsatur[y] += neighbor_colors[y].size();
        //     avg_dsatur[y] /= 2;
        // }
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

        // if (avg_dsatur.size() > 0) {
        //     avg_dsatur[y] += neighbor_colors[y].size();
        //     avg_dsatur[y] /= 2;
        // }
    }

    void clear()
    {
        last_vertex.clear();
        color.clear();
        for (auto v : order)
            neighbor_colors[v].clear();
        order.clear();
        ncolor.clear();
    }

    template <class graph_struct>
    void init_local_search(graph_struct& g, std::vector<int>& isol)
    {
        full = true;

        auto colored{true};
        for (auto it{rbegin(order)}; colored and it != rend(order); ++it)
            colored = (color[*it] >= 0);

				int ub{0};
        for (auto v : order) {
            color[v] = isol[v];
						ub = std::max(ub, color[v]);
        }
				++ub;
        if (use_recolor or !colored)
            for (auto v : order) {
								neighbor_colors[v].resize(ub);
                neighbor_colors[v].clear();
						}

        // update the color neighborhood ()
        for (auto it{rbegin(order)}; it != rend(order); ++it) {
            auto v{*it};
            for (auto u : g.matrix[v]) {
                if (rank[u] < rank[v] or use_recolor or !colored) {
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

        // for (auto v : order)
        // std::cout << v << " " << color[v] << " " << degree[v] <<
        // neighbor_colors[v] << std::endl;

        // last_update.resize(g.capacity(), 0);
        // avg_dsatur.resize(g.capacity());
        for (auto vptr{begin(order)}; vptr != end(order); ++vptr) {
            rank[*vptr] = vptr;
            // avg_dsatur[*vptr]
            //     = static_cast<double>(neighbor_colors[*vptr].size());
        }

        // std::cout << "size = " << last_vertex.size() << std::endl;

				last_vertex.resize(ub + 1);
        auto d{last_vertex.size() - 1};
        for (auto it{begin(order)}; it != end(order); ++it) {
            auto v{*it};
            // std::cout << "v = " << v << std::endl;
            while (d > 0 and neighbor_colors[v].size() < d) {

                // std::cout << "d = " << d << std::endl;

                last_vertex[d--] = it;
            }
        }
        for (auto i{0}; i <= d; ++i) {
            last_vertex[i] = end(order);
        }

        visited.initialise(0, g.capacity() - 1, gc::bitset::empt);
        prev.resize(g.capacity(), -1);

        assert(numcolors <= ub);
        // color_bag.resize(g.capacity(), ub);
        color_bag.resize(g.capacity(), numcolors);
        // std::vector<size_t> idx(order.size(), sparseset::NOVAL);
        // bag_idx.resize(order.size(), NOINDEX);
        // std::vector<sparseset> bags(numcolors);
        // for (auto it{begin(color_bag)}; it != end(color_bag); ++it)
        //     it->binds(&bag_idx);

        for (auto it{begin(order)}; it != end(order); ++it) {
            assert(color[*it] >= 0);
            color_bag.add_elt(*it, color[*it]);
        }

        // std::cout << color_bag << std::endl;
    }

    // template <class graph_struct> void percolate(graph_struct& g)
    // {
    //     compute_color_bags();
    //     for (auto r{begin(col_bag)}; r != end(col_bag); ++r) {
    //         auto v{*r};
    //         auto c{color[v]};
    //         auto b{neighbor_colors[v].get_first_allowed()};
    //         if (b < c) {
    //             unassign_color(g, v, c);
    //             assign_color(g, v, b);
    //         }
    //     }
    //     compute_color_bags();
    // }

    template <class graph_struct> bool dsat_move(graph_struct& g, const int limit)
    {
        int moves{0};
        auto progress{true}, improvement{false};
        int change;
				int iter{0};
        while (progress and iter < limit) {
            progress = false;
            for (auto xp{rbegin(order)}; xp != rend(order); ++xp) {
                auto x{*xp};
                int c{0};
								
								++iter;

                for (; c < color_bag.size(); ++c) {
                    if (c != color[x] and !neighbor_colors[x].contain(c)) {
                        change = 0;
                        for (auto y : g.matrix[x]) {
                            if (neighbor_colors[y].num_neighbors_of_color(
                                    color[x])
                                    <= 1
                                and neighbor_colors[y].contain(c)) {
                                ++change;
                            } else if (neighbor_colors[y]
                                           .num_neighbors_of_color(color[x])
                                    > 1
                                and !neighbor_colors[y].contain(c)) {
                                --change;
                            }
                        }
                        if (change > 0) {
                            break;
                        }
                    }
                }
                if (c < color_bag.size()) {
                    ++moves;
                    auto o{color[x]};
                    re_assign(g, x, c);
                    progress = true;

                    if (color_bag[o].empty()) {
                        remove_color(g, o);
                        improvement = true;
                    }

                    break;
                }
            }
        }
				
				// std::cout << moves << " / " << iter << std::endl;

        return improvement;
    }

    // explore randomly a path from x
    template <class graph_struct>
    bool randpath_out(
        graph_struct& g, const int x, const int taboo, const int depth_limit)
    {
        single.clear();
        single.resize(color_bag.size(), -1);
        for (auto y : g.matrix[x]) {
            if (single[color[y]] == -1)
                single[color[y]] = y;
            else
                single[color[y]] = -2;
        }
        single[color[x]] = -2;

        auto success{-1};
        stack.clear();
        for (int b{0}; b < color_bag.size() and success < 0; ++b)
            if (single[b] >= 0
                // and (!visited.fast_contain(single[b]) or prev[single[b]]
                // != x)
                )
                stack.push_back(single[b]);
            else if (single[b] == -1 and b != taboo)
                success = b;

        // assert(depth_limit > 0);

        if (stack.size() == 0 or depth_limit == 0) {
            while (trail.size() > 0) {
                auto y{trail.back()};
                trail.pop_back();
                re_assign(g, y, trail.back());
                trail.pop_back();
            }

            return false;

        } else if (success >= 0) {
            re_assign(g, x, success);
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

            prev[y] = x;

            return randpath_out(g, y, taboo, depth_limit - 1);
        }
    }

    // find a path from the vertex in stack to a free vertex
    template <class graph_struct>
    bool findpath_out(graph_struct& g, const int col)
    {

        assert(trail.size() == 0);
        assert(visited.empty());
        assert(stack.size() == 1);

        for (auto x : stack) {
            visited.fast_add(x);
            prev[x] = x;
        }

        while (stack.size() > 0) {
            auto x = stack.back();
            stack.pop_back();

            // std::cout << "backtrack\n";
            // unroll that branch
            while (trail.size() > 0 and trail.back() != prev[prev[x]]) {
                auto p{trail.back()};
                trail.pop_back();
                auto c{trail.back()};
                trail.pop_back();
                re_assign(g, p, c);
            }
            // std::cout << "ok\n";

            // std::cout << "[";
            // for(int i=0; i<trail.size(); i+=2)
            // 	std::cout << " " << trail[i+1] << ":" << trail[i];
            // std::cout << " ] - " << x << " |";
            // for(auto x_ : stack)
            // std::cout << " " << x_ ;
            // 	std::cout << std::endl;

            if (neighbor_colors[x].size() - neighbor_colors[x].contain(col)
                < color_bag.size() - 2) {
                // std::cout << "success\n";

                auto a{neighbor_colors[x].get_first_allowed()};
                while (a == col or a == color[x])
                    a = neighbor_colors[x].get_first_allowed(a);
                assert(a != col);
                assert(a != color[x]);
                re_assign(g, prev[x], color[x]);

                re_assign(g, x, a);
                trail.clear();
                visited.clear();
                // std::cout << "ok\n";

                return true;
            }

            // std::cout << "branches (1)\n";
            single.clear();
            single.resize(color_bag.size(), -1);
            for (auto y : g.matrix[x]) {
                if (single[(color)[y]] == -1)
                    single[(color)[y]] = y;
                else
                    single[(color)[y]] = -2;
            }
            // std::cout << "branches (2)\n";

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

                        assert(color[w] == b);
                        assert(neighbor_colors[w].b[b] == 0);

                        stack.push_back(w);

                        // std::cout << " + " << w;
                    }
                }
            }
            // std::cout << std::endl;

            if (!backtrack and prev[x] != x) {

                // std::cout << "assign predecessor of " << x << ": " <<
                // prev[x] << " <- " << color[x] << std::endl;

                // std::cout << prev[x] << " was " << color[prev[x]] <<
                // std::endl;
                if (trail.size() == 0
                    or trail.back()
                        != prev[x]) { // store only at the first change
                    trail.push_back(color[prev[x]]);
                    trail.push_back(prev[x]);
                }
                // std::cout << "here\n";

                re_assign(g, prev[x], color[x]);

                // std::cout << "ok\n";
            }
        }

        // std::cout << "fail\n";
        while (trail.size() > 0) {
            auto p{trail.back()};
            trail.pop_back();
            auto c{trail.back()};
            trail.pop_back();
            re_assign(g, p, c);
        }
        visited.clear();
        // std::cout << "ok\n";

        return false;
    }

    template <class graph_struct>
    void remove_color(graph_struct& g, const int c)
    {
        const std::vector<int>& bag(color_bag[color_bag.size() - 1]);
        while (!bag.empty()) {
            re_assign(g, bag.back(), c);
        }
        color_bag.remove(color_bag.size() - 1);
        --numcolors;
    }

    template <class graph_struct> bool descent(graph_struct& g, int& npath)
    {
        auto improvement{false};
        auto progress{true};

        visited.clear();
        while (progress) {
            progress = false;
            int c{color_bag.size() - 1};
            while (c >= 0) {

#ifdef _DEBUG_DSATUR
                check_full_consistency(g, "start findpath loop");
#endif

                while (!color_bag[c].empty()) {

                    stack.clear();
                    stack.push_back(color_bag[c].back());

                    if (!findpath_out(g, c))
                        break;

                    ++npath;

#ifdef _DEBUG_DSATUR
                    check_full_consistency(g, "successful findpath");
#endif
                }

                if (color_bag[c].empty()) {

#ifdef _DEBUG_DSATUR
                    check_full_consistency(g, "before rm color (d)");
#endif

                    remove_color(g, c);

#ifdef _DEBUG_DSATUR
                    check_full_consistency(g, "after rm color (d)");
#endif

                    improvement = progress = true;

                    assert(numcolors == color_bag.size());

                    // std::cout << "ub = " << color_bag.size() << "
                    // (descent)\n";
                    // std::cout << "\ndescent:\n" << color_bag <<
                    // std::endl;
                }

                --c;
            }
        }

        return improvement;
    }

    template <class graph_struct>
    bool randomwalk(
        graph_struct& g, const int limit, const int target, int& npath)
    {
        int iter{0};
        while (iter < limit) {
            auto x{(rand() % 2
                    ? order[rand() % order.size()]
                    : color_bag[target][rand() % color_bag[target].size()])};

            assert(trail.size() == 0);

            visited.clear();

#ifdef _DEBUG_DSATUR
            check_full_consistency(g, "before randpath ");
#endif

            randpath_out(g, x, target, 100);

#ifdef _DEBUG_DSATUR
            check_full_consistency(g, "after randpath");
#endif

            ++iter;

            if (color_bag[target].size() == 0) {

#ifdef _DEBUG_DSATUR
                check_full_consistency(g, "before rm color (rw)");
#endif

                remove_color(g, target);

#ifdef _DEBUG_DSATUR
                check_full_consistency(g, "after rm color (rw)");
#endif

                // std::cout << "ub = " << color_bag.size() << "
                // (descent)\n";
                // std::cout << "\nrandwalk:\n" << color_bag << std::endl;
                npath += iter;

                return true;
            }
        }

        npath += iter;

        return false;
    }

    template <class graph_struct>
    void local_search(
        graph_struct& g, std::vector<int>& isol, gc::statistics& stat)
    {
				// assert(g.size() == isol.size());


        assert(isol.size() == color.size());

        init_local_search(g, isol);

        assert(stat.best_ub == color_bag.size());

        int num_rw_iter = 100, num_iter = 100000;

        int num_rp{0};
        int num_fp{0};

        std::vector<int> dsat_order;

        for (int i = 0; stat.best_lb < stat.best_ub and i < num_iter; ++i) {

            if (descent(g, num_fp)) {
                stat.notify_ub(color_bag.size());
                stat.display(std::cout);

                isol = color;
            }
						

            auto t{rand() % color_bag.size()};

            if (randomwalk(g, num_rw_iter, t, num_rp)) {
                stat.notify_ub(color_bag.size());
                stat.display(std::cout);

                isol = color;
            }
						
						
            if (dsat_move(g, 1000)) {
                stat.notify_ub(color_bag.size());
                stat.display(std::cout);

                isol = color;
            }
						
						// if(i % 100 == 0) {
						// std::cout << "-- " << num_reassign << "\n";
						// }

            // if(i % 100 == 0) {
            // 	dsat_order = order;
            // 	std::sort(begin(dsat_order), end(dsat_order));
            // 	auto m{avg_dsatur[*std::min_element(begin(dsat_order),
            // end(dsat_order), [&](const int x, const int y) {return
            // avg_dsatur[x] < avg_dsatur[y];})]};
            // 	int count = 0;
            // 	// std::cout << m << std::endl;
            // 	for(auto x : dsat_order) {
            // 		if(++count == 25)
            // 			break;
            // 		std::cout << std::setw(4) <<
            // std::setprecision(1)
            // << (avg_dsatur[x] - m);
            // 	}
            // 	std::cout << std::endl;
            // }
        }
    }

    // void compute_color_bags()
    // {
    //     first_of_color.clear();
    //     first_of_color.resize(last_vertex.size(), 0);
    //     col_bag.resize(order.size());
    //
    //     for (auto r{begin(order)}; r != end(order); ++r) {
    //         auto v{*r};
    //         ++first_of_color[color[v]];
    //     }
    //
    //     for (auto d{begin(first_of_color) + 1}; d != end(first_of_color);
    //     ++d) {
    //         *d += *(d - 1);
    //     }
    //
    //     for (auto r{begin(order)}; r != end(order); ++r) {
    //         auto v{*r};
    //         col_bag[--first_of_color[color[v]]] = v;
    //     }
    // }

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

                // std::cout << d << "/" << last_vertex.size() << ": "
                // 	<< (last_vertex[d + 1] - begin(order)) << ".."
                // 		<< (r - begin(order)) << ".."
                // 		<< (last_vertex[d] - begin(order)) << "\n";

                assert(last_vertex[d] > r);
                assert(d + 1 == last_vertex.size() or last_vertex[d + 1] <= r);
            }
        }
    }

    template <class graph_struct>
    void check_full_consistency(graph_struct& g, const char* msg)
    {
        // print(g);
        for (auto r{begin(order)}; r != end(order); ++r) {
            assert(rank[*r] == r);
        }

        for (size_t d{last_vertex.size() - 1}; d > 0; --d) {
            assert(last_vertex[d] <= last_vertex[d - 1]);
            for (auto r{last_vertex[d]}; r != last_vertex[d - 1]; ++r) {

                if (neighbor_colors[*r].size() != (d - 1)) {

                    std::cout << msg << ": " << *r << " has satur degree "
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

                    std::cout << msg << ": "
                              << "problem in colvector of " << v << " (" << c
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
                        std::cout << msg << ":\n"
                                  << "N(" << v << ") = " << g.matrix[v]
                                  << std::endl;
                        std::cout << "ERROR: " << u << ":=" << color[u]
                                  << " and " << v << ":=" << color[v]
                                  << std::endl;
                    }

                    assert(color[u] != color[v]);

                    if (color[u] < 0) {
                        if (!neighbor_colors[u].contain(color[v]))
                            std::cout << msg << ": "
                                      << "ERROR: NC(" << u
                                      << ") = " << neighbor_colors[u] << " - c["
                                      << v << "] = " << color[v] << std::endl;

                        assert(neighbor_colors[u].contain(color[v]));
                    }
                }
            } else {

                if (last_vertex[d] <= r or last_vertex[d + 1] > r)
                    std::cout
                        << msg << ": " << v << " @" << (rank[v] - begin(order))
                        << " d=" << d << " start[" << d
                        << "]=" << (last_vertex[d + 1] - begin(order))
                        << " start[" << ((int)d - 1)
                        << "]=" << (last_vertex[d] - begin(order)) << std::endl;

                assert(last_vertex[d] > r);
                assert(last_vertex[d + 1] <= r);
            }
        }
    }
};

} // namespace gc

#endif // __DSATUR_HPP


