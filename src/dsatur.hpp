#include <iostream>

#include "graph.hpp"


#ifndef __DSATUR_HPP
#define __DSATUR_HPP


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
    int get_first_allowed() const
    {
        for (auto i{begin(b)}; i != end(b); ++i)
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

    // return true is recoloring was succesful (the next color to use is <
    // numcolors) and false otherwise. numcolors is set to the color to be used,
    // and the first vertex in the order list might change
    template <class graph_struct>
    bool recolor(graph_struct& g, const int x, int& numcolors)
    {		
        auto rx{rank[x]};

        // count the number of x's neighbors in each color bag
        single.clear();
        single.resize(numcolors, -1);
        for (auto y : g.matrix[x]) {
            if (color[y] >= 0 and color[y] < numcolors - 1) {
                if (single[color[y]] == -1)
                    single[color[y]] = y;
                else
                    single[color[y]] = -2;
            }
        }

        // find a color bag where x has only one neighbor
        for (int b{0}; b < numcolors - 1; ++b) {
            if (single[b] >= 0) {
                auto w{single[b]};
                // w is the only neighbor of x colored with c
						
                colorbag.initialise(numcolors);

                // make sur that we do not recolor vertices to smaller colors,
                // as it could entail infinite loops
                for (int c = 0; c <= b; ++c)
                    colorbag.add(c);

                for (auto y : g.matrix[w])
                    if (color[y] >= 0) {

                        // std::cout << color[y] << " in N(" << w << ")\n";

                        colorbag.add(color[y]);
                        if (colorbag.size() == numcolors)
                            break;
                    }

                if (colorbag.size() < numcolors) {
                    auto a{colorbag.get_first_allowed()};
                    // color w with a instead of b and x with b

                    unassign_color(g, w, b);
                    assign_color(g, w, a);

#ifdef _DEBUG_DSATUR
                    check_consistency(g);
#endif
										
                    auto y{*rx};
                    if (y != x) {
												// recoloring the neighbors of x decreased its sat-degree, but increased that of y (or it was already high)

                        if (neighbor_colors[y].size() == numcolors) {
													// y's sat-degree is not satisfying, let's try to recolor it as well
                            return recolor(g, y, numcolors);
                        } else {
													// ok, but we switch to assigning y before x
                            numcolors = neighbor_colors[y].get_first_allowed();
                        }

                    } else {						
											// everything went as planned
                        numcolors = b;
                    }

                    return true;
                }
            }
        }

        return false;
    }

    template <class graph_struct>
    int brelaz_color(
        graph_struct& g, const int ub, const int limit = 1, const int seed = 1)
    {
        if (g.nodes.empty())
            return 0;
				
        rank.resize(g.capacity());
        color.resize(g.capacity(), -1);
        degree.resize(g.capacity());
        // ncolor.reserve(g.size());

        int c, d, numcolors{0};

        // nodes are stored by non-increasing saturation
        // degree
        for (auto v : g.nodes) {
            order.push_back(v);
            degree[v] = g.matrix[v].size();
        }
        if (seed > 0) {
            std::mt19937 s(rd());
            std::shuffle(begin(order), end(order), s);
        }

        std::sort(begin(order), end(order), [&](const int x_, const int y_) {
            return (degree[x_] > degree[y_]);
        });

        for (auto vptr{begin(order)}; vptr != end(order); ++vptr)
            rank[*vptr] = vptr;

        neighbor_colors.resize(g.capacity(), colvector(ub));

        last_vertex.resize(ub + 1, begin(order));
        *begin(last_vertex) = end(order);

        std::vector<int>::iterator candidate{begin(order)};

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
								std::swap(rank[*best], rank[*candidate]); // not sure this is useful
                std::swap(*best, *candidate);
            }

            // use the first possible color for x
            c = neighbor_colors[*candidate].get_first_allowed();

            if (c == numcolors) {
                if ((last_vertex[d] - last_vertex[d + 1]) > 2
                    or !recolor(g, *candidate, c)) {
                    ++numcolors;
                    frontier = candidate;
                } else {
                    --d;
                }
            }

            ncolor.push_back(numcolors);

            if (numcolors > ub) {
                frontier = end(order);
                return g.size();
            }

            // move all the pointers >= d
            while (++d < last_vertex.size())
                ++last_vertex[d];
						
            assign_color(g, *candidate, c);

            // update degrees
            if (limit > 1)
                for (auto y : g.matrix[*candidate])
                    --degree[y];

            ++candidate;
        }

        std::cout << "END DSATUR " << numcolors << std::endl;
        std::cout << ncolor.size() << " " << order.size() << " " << g.size()
                  << std::endl;
        assert(ncolor.size() == order.size());
        for (int i = 0; i < ncolor.size(); ++i) {
            assert(ncolor[i] >= 0);
            assert(rank[order[i]] == begin(order) + i);
        }

        gc::bitset visited_vertex(0, g.capacity() - 1, gc::bitset::empt);
        gc::bitset visited_color(0, numcolors - 1, gc::bitset::empt);
        std::vector<int> core;
        core.reserve(g.size());
        core.push_back(*frontier);

        for (auto vi{begin(core)}; vi != end(core); ++vi) {
            auto v{*vi};
            auto r{rank[v]};
            auto c{ncolor[r - begin(order)]};
            visited_color.clear();

            std::cout << " reason for " << v << " = " << color[v] << " (" << c
                      << ")\n";

            for (auto u : g.matrix[v])
                if (rank[u] < r and color[u] < c
                    and !visited_color.fast_contain(color[u])) {
                    visited_color.fast_add(color[u]);

                    std::cout << "  - " << u << " (" << color[u] << ")";

                    if (!visited_vertex.fast_contain(u)) {
                        visited_vertex.fast_add(u);
												core.push_back(u);
                        std::cout << " (add to list)";
                    }

                    std::cout << std::endl;
                }
        }

        std::cout << "core size = " << core.size() << std::endl;

        return numcolors;
    }

    template <class graph_struct>
    void assign_color(graph_struct& g, const int x, const int c)
    {
        color[x] = c;

#ifdef _DEBUG_DSATUR
        std::cout << "assigns " << c << " to " << x << "\n";
#endif
				
        // update the saturation degree of x's neighbors
        for (auto y : g.matrix[x])
            if (color[y] < 0) {
                if (neighbor_colors[y].add(c)) {
                    auto d{neighbor_colors[y].size()};
                    // move y one partition up in the saturation degree
                    // list
                    move_up(y, d);
                }
            } // else neighbor_colors[y].add(c);
    }

    template <class graph_struct>
    void unassign_color(graph_struct& g, const int x, const int c)
    {
			
#ifdef _DEBUG_DSATUR
        std::cout << "unassigns " << x << " (" << c << ")\n";
#endif
				
        // update the saturation degree of x's neighbors
        for (auto y : g.matrix[x])
            if (color[y] < 0) {
                if (neighbor_colors[y].remove(c)) {
										// move y one partition down in the saturation degree
                    // list
                    move_down(y, neighbor_colors[y].size() + 1);
                }
            } else
                assert(color[y] != c);
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
    }

    template <class graph_struct> void print(graph_struct& g)
    {
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

            std::cout << v << ": ";
            if (color[v] >= 0) {
                std::cout << color[v] << " | N(" << v << ")";
                for (auto u : g.matrix[v]) {
                    if (color[u] < 0)
                        std::cout << " " << u;
                }
                std::cout << std::endl;
            } else {
                std::cout << " (" << neighbor_colors[v].size() << ")"
                          << neighbor_colors[v] << std::endl;
            }
        }
    }

    template <class graph_struct> void check_consistency(graph_struct& g)
    {
        print(g);
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

            for (auto r{begin(order)}; r != end(order); ++r) {
                auto v{*r};
                auto d{neighbor_colors[v].size()};

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
                                          << ") = " << neighbor_colors[u]
                                          << " - c[" << v << "] = " << color[v]
                                          << std::endl;

                            assert(neighbor_colors[u].contain(color[v]));
                        }
                    }
                } else {
                    assert(last_vertex[d] > r);
                    assert(last_vertex[d + 1] <= r);
                }
            }
    }
};

} // namespace gc

#endif // __DSATUR_HPP


