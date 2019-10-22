#include <iomanip>
#include <iomanip>
#include <iostream>

#include "graph.hpp"
#include "options.hpp"
#include "partition.hpp"
#include "statistics.hpp"

#include <minicsp/core/utils.hpp>

#ifndef __COLORING_HEURISTIC_HPP
#define __COLORING_HEURISTIC_HPP

// #define _DEBUG_DSATUR


namespace gc
{

struct coldomain {

    std::vector<int> b;
    size_t size_;

    inline size_t size() { return size_; }

    coldomain() {}
    coldomain(const int ub)
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

std::ostream& operator<<(std::ostream& os, const coldomain& x)
{
	x.display(os);
	return os;
}

std::ostream& operator<<(std::ostream& os, const coldomain* x)
{
	x->display(os);
	return os;
}



//..............{D+1}{-}(d+1)(-)x....(d)......(d-1)(d-2)(d-3).....(d-4)...(-)
//[-----------colored----------][---------------uncolored-------------------]

struct coloring_heuristic {

		// store the best full coloring
    std::vector<int> best;
		// number of distinct colors in 'best_coloring'
		size_t best_numcolors;

		// current (partial)  coloring
    std::vector<int> color;
		// number of distinct colors in 'color'
		size_t numcolors;
		
		// degree of the vertices in the graph
    std::vector<int> degree;

		// store info about the colors of the niehgbors in 'color'
    std::vector<coldomain> neighbor_colors;

		// ordered list of the vertices of the graph 
    std::vector<int> order;
		
		// rank of every vertex in 'order'
    std::vector<std::vector<int>::iterator> rank;
		
		// 'last_vertex[d]' points just after the last unassigned vertex of saturation degree d
    std::vector<std::vector<int>::iterator> last_vertex;

		// 1 assigning, -1 unassigning
		// int status;


		std::vector<int>::iterator beg_update;
		std::vector<int>::iterator end_update;


    std::mt19937 random_generator;


    // initialise the structures w. r. t. graph g and upper bound ub
    template <class graph_struct, class compare>
    void initialise(graph_struct& g, const int ub, compare comp) //=[&](const int, const int) { return 1; }
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
        for (auto it{begin(neighbor_colors)}; it != end(neighbor_colors);
             ++it) {
            it->initialise(ub);
        }
				// allocate as much as needed
        neighbor_colors.resize(g.capacity(), coldomain(ub));

				// the saturation degree is null for every vertex
        last_vertex.resize(ub + 2, begin(order));
				// hence last_vertex[0] oints at the end of order
        *begin(last_vertex) = end(order);

				// status = 1;
    }
		
		
		
    template <class graph_struct>
    void assign_color(graph_struct& g, const int x, const int c)
    {			
        color[x] = c;

        // update the saturation degree of x's neighbors
        for (auto y : g.matrix[x])
            // if (rank[y] > rank[x])
					// if(color[y] < 0)
					if(rank[y] >= beg_update and rank[y] < end_update)
                if (neighbor_colors[y].add(c)) 
										move_up(y, neighbor_colors[y].size());
    }

    template <class graph_struct>
    void unassign_color(graph_struct& g, const int x)
    {
			
			// std::cout << "unassign " << x << std::endl;
			
			auto c{color[x]};
        color[x] = -1;
				
        // update the saturation degree of x's neighbors
        for (auto y : g.matrix[x])
            // if (rank[y] > rank[x])
					// if(color[y] >= 0)
					if(rank[y] >= beg_update and rank[y] < end_update)
                if (neighbor_colors[y].remove(c)) 
                        move_down(y, neighbor_colors[y].size() + 1);
    }
		
		
		
    template <class graph_struct>
    void reassign_color(graph_struct& g, const int x, const int c)
    {
			
			std::cout << "reassign " << x << std::endl;
			
			auto oc{color[x]};
        color[x] = c;
				
        // update the saturation degree of x's neighbors
        for (auto y : g.matrix[x])
            if (color[y] >= 0) {
                if (neighbor_colors[y].remove(oc)) 
                        move_down(y, neighbor_colors[y].size() + 1);
                if (neighbor_colors[y].add(c)) 
                        move_up(y, neighbor_colors[y].size() + 1);
								
							}
    }

    // template <class graph_struct>
    // void re_assign(graph_struct& g, const int x, const int c)
    // {
    //     // ++num_reassign;
    //     // auto o{color[x]};
    //
    //     // color_bag.move(x, o, c);
    //     unassign_color(g, x);
    //     assign_color(g, x, c);
    // }

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
				
				std::cout << "swap " << l << " and " << y << std::endl;

        rank[l] = rank[y];
        rank[y] = last_vertex[d];

        *rank[y] = y;
        *rank[l] = l;
    }
		


		
    template <class graph_struct, typename tiebreaker>
    int brelaz_greedy(graph_struct& g, const int ub,
        std::vector<int>::iterator start, const int limit, tiebreaker criterion
        )
    {
			beg_update = start;
			end_update = end(order);
			
        int potential_colors = begin(neighbor_colors)->b.size();

        int c, d;

        std::vector<int>::iterator candidate{start};
				
				

        while (candidate != end(order)) {

#ifdef _DEBUG_DSATUR
            check_consistency(g);
#endif
						
            // get the highest saturation degree
            d = neighbor_colors[*candidate].size();
						
            if (limit > 1) {
                auto best{std::min_element(candidate,
                    std::min(last_vertex[d], candidate + limit), criterion
                    )};

                std::swap(
                    rank[*best], rank[*candidate]); // not sure this is useful
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
    int dsatur(
        graph_struct& g, const int ub, const int limit = 1, const int seed = 1)
    {
        if (g.nodes.empty())
            return 0;
				
				random_generator.seed(seed);

        initialise(g, ub, [&](const int x_, const int y_) {
            return (degree[x_] > degree[y_]);
        });	

        return brelaz_greedy(g, ub, begin(order), limit,
            [&](int x, int y) { return degree[x] > degree[y]; });
    }
		
		// update the saturation degrees; reorder by increasing SD; fix the last_vertex pointers
    template <class graph_struct>
    void close(graph_struct& g)
    {
			
        // update the color neighborhood ()
        for (auto v : order) {
            for (auto u : g.matrix[v]) 
                if (rank[u] < rank[v]) 
                    neighbor_colors[u].add(color[v]);
            degree[v] = g.matrix[v].size();
        }

        std::sort(begin(order), end(order), [&](const int x_, const int y_) {
            return (neighbor_colors[x_].size() > neighbor_colors[y_].size()
                or (neighbor_colors[x_].size() == neighbor_colors[y_].size()
                       and degree[x_] < degree[y_]));
        });

        for (auto vptr{begin(order)}; vptr != end(order); ++vptr) {
            rank[*vptr] = vptr;
        }


        last_vertex.resize(numcolors+1);
				auto d{0};
        for (auto it{end(order)}; it != begin(order); --it) {
            auto v{*(it-1)};
            while (d < numcolors and neighbor_colors[v].size() >= d) {
                last_vertex[d++] = it;
            }
        }
				last_vertex.back() = begin(order);
				
				
				// status = -1;
	}
	
			
  template <class graph_struct>
  int degeneracy(graph_struct& g)
  {		
		beg_update = begin(order);
		end_update = end(order);
		
		int dgn{0};
		for(auto vptr{rbegin(order)}; vptr!=rend(order); ++vptr)
		{

			print(g);
			std::cout << std::endl;
				
				int d = neighbor_colors[*vptr].size();
				
				dgn = std::max(dgn, d);
				
        // move all the pointers >= d
        while (d >= 0)
            --last_vertex[d--];
				
				--end_update;
				unassign_color(g, *vptr);
			}

			print(g);
			std::cout << std::endl;
				
			return dgn;
    }
		
		
		// try to find a color for x within y's neighborhood
    template <class graph_struct>
    bool recolor(graph_struct& g, const int x, const int y)
    {
			for(auto c{0}; c<numcolors; ++c)
			if(!neighbor_colors[y].contain(c) and neighbor_colors[x].contain(c)) {
					reassign_color(g, x, c);
				break;	
				}
		}
		
		
    void clear()
    {
        last_vertex.clear();
        color.clear();
        if (neighbor_colors.size() > 0)
            for (auto v : order)
                neighbor_colors[v].clear();
        order.clear();
    }


    template <class graph_struct> void print(graph_struct& g)
    {

        std::cout << std::endl;
				
				
				
        // dd = neighbor_colors[*begin(order)].size()+1;
        //
        //  for (auto v{begin(order)}; v != end(order); ++v) {
        //
        //      if (v == last_vertex[dd]) {
        //          std::cout << "-------(" << dd << ")" << std::endl;
        //  								--dd;
        //      }
        //      std::cout << std::setw(3) << *v << " " << std::setw(3)
        //                << neighbor_colors[*v].size() << " " << std::setw(3)
        //                << degree[*v] << std::endl;
        //  }
				
				
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

            std::cout << "(" << neighbor_colors[v].size();

            std::cout << "|" << degree[v] << ") " ; //<< neighbor_colors[v];

            if (color[v] >= 0) {
                std::cout << " ** " << color[v] << " **";
            }

            std::cout << std::endl;
        }
    }

    template <class graph_struct> void check_consistency(graph_struct& g)
    {

        // assert(!full);

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


                for (auto c{0}; c < numcolors; ++c) {

                    if (colv[c] != 0) {

                        std::cout << "problem in coldomain of " << v << " ("
                                  << c << ") ";
                        neighbor_colors[v].display(std::cout);
                        std::cout << std::endl;

                        for (auto b{0}; b < numcolors; ++b) {
                            std::cout << " " << neighbor_colors[v].b[b];
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
                                  << "ERROR: " << u << ":=" << color[u]
                                  << " and " << v << ":=" << color[v]
                                  << std::endl;
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


};

} // namespace gc

#endif // __COLORING_HEURISTIC_HPP


