

#include "sparse_dynamic_graph.hpp"


// #define _VERIFY_MCGRAPH

using std::swap;

using namespace gc;


// remove d from the bag 'd'
void coloring::remove(const int y, const int d) {
	
	int idy = rank[y];
	
	// std::cout << "remove " << y << " @ " << idy << " from bag " << d << std::endl;
	
	
	// swap y with first[satur[y].size-1], and increment first[satur[y].size-1]
	int idf = first[d];
	int f = order[idf];
	
	
	
	// std::cout << " -> swap with " << f << " @ " << idf << std::endl;
	
	rank[y] = idf;
	rank[f] = idy;
	
	order[idf] = y;
	order[idy] = f; 

	++first[d];
	
	
	// exit(1);
}


void coloring::brelaz_color(dyngraph& g) {	
		if(g.node.empty()) return;

		// first order the nodes by degree
		// degree.resize(g.node.size());
		for(auto v : g.node) {
				order.push_back(v);
				// degree[v] = degree(v);
		}

    // std::sort(begin(order), end(order),
    //     [&](const int x, const int y) { return (g.degree(x) > g.degree(y)); });

		int iter = 0;
		rank.resize(g.capacity);
		for(auto v : order) {
			rank[v] = iter++;
			// std::cout << v << " " << g.degree(v) << std::endl;
		}


		iter = 0;	
		int c, d;
		int x; // = order[iter]; // the next vertex to color
		// int num_colors = 0; // the number of colors used so far

		color.resize(g.capacity); 
		satur.resize(g.capacity); 
		first.push_back(0); // every node has saturation degree 0


		
		do {
			if(g.node.empty()) break;
			
			
			
			d = satur[order[iter]].size;
			
			while(first.size() > d + 1) {
				first.pop_back();
			}
						
			int i, prev;
			//
			//
			// for(auto v : order) {
			// 	if(rank[v] == iter)
			// 		std::cout << "| ";
			// 		std::cout << std::setw(3) << v << " ";
			// }
			// std::cout << std::endl;
			// for(auto v : order) {
			// 	if(rank[v] == iter)
			// 		std::cout << "| ";
			// 		std::cout << std::setw(3) << satur[v].size << " ";
			// }
			// std::cout << std::endl;
			// for(auto v : order) {
			// 	if(rank[v] == iter)
			// 		std::cout << "| ";
			// 		std::cout << std::setw(3) << g.degree(v) << " ";
			// }
			// std::cout << std::endl;
			// for(auto r : first) {
			// 	std::cout << std::setw(3) << order[r] << " ";
			// }
			// std::cout << std::endl;
			// std::cout << order[iter] << std::endl;
			
			
			
			
			i = 0;
			for(auto v : order) {
				assert(rank[v] == i++);
			}
			prev = -1;
			for(auto v : order) {
				assert(rank[v] <= iter or satur[v].size <= satur[prev].size);
				prev = v;
			}
			// prev = -1;
			// for(auto v : order) {
			// 	assert(prev < 0 or satur[v].size < satur[prev].size or g.degree(v) <= g.degree(prev));
			// 	prev = v;
			// }
			i=0;
			for(auto r : first) {
				assert(r>=0);
				assert(r<=order.size());
				if(r<order.size()) {
					// std::cout << "assert the satur degree of the first elt of bag " << i << " (" << order[r] << ") [or first["<< i-1 << "] == " << order[r] << "]\n";
					// std::cout << "r = " << r << " iter = " << iter << std::endl;
					
					if(r>=iter) {
					assert(satur[order[r]].size == i or (i>0 and first[i-1] == r));
					}
						// std::cout << "assert that satur degree of the predecessor " << order[r-1] << " of " << order[r] << " is " << (i+1) << " or empty " << std::endl;

						
					if(r>iter) {	
						assert(satur[order[r-1]].size == i + 1 or (first.size() > i + 1 and first[i + 1] == r));
					}
				}				
				++i;
			}
			

			
			
			
			// std::cout << "d(" << order[iter] << ") = " << d << std::endl;
			//
			// if(d > 0)
			// 	std::cout << "search in ["<< *(begin(order)+first[d]) << ".." <<  *(begin(order)+first[d-1]) << "[\n";
			// else
			// 	std::cout << "search from " << order[iter] << std::endl;
			
			
			d = satur[order[iter++]].size;
			
			x = *std::max_element(begin(order)+first[d], (d > 0 ? begin(order)+first[d-1] : end(order)), [&](const int x_, const int y_) { return (g.degree(x_) < g.degree(y_)); });
			
			
			// std::cout << "assert that satur[" << x << "] = " << d << std::endl;
			assert( satur[x].size == d );
			
			remove(x, d);
		
				c = satur[x].get();
				
				// std::cout << std::endl << "color " << x << " " << satur[x] << " with " << c << std::endl;
				// std::cout << std::endl;
				
				
				color[x] = c;
				// if(c >= num_colors)
				// 		++num_colors;
				
				
				assert(g.node.contain(x));
				assert(g.nodeset.contain(x));
				
				g.rem_node(x);
				
				
				// for(auto v : order) {
				// 	if(rank[v] == iter)
				// 		std::cout << "| ";
				// 		std::cout << std::setw(3) << v << " ";
				// }
				// std::cout << std::endl;
				// for(auto v : order) {
				// 	if(rank[v] == iter)
				// 		std::cout << "| ";
				// 		std::cout << std::setw(3) << satur[v].size << " ";
				// }
				// std::cout << std::endl;
				// for(auto v : order) {
				// 	if(rank[v] == iter)
				// 		std::cout << "| ";
				// 		std::cout << std::setw(3) << g.degree(v) << " ";
				// }
				// std::cout << std::endl;
				// for(auto r : first) {
				// 	std::cout << std::setw(3) << order[r] << " ";
				// }
				// std::cout << std::endl;
				// std::cout << order[iter] << std::endl;
				//
				// std::cout << "update saturation degree:\n";
				
				
				for(auto y : g.neighbor[x])
						if(satur[y].add(c)) {
								if(first.size() <= satur[y].size) {
										assert(first.size() == satur[y].size);
										first.push_back(*rbegin(first));
								}
								
								d = satur[y].size-1;
								
								remove(y,d);

								// int i = rank[y];
								// int l = first[satur[y].size];
								// assert(l > iter);
								//
								// std::cout << "re-rank " << y << " in [" << l << ".." << i << "]\n";
								//
								// // reorder y w.r.t. its degree
								// while(i > l && g.degree(y) > g.degree(order[i-1])) {
								// 		// std::cout << " -> " << (i-1) << std::endl;
								//
								// 		order[i] = order[i-1];
								// 		rank[order[i]] = i;
								// 		order[i-1] = y;
								// 		rank[y] = --i;
								//
								// }
								
						}
						

						// // clean the "first" pointer struct
						// while(first.back() < iter) {
						// 	first.pop_back();
						// }
						// while(first.back() == iter and first.size() > 0 and first.back() == *(rbegin(first)-1)) {
						// 	first.pop_back();
						// }
						
						// while(first.back() < iter) {
						// 	first.pop_back();
						// }
				
				// x = order[++iter];
				// ++first[satur[x].size];

		} while( true );
		


		// std::cout << "num_colors = " << num_colors << std::endl;

}



dyngraph::dyngraph(const int n)
{
	  capacity = n;

	  node.reserve(capacity);
		node.fill();

	  nodeset.initialise(0, capacity - 1, gc::bitset::full);

	  neighbor.resize(capacity); 
	  nb_index.resize(capacity); 

	  num_edges = 0;
}

dyngraph::dyngraph(const dyngraph& g) { 
	  capacity = g.capacity;

	  nodeset.initialise(0, capacity - 1, gc::bitset::empt);
		for(auto v : g.node) {
			declare_node(v);
		}

	  neighbor.resize(capacity); 
	  nb_index.resize(capacity); 

	  num_edges = g.num_edges;
	  nodeset.copy(g.nodeset);
	  for (unsigned i = 0; i < g.ranks.size(); ++i) {
	      edges.push_back(g.edges[i]);
	      ranks.push_back(g.ranks[i]);
	  }

	  for (int x = 0; x < capacity; ++x) {

	      for (auto it = begin(g.neighbor[x]); it != end(g.neighbor[x]); ++it)
	          neighbor[x].push_back(*it);

	      for (auto it = begin(g.nb_index[x]); it != end(g.nb_index[x]); ++it)
	          nb_index[x].push_back(*it);
	  }
}

int dyngraph::size() const { return node.size(); }

bool dyngraph::null() const { return node.size() == 0; }

bool dyngraph::empty() const { return num_edges == 0; }

bool dyngraph::full() const
{
    return node.size() * (node.size() - 1) == 2 * num_edges;
}

double dyngraph::get_density() const
{
    return (double)num_edges / (double)(size() * (size() - 1) / 2);
}

void dyngraph::sort(bool non_decreasing)
{
    verify("before sort");

    std::vector<int> sorted;
    for (int i = 0; i < size(); ++i) {
        sorted.push_back(node[i]);
    }

    std::sort(sorted.begin(), sorted.end(),
        [&](const int x, const int y) { return (degree(x) > degree(y))^non_decreasing; });

    std::vector<int> srank(capacity);
    for (int i = 0; i < size(); ++i) {
        srank[sorted[i]] = i;
    }

    ranks.clear();
    num_edges = 0;
    std::vector<Edge> old_edges;
    std::swap(edges, old_edges);

    for (auto x = 0; x < capacity; ++x) {
        neighbor[x].clear();
        nb_index[x].clear();
    }

    for (auto p : old_edges) {
        add_edge(srank[p[0]], srank[p[1]]);
    }

    verify("after sort");
}


void dyngraph::clear()
{
    int x;
    while (!node.empty()) {
        x = node.back();
				node.pop_back();
        neighbor[x].clear();
        nb_index[x].clear();
    }
    ranks.clear();
}

void dyngraph::declare_node(const int x)
{
    node.add(x);
    nodeset.add(x);

		// if(x >= capacity) {
		// 		capacity = x+1;
		// 	  		neighbor.resize(capacity);
		// 	  		nb_index.resize(capacity);
		// }
}

void dyngraph::add_node(const int x)
{
#ifdef _VERIFY_MCGRAPH
    verify("before add node");
#endif
		
    declare_node(x);

    int i = degree(x), y, e, pos;
    num_edges += i;

    while (i--) {
        y = neighbor[x][i];
        e = nb_index[x][i];

        // we add x at the back of y's neighbor list
        pos = (e & 1);

        // we change the position of x in y's neighbor list
        ranks[e / 2][1 - pos] = neighbor[y].size();

        // points to the edge from y's perspective
        nb_index[y].push_back(e ^ 1);

        // add x in y's neighbors
        neighbor[y].push_back(x);

    }

#ifdef _VERIFY_MCGRAPH
    verify("after add node");
#endif
}

void dyngraph::rem_node(const int x)
{	
		// std::cout << "remove " << x << std::endl;
#ifdef _VERIFY_MCGRAPH
    verify("before rem node");
#endif
	
    node.remove(x);
    nodeset.remove(x);
    auto i = degree(x);
    int y, z, rx, ex, posx, ey, posy;
    num_edges -= i;

    while (i--) {

        y = neighbor[x][i];

        ex = nb_index[x][i];
        posx = ex & 1;

        // store the position of x in y's neighborhood
        rx = ranks[ex / 2][1 - posx];

				assert(edges[ex / 2][posx] == x);
				assert(edges[ex / 2][1- posx] == y);
				assert(neighbor[y][rx] == x);
				assert(ranks[ex / 2][posx] == i);
				

				// if(y == 73 or y == 74)
				// std::cout << " - from N(" << y << "): " << edges[ex / 2] << " where it is at rank " << rx << "\n";


        // replace x by z
        z = neighbor[y].back();
        neighbor[y].pop_back();
				
				// if(y == 73 or y == 74)
				// std::cout << "replace edge " << edges[ex / 2] << " by edge " << edges[ nb_index[y].back() / 2 ] << std::endl;
				
				
				
				if(z != x) {
						neighbor[y][rx] = z;

		        // set the new position of z in y's neighborhood
		        ey = nb_index[y].back();
		        posy = (ey & 1);
		        ranks[ey / 2][posy] = rx;

		        //
		        nb_index[y][rx] = ey;
				} 
				
				nb_index[y].pop_back();
    }

#ifdef _VERIFY_MCGRAPH
    verify("after rem node");
#endif
}

bool dyngraph::has_node(int x) const { return nodeset.fast_contain(x); }

int dyngraph::add_edge(const int x, const int y)
{
    assert(node.contain(x));
    assert(node.contain(y));

    nb_index[x].push_back(2 * ranks.size());
    nb_index[y].push_back(2 * ranks.size() + 1);

    Edge r(neighbor[x].size(), neighbor[y].size());
    ranks.push_back(r);

    Edge e(x, y);
    edges.push_back(e);

    neighbor[x].push_back(y);
    neighbor[y].push_back(x);

    ++num_edges;

#ifdef _VERIFY_MCGRAPH
    verify("after add edge");
#endif

    return ranks.size() - 1;
}

void dyngraph::rem_edge(const int i)
{
    Edge edge = edges[i];
    rem_edge(edge[0], edge[1], i);
}
void dyngraph::rem_edge(const int x, const int y, const int e)
{

#ifdef _VERIFY_MCGRAPH
    assert(node.contain(x));
    assert(node.contain(y));
#endif

    int ry = ranks[e][0];
    int rx = ranks[e][1];

    int sx = neighbor[y].back();
    neighbor[y].pop_back();

    int ey = nb_index[y].back();
    nb_index[y].pop_back();
    if (x != sx) {
        int posy = (ey & 1);
        ranks[ey / 2][posy] = rx;

        neighbor[y][rx] = sx;
        nb_index[y][rx] = ey;
    }

    int sy = neighbor[x].back();
    neighbor[x].pop_back();

    int ex = nb_index[x].back();
    nb_index[x].pop_back();
    if (y != sy) {
        int posx = (ex & 1);

        ranks[ex / 2][posx] = ry;

        neighbor[x][ry] = sy;
        nb_index[x][ry] = ex;
    }

    Edge ez = edges.back();
    edges.pop_back();
    Edge rz = ranks.back();
    ranks.pop_back();

    if (static_cast<unsigned>(e) != edges.size()) {

        int s;
        for (int i = 0; i < 2; ++i) {
            s = (nb_index[ez[i]][rz[i]] & 1);
            nb_index[ez[i]][rz[i]] = 2 * e + s;
        }

        edges[e] = ez;
        ranks[e] = rz;
    }

    --num_edges;

#ifdef _VERIFY_MCGRAPH
    verify("after rem edge");
#endif
}

void dyngraph::maximal_matching(
    std::vector<int>& matching, int& nmatch, std::vector<int>& ranklist)
{
    ranklist.clear();
    for (int i = 0; i < size(); i++)
        ranklist.push_back(i);
		
		// std::random_device rd;
		// std::mt19937 g(rd());
		//     std::shuffle(ranklist.begin(), ranklist.end());

    int u, v;
    matching.assign(capacity, -1);
    nmatch = 0;
    for (auto i = 0; i < size(); i++) {
        u = node[ranklist[i]];
        if (matching[u] == -1) {
            for (auto j = 0; j < neighbor[u].size(); j++) {
                v = neighbor[u][j];
                if (matching[v] == -1) {
                    matching[u] = v;
                    matching[v] = u;
                    nmatch++;
                    break;
                }
            }
        }
    }
}

std::ostream& dyngraph::display(std::ostream& os) const
{

    for (auto i = 0; i < size(); ++i) {
        int x = node[i];

        os << x << ": ";
        vecdisplay(neighbor[x], os);
				os << " (" << degree(x) << ")" << std::endl;
				
				// os << "   [";
				// for(auto e : nb_index[x])
				// 		os << " " << (e/2);
				// os << " ]\n";
    }
		// for(auto e : edges)
		// 		std::cout << e << std::endl;
			
    return os;
}

void dyngraph::print_dimacs(std::ostream& os) const
{
    os << "c Generated by minicsp" << std::endl
       << "p edge " << size() << " " << edges.size() << std::endl;

    // not true during search
    assert(static_cast<int>(edges.size()) == num_edges);

    for (auto edge : edges) {
        os << "e " << (edge[0] + 1) << " " << (edge[1] + 1) << std::endl;
    }
}

void dyngraph::verify(const char* msg)
{
    for (unsigned i = 0; i < edges.size(); ++i) {
        Edge e = edges[i];
        int x = e[0];
        int y = e[1];

        Edge r = ranks[i];
        int ry = r[0];
        int rx = r[1];


				// std::cout << "rank of " << y << " N(" << x << ") = " << ry << std::endl;
				// for()
				

        if (node.contain(y) && neighbor[x][ry] != y) {
            std::cout << msg << " " << i << "-th edge " << e << " points to "
                      << r << ", however the " << ry << "-th element of ";
            vecdisplay(neighbor[x], std::cout);
            std::cout << " is not " << y << std::endl;
            assert(0);
        }
        if (node.contain(y) && static_cast<unsigned>(nb_index[x][ry]) != 2 * i) {
            std::cout << msg << " " << i << "-th edge " << e << " points to "
                      << r << ", however the " << ry << "-th element of ";
            vecdisplay(nb_index[x], std::cout);
            std::cout << " is not " << (2 * i) << std::endl;
            assert(0);
        }

        if (node.contain(x) && neighbor[y][rx] != x) {
            std::cout << msg << " " << i << "-th edge " << e << " points to "
                      << r << ", however the " << rx << "-th element of ";
            vecdisplay(neighbor[y], std::cout);
            std::cout << " is not " << x << std::endl;
            assert(0);
        }
        if (node.contain(x)
            && static_cast<unsigned>(nb_index[y][rx]) != 2 * i + 1) {
            std::cout << msg << " " << i << "-th edge " << e << " points to "
                      << r << ", however the " << ry << "-th element of ";
            vecdisplay(nb_index[y], std::cout);
            std::cout << " is not " << (2 * i + 1) << std::endl;
            assert(0);
        }
    }

    int ecount = 0;
    for (int i = 0; i < size(); ++i) {
        int x = node[i];

        for (unsigned j = 0; j < neighbor[x].size(); ++j) {

            ++ecount;
            int y = neighbor[x][j];
            int ey = nb_index[x][j];

            Edge e = edges[ey / 2];

            int posx = (ey & 1);

            if(e[posx] != x) {
		            std::cout << msg << " " << j << "-th neighbor of " << x << " is " << y << " but " << x << " is not the " << posx << "-th node of edge " << e << std::endl;
		            assert(0);
            }
            if(e[1 - posx] != y) {
		            std::cout << msg << " " << j << "-th neighbor of " << x << " is " << y << " but " << y << " is not the " << 1-posx << "-th node of edge " << e << std::endl;
		            assert(0);
            }
						

            Edge r = ranks[ey / 2];

            if(static_cast<unsigned>(r[posx]) != j) {
	            std::cout << msg << " " << x << "'s rank in N(" << x << ") is " << r[posx] << ", should be " << j << std::endl;
							vecdisplay(nb_index[y], std::cout);
							std::cout << std::endl;
	            assert(0);
            }
        }
    }
    ecount /= 2;

    if (ecount != num_edges) {
        std::cout << " " << msg << ": num_edges = " << num_edges
                  << " real count = " << ecount << std::endl
            ;

        assert(0);
    }
}

// void dyngraph::brelaz_color(coloring& col) {
// 		if(node.empty()) return;
//
// 		// first order the nodes by degree
// 		col.degree.resize(node.size());
// 		for(auto v : node) {
// 				col.order.push_back(v);
// 				col.degree[v] = degree(v);
// 		}
//
// 		    std::sort(col.order.begin(), col.order.end(),
// 		        [&](const int x, const int y) { return (degree(x) > degree(y)); });
//
//
// 		for(auto v : col.order) {
// 			std::cout << v << " " << degree(v) << std::endl;
// 		}
//
// 		int c;
// 		int x = order[0]; // the next vertex to color
// 		int num_colors = 0; // the number of colors used so far
//
// 		col.satur.resize(node.size());
// 		col.first.push_back(0); // every node has saturation degree 0
//
// 		do {
// 				// forbidden_col.clear();
// 				// for(int i=neighbor[x].size(); i<col.degree[i]; ++i) {
// 				// 		forbidden_col.add(col.color[neighbor[x][i]]);
// 				// }
// 				// for(c=0; c<num_colors; ++c)
// 				// 		if(!forbidden_col.contain(c)) break;
// 				c = col.satur[x].get();
// 				col.color[x] = c;
// 				for(auto y : neighbor[x])
// 						if(col.satur[y].add(c)) {
// 							if(col.first.size() <= col.satur[y].size)
// 							// swap
//
// 						}
//
//
// 		} while( !node.empty() );
//
//
//
// }





std::ostream& gc::operator<<(std::ostream& os, const dyngraph& x)
{
    return x.display(os);
}

std::ostream& gc::operator<<(std::ostream& os, const dyngraph* x)
{
    return (x ? x->display(os) : os);
}
