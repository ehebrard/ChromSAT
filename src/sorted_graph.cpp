
#include <iostream>

#include "sorted_graph.hpp"




namespace gc
{
	
	sorted_graph::sorted_graph(const int n) {
		nodes.reserve(n);
		matrix.resize(n);
		nodeset.resize(n, false);
		for(auto i{0}; i<n; ++i)
			add_vertex(i);
		// dtrail.push_back(current.size());
	}
	
	void sorted_graph::reserve(const int n) {
		nodes.reserve(n);
		nodeset.resize(n);
	}
	
	void sorted_graph::add_vertex(const int v) {
		if(not contain(v)) {
			nodeset.set(v);
			nodes.push_back(v);
		}
	}
	void sorted_graph::add_edge(const int u, const int v) {
		_num_edge += 2;
		matrix[u].push_back(v);
		matrix[v].push_back(u);
	}
	void sorted_graph::canonize() {
		std::sort(nodes.begin(), nodes.end());
		for(auto v : nodes) 
				std::sort(matrix[v].begin(), matrix[v].end());
	}
	

		bool sorted_graph::contain(const int v) const { return nodeset[v]; }
		int sorted_graph::degree(const int v) const { return matrix[v].size(); }
		int sorted_graph::num_edge() const { return _num_edge / 2; }
		int sorted_graph::size() const { return nodes.size(); }		
		int sorted_graph::capacity() const { return nodes.capacity(); }	
		bool sorted_graph::empty() const { return size() == 0; }	
		
		// std::vector<int>::iterator sorted_graph::begin_neighbor(const int v) { return neighbor[v].begin(); }
		// std::vector<int>::iterator sorted_graph::end_neighbor(const int v) { return neighbor[v].begin() + num_neigbor[v]; }
		
		
		// void sorted_graph::remove(const int v) {
		//
		//
		// }
		
		
		

		
		void sorted_graph::intersect_neighbors(const int v) {

			inbuffer.clear();
			outbuffer.clear();
			
			for(auto u : matrix[v]) {
				if(nodeset[u])
					inbuffer.push_back(u);
				else
					outbuffer.push_back(u);
			}
				
			auto n{matrix[v].begin()};
			auto i{inbuffer.begin()};
			while(i!=inbuffer.end())
				*n++ = *i++;
			auto o{outbuffer.begin()};
			while(o!=outbuffer.end())
				*n++ = *o++;
			
			trail.push_back(matrix[v].size());
			
			_num_edge -= matrix[v].size();
			_num_edge += inbuffer.size();
			
			matrix[v].set_size(inbuffer.size());
			
			trail.push_back(v);
		}
		
		

		
		
		void sorted_graph::undo() {
			// verify("before undo");
			
			
			// auto v{trail.back()}; // vertex added to the clique
			// current.add(v);
			// current.remove_back(v);
			
			
			
			// for(auto i{current.begin_front()}; i!=current.end_front(); ++i)
			// 	std::cout << " " << *i;
			// std::cout << " tried before:";
			// for(auto i{current.begin_back()}; i!=current.end_back(); ++i)
			// 	std::cout << " " << *i;
			// std::cout << std::endl;
			
			auto v{trail.back()}; // vertex added to the clique

#ifdef TRACE			
			std::cout << "undo (" << v << ")\n";
#endif
			
			trail.pop_back();
			
			auto s{trail.back()}; // old |V|
			trail.pop_back();
			
			auto j{nodes.begin()+s};
			for(auto i{nodes.end()}; i!=j; ++i) {
				nodeset.set(*i);
				_num_edge += matrix[*i].size();
			}
			
			inbuffer.clear();

			nodes.undo(s, inbuffer);
			
			while(true) {
				v = trail.back(); // vertex 
				trail.pop_back();
				
				if(v < 0)
					break;
			
				s = trail.back(); // old |N(v)|
				trail.pop_back();
			
				_num_edge += s;
				_num_edge -= matrix[v].size();
			
				inbuffer.clear();
				matrix[v].undo(s, inbuffer);
			}
			
			VERIFY("after undo");	
		}
		
		
		std::ostream& sorted_graph::describe(std::ostream& os, const int verbosity) const {
	    switch (verbosity) {
	    case 0:
	        os << "n=" << size() << "; m=" << num_edge() << std::endl;
	        break;

	    case 1:
	        os << "V=" ;
					for(auto v : nodes) 
						os << " " << v;
					os << "; m=" << num_edge();
	        break;

	    case 2:
	        for(auto v : nodes) {
	            os << v << ":" << matrix[v] << std::endl;
	        }
	        break;

	    case 3:
      for(auto v : nodes) {
          os << v << ":";
					matrix[v].describe(os, 3);
          os << std::endl;
      }
			os << std::endl;
			for(auto v{nodes.end()}; v!=nodes.begin()+nodes.capacity(); ++v) {
				if(not nodeset[*v])
					os << "("<< *v << "):";
				else
					os << *v << ":";
				matrix[*v].describe(os, 3);
        os << std::endl;
			}
	        for (auto t : trail)
	            os << " " << t;
	        os << std::endl;
	    }

	    return os;
		}
		

void sorted_graph::verify(const char* msg) {
	
	// std::cout << "verify " << msg << "\nN(34)=";
	// for(auto v : matrix[34]) {
	// 	std::cout << " " << v;
	// }
	// std::cout << std::endl << nodes << std::endl << nodeset[0] << std::endl;
	//
	
	if(nodeset.count() != nodes.size()) {
		std::cout << msg << ": node count (" << nodeset.count() << "/" << nodes.size() << ")\n";
		
		for(auto i{0}; i<nodes.capacity(); ++i)
			if(nodeset[i])
				std::cout << " " << i;
		std::cout << std::endl << nodes << std::endl;
		
		
		
		exit(1); 
	} 
	
	auto last{-1};
	auto edge_count{0};
	auto vptr{nodes.begin()};
	for(auto v{0}; v<capacity(); ++v) {
		if(vptr!=nodes.end() and v == *vptr) {
			++vptr;
			edge_count += matrix[v].size();
		
			if( not nodeset[v] ) {
				std::cout << msg << ": " << v << " is in nodes but not in nodeset\n";
				exit(1); 
			}
			if( last >= v ) {
				std::cout << msg << ": nodes is not ordered (..." << last << ", " << v << "...\n";
				exit(1); 
			}
		
			for(auto u : matrix[v]) {
				if( not nodeset[u] ) {
					std::cout << msg << ": neighbor[" << v << "] contains " << u << " which is not in nodeset\n";
					exit(1); 
				}
			}
		
			last = v;
		} else {
			if( nodeset[v] ) {
				std::cout << msg << ": " << v << " is in nodeset but in nodes\n";
				exit(1); 
			}
		}
	}
	
	if( _num_edge != edge_count ) {
		std::cout << msg << ": wrong edge count (" << _num_edge << " / " << edge_count << ")\n";
		exit(1); 
	}
}


std::ostream& operator<<(std::ostream& os, const sorted_graph& g) {
	g.describe(os, 3);
	return os;
}

std::ostream& operator<<(std::ostream& os, const vertex_list& v) {
	v.describe(os, 2);
	return os;
}

} // namespace gc
