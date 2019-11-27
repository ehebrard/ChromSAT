#ifndef __CG_SORTED_GRAPH_HH
#define __CG_SORTED_GRAPH_HH


#include <boost/dynamic_bitset.hpp>

#include <vector>


// #define _VERIFIED_

#ifdef _VERIFIED_
#define VERIFY(X) verify(X);
#else
#define VERIFY(X) 
#endif


namespace gc
{
	
	class sorted_graph;

	
		class vertex_list {

		public:

			size_t capacity() const { return vertices.size(); }
			int size() const { return num_vertex; }
			void reserve(const int n) { vertices.reserve(n); }
			void push_back(const int v) { vertices.push_back(v); ++num_vertex; }
			void set_size(const int n) { num_vertex = n; }
			
			
			void undo(const int n, std::vector<int>& buffer) {
				auto x{begin()};
				auto y{end()};
				auto ex{end()};
				auto ey{begin() + n};
				
				while(x != ex and y != ey) {
					if(*x < *y) {
						buffer.push_back(*x++);
						// if(++x == ex) break;
					} else {
						buffer.push_back(*y++);
						// if(++y == ey) break;
					}
				}
				while(x != ex)
					buffer.push_back(*x++);
				while(y != ey)
					buffer.push_back(*y++);
				x = begin();
				y = buffer.begin();
				
				while(x != ey)
					*x++ = *y++;
						
				set_size(n);
			}

			std::vector<int>::iterator begin() { return vertices.begin(); }
			std::vector<int>::iterator end() { return vertices.begin() + size(); }
			
			std::vector<int>::const_iterator begin() const { return vertices.begin(); }
			std::vector<int>::const_iterator end() const { return vertices.begin() + size(); }
			
			int& operator[](const int i) { return vertices[i]; }
			
			std::ostream& describe(std::ostream& os, const int verbosity) const {
				for(auto v : *this)
					os << " " << v;
				if(verbosity >= 3) {
					os << " |";
					for(auto v{end()}; v!=vertices.end(); ++v)
						os << " " << *v;
				}
				return os;
			}

		private:
			int num_vertex{0};
			std::vector<int> vertices;

		};
	

	class sorted_graph {
		
	public:
		sorted_graph() {}
		sorted_graph(const int n);
		void reserve(const int n);
		void add_vertex(const int v);
		void add_edge(const int u, const int v);
		void canonize();
			
		bool contain(const int v) const;
		int degree(const int v) const;
		int num_edge() const;
		int size() const;
		int capacity() const;
		bool empty() const;		
	
		template< class vertex_iterator >
		void remove(vertex_iterator b, vertex_iterator e);
		
		void remove(const int v);
		
		// template< class vertex_iterator >
		// void intersect(vertex_iterator b, vertex_iterator e);
		
		
		template<class container>
		void add_to_clique(const int v, container C);
		
		
		template<class container>
		void intersect(const int x, container C);
		
		std::ostream& describe(std::ostream& os, const int verbosity) const;
		
		vertex_list nodes;
		
		std::vector<vertex_list> matrix;
		
		void add_to_clique(const int v);
		
		void backtrack();
		void undo();
		void commit();
		
		
				
	private:
		boost::dynamic_bitset<> nodeset;
		
		int _num_edge{0};
				
		std::vector<int> inbuffer;
		std::vector<int> outbuffer;
		
		std::vector<int> trail;
		
		void intersect_neighbors(const int v);
		
		void verify(const char* msg);
		
	};
	
	
	// template< class vertex_iterator >
	// void sorted_graph::remove(vertex_iterator b, vertex_iterator e) {
	//
	// }
	
	// template< class vertex_iterator >


	std::ostream& operator<<(std::ostream& os, const sorted_graph& g);

	std::ostream& operator<<(std::ostream& os, const vertex_list& v);




	template<class container>
	void sorted_graph::add_to_clique(const int v, container C) {
		
		// VERIFY("before add clique");
		
		// assert(current.contain(v));
		// current.remove_front(v);
		
		
		trail.push_back(-1);
		
		


		// for(auto i{current.begin_back()}; i!=current.end_back(); ++i) {
		// 	nodeset.reset(*i);
		// }
		
		auto s{nodes.size()};
		// intersect(matrix[v].begin(), matrix[v].end());
		
		

		
		intersect(v, C);
		
		
		
		
		trail.push_back(s);
		trail.push_back(v);
		
		VERIFY("after add clique");
		// verify("after add clique");
	}
	
	
	template<class container>
	void sorted_graph::intersect(const int x, container C) { //, vertex_iterator b, vertex_iterator e) {
		outbuffer.clear();
	
		auto b{matrix[x].begin()};
		auto e{matrix[x].end()};
	
	
		// std::cout << " inter      " ;
		// nodes.describe(std::cout, 2);
		// std::cout << "\n with N(" << x << ") =";
		// for(auto x{b}; x!=e; ++x)
		// 	std::cout << " " << *x;
		// std::cout << std::endl;
	
		auto w{nodes.begin()}, u{nodes.begin()};
		for(auto v{b}; v != e; ++v)
		{
			// std::cout << *u << " <> " << *v << (C.contain(*u) ? "" : "**") << std::endl;
			
			assert(*u <= *v);
			
			while(*u < *v) 
				outbuffer.push_back(*u++);
		
			assert(*u == *v);
			
			if(not C.contain(*u)) {
				// std::cout << " --" << *u << std::endl;
				outbuffer.push_back(*u++);
			} else {
				*w++ = *u++;
			}
			// *w++ = *u++;
		}
	
		auto nsize{(w - nodes.begin())};
		
		for(auto v : outbuffer) {
			nodeset.reset(v);
			*w++ = v;
			_num_edge -= matrix[v].size();
		}
	
		while(u != nodes.end()) {
			_num_edge -= matrix[*u].size();
			nodeset.reset(*u);
			*w++ = *u++;
		}
	
		// std::cout << " ==> " << nodes << std::endl;
	
	
		// std::cout << "cur size = " << size() << " #removed = " << (size() - nsize) << " #remaining " << nsize << std::endl;
	
	
		nodes.set_size(nsize);
	
		for(auto v : nodes) 
			intersect_neighbors(v);
	
		// trail.push_back(nsize);
		// trail.push_back(v);
	}


} // namespace gc

#endif
