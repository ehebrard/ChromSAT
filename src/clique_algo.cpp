
#include "clique_algo.hpp"
#include "cliquesampler.hpp"

// #define TRACE

namespace gc 
{
	class EndSearch : public std::exception {
	};


	int clique_algo::clique_size() const { return current.num_front(); }
	
	
	void clique_algo::initialise_lower_bound() {
		
		
    gc::clique_sampler<int> cs;
    cs.set_seed(12345);

    cs.set_domain(begin(df.order), end(df.order), G.capacity(), true);

    lb = (G.num_edge() > 0 ? 2 : 1);
    size_t width{1};

    gc::no_weight<int> no_weights{1};

    while (width <= 32) {

        auto nlb{cs.find_clique(
            G, lb, end(df.order), end(df.order), 128, width, no_weights)};

        if (nlb > lb) {
            lb = nlb;
        } else {
            width *= 2;
        }
    }	
	}
	
	
	
	int clique_algo::bound() {
		

		
		
		brelaz.clear();
		
		auto ncolor{brelaz.dsatur(G, ub - clique_size() + 1)};
	
		auto local_ub{std::min(ub, clique_size() + ncolor)};
		
	
		if(local_ub > lb) {
			brelaz.close(G);
	
			coloring = brelaz.get_coloring();
			
			// std::cout << "Required degeneracy: " << (lb - clique_size()) << std::endl;
			
			auto chrom_deg{brelaz.degeneracy(G, false, lb - clique_size())};
	
		
			local_ub = std::min(ub, clique_size() + chrom_deg + 1);
			
		
			// for(auto v : brelaz.order) {
			// 	std::cout << std::setw(3) << v << ": " << std::setw(3)
			// 		<< G.degree(v) << " / " << std::setw(3)
			// 			<< std::left << brelaz.neighbor_colors[v].count() << std::right
			// 				<< " => " << std::setw(2) << coloring[v]
			// 					<< (v == *(brelaz.core) ? "**" : "")
			// 						<< std::endl;
			// }
			// //
			// // exit(1);
			
		}
		
#ifdef TRACE		
		for(auto i{0}; i<dtrail.size(); ++i)
			std::cout << " "; 
		std::cout << "bounds: [" << lb << "," << local_ub << "]\n";
#endif
				
		return local_ub;
	}
	
	void clique_algo::get_branches() {
		
		commit();
		
#ifdef TRACE
		for(auto i{0}; i<dtrail.size(); ++i)
			std::cout << " "; 
		std::cout << "branches:";	
#endif
				
		for(auto v{brelaz.core}; v!=brelaz.order.rend(); ++v) {
			if(coloring[*v] + clique_size() >= lb) {
				decision_stack.push_back(*v);

// #ifdef TRACE
// 				std::cout << " x" << *v ;
// #endif
				
				assert(current.contain(*v));
				assert(G.contain(*v));
								
			}
		}
		
#ifdef TRACE
		// std::cout << "\n";
		
		int z{3};
		for(int i{0}; i<decision_stack.size(); ++i) {
			if(z < dtrail.size() and i == dtrail[z]) {
				std::cout << " [" << (z/2) << "]";
				z += 2;
			}
			std::cout << " " << decision_stack[i];
		}
		std::cout << "\n";
#endif		
		
	}
	
	bool clique_algo::dead_end() {
		
		return dtrail.back() == decision_stack.size();
		
	}
	
bool clique_algo::root_node() {
	
	// std::cout << "root? " << dtrail.size() << std::endl;
	
	return dtrail.size() <= 2;
}
	
	void clique_algo::branch() {
		
		add_to_clique(decision_stack.back());
		decision_stack.pop_back();
		
	}
		
	
	

void clique_algo::search() {
	
	
	try {
	while(true) {
		
		// for(int v{0}; v<G.capacity(); ++v) {
		// 	assert(current.contain(v) >= G.contain(v));
		// }
		
		if(G.empty())
			std::cout << "solution?\n";
		
				
		auto local{bound()};	
		
		if(root_node()) {
			ub = local;
			std::cout << "bounds: [" << lb << "," << ub << "]\n";
		}
		
		if(local > lb) {
			get_branches();
		} else {
			++num_fails;
			undo();
		}
		while(dead_end()) { 
			backtrack();
		} 
		branch();	
	}
} catch(std::exception& e) {
	ub = lb;
	std::cout << "clique max = " << ub << " (" << num_fails << ")\n";
}
	
	
}





void clique_algo::backtrack() {
	
	// decision_stack
	dtrail.pop_back();
	
	current.restore_back(dtrail.back());
	dtrail.pop_back();
	
	if(root_node())
		throw EndSearch();
	else
		undo();

#ifdef TRACE	
	for(auto i{0}; i<dtrail.size(); ++i)
		std::cout << " "; 
	std::cout << "backtrack\n";
#endif
	
}

void clique_algo::commit() {
	
	size_t st;
	current.save_back(st);
	dtrail.push_back(st);
	dtrail.push_back(decision_stack.size());
	
#ifdef TRACE
	for(auto i{0}; i<dtrail.size(); ++i)
		std::cout << " "; 
	std::cout << "commit\n";
#endif
	
}

void clique_algo::undo() {
	
	
	// for(auto i{0}; i<current.capacity(); ++i)
	// 	std::cout << " " << current[i-current.num_front()] ;
	// std::cout << std::endl;
	//
	// for(auto i{current.begin_front()}; i!=current.end_front(); ++i)
	// 	std::cout << " " << *i;
	// std::cout << " |";
	// for(auto i{current.begin()}; i!=current.end(); ++i)
	// 	std::cout << " " << *i;
	// std::cout << " |";
	// for(auto i{current.begin_back()}; i!=current.end_back(); ++i)
	// 	std::cout << " " << *i;
	// std::cout << "\n";
	//
	//
	// for(auto i{current.rbegin_back()}; i!=current.rend_back(); ++i)
	// 	std::cout << " " << *i;
	// std::cout << " |";
	// for(auto i{current.rbegin()}; i!=current.rend(); ++i)
	// 	std::cout << " " << *i;
	// std::cout << " |";
	// for(auto i{current.rbegin_front()}; i!=current.rend_front(); ++i)
	// 	std::cout << " " << *i;
	// std::cout << "\n";
	
	
	auto v{*(current.rbegin_front())}; // vertex added to the clique
	current.add(v);
	current.remove_back(v);
	
#ifdef TRACE
	for(auto i{0}; i<dtrail.size(); ++i)
		std::cout << " "; 
	std::cout << "undo " << v << " @" << (dtrail.size() / 2 - 1) << "\n" ;	
#endif
	
	G.undo();

}

void clique_algo::add_to_clique(const int v) {
	
#ifdef TRACE
	for(auto i{0}; i<dtrail.size(); ++i)
		std::cout << " "; 
	std::cout << "try " << v << " @" << (dtrail.size() / 2 - 1) << " (";
	
	for(auto v{current.begin_front()}; v!=current.end_front(); ++v)
		std::cout << " " << *v;
	std::cout << " /";
	for(auto v{current.begin_back()}; v!=current.end_back(); ++v)
		std::cout << " " << *v;
	
	std::cout << ")\n" ;	
#endif 

	assert(current.contain(v));
	
	assert(G.contain(v));
	
	current.remove_front(v);

	
	G.add_to_clique(v, current);

}


} // namespace gc

