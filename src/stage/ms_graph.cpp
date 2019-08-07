//============================================================================//
//====// Includes //==========================================================//

#include "ms_graph.hpp"


#include <chrono>
#include <thread>

//============================================================================//
//====// Namespace //=========================================================//
using namespace gc;

//============================================================================//
//====// Namespace //=========================================================//

//>> The threshold beyond which we consider that there's no improvement.
//>> (ie: [-RC_EPS < lambda < RC_EPS] <=> [lambda = 0])
#define RC_EPS 1e-6

//>> DEFINE
#define NOTHING_FOUND -1
#define HAS_BEEN_FOUND >=0

//============================================================================//
//====// Public Methods //====================================================//

int ms_graph::ms_select_node() const { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	if (this->_ms_sn_mode == SN_MODE_MIN_SWAP) {
		return this->_ms_select_min_swap();
	} else if (this->_ms_sn_mode == SN_MODE_SORT) {
		return this->_ms_select_by_sorting();
	} else if (this->_ms_sn_mode == SN_MODE_NEAR) {
		return this->_ms_select_nearest_improvement();
	} else {
		return this->_ms_select_min_swap();
	}
}

/// //////////////////////////////////////////////////////////////////////// ///

std::vector<int> ms_graph::ms_find_neighbors(const int u) const { //<<<<<<<<<<// 
    return this->matrix[u];
}

/// //////////////////////////////////////////////////////////////////////// ///

void ms_graph::ms_remove_node(const int u) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//    
	this->nodes.remove(u);
}

/// //////////////////////////////////////////////////////////////////////// ///

void ms_graph::ms_remove_nodes(const std::vector<int> su) { //<<<<<<<<<<<<<<<<//
	for (int u : su) {
		this->ms_remove_node(u);
	}
}

/// //////////////////////////////////////////////////////////////////////// ///

int ms_graph::ms_forward() { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	//>> Init
	int u = NOTHING_FOUND;
	
	//>> Select a node
	u = this->ms_select_node();
	
	//>> Update the graph if a node is found
	if (u HAS_BEEN_FOUND) {
		//>> Create a node
		ms_pChoice choice = std::make_shared<ms_Choice>();
		choice->id            = u;
		choice->is_selected   = true;
		choice->node_restore_point = this->nodes.size();
		
		//>> Find its neighbors
		std::vector<int> neighbors = this->ms_find_neighbors(u);
		
		//>> Remove u from nodes
		this->ms_remove_node(u);		
		
		//>> Update the neighbors restore point
		choice->neighbors_restore_point = this->nodes.size();
		
		//>> Remove u's neighbors from nodes
		this->ms_remove_nodes(neighbors);		
		
		//>> Add it to the current set
		this->_ms_current_set.push_back(u);

		//>> Add it to the current branch
		this->_ms_choices.push_back(choice);
	}
	
	//>> Return the node or NOTHING_FOUND
	return u;	  
}

/// //////////////////////////////////////////////////////////////////////// ///

void ms_graph::ms_deep_forward() { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	int node = this->ms_forward();
	while (node HAS_BEEN_FOUND) {
		//if ( this->ms_score()-this->_ms_obj >= -RC_EPS ) break;
		/* /!\ |^|-------------------------------------------|^| /!\ */
		/* /!\ |^| It looks like it slows the process a lot, |^| /!\ */
		/* /!\ |^| I don't know why.                         |^| /!\ */
		/* /!\ |^| I have not investigated it because search |^| /!\ */
		/* /!\ |^| is quite fast compared to the rest of the |^| /!\ */
		/* /!\ |^| column generation methods, so far         |^| /!\ */
		/* /!\ |^|-------------------------------------------|^| /!\ */

		node = this->ms_forward();
	}
}

/// //////////////////////////////////////////////////////////////////////// ///

int ms_graph::ms_backward() { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	//>> Init
	int u = NOTHING_FOUND;
	
	//>> Restore if not empty
	if (not(this->_ms_choices.empty())) {
		//>> Restore to the previous node restore point
		this->nodes.restore(this->_ms_choices.back()->node_restore_point);
		
		//>> Remove the node of the set
		this->_ms_current_set.pop_back();
		
		//>> Remove the choice from the tree
		this->_ms_choices.pop_back();
	}
	
	//>> Return the node removed from the branch or NOTHING_FOUND
	return u;	
}

/// //////////////////////////////////////////////////////////////////////// ///

int ms_graph::ms_backtrack() { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	//>> Init
	int u = NOTHING_FOUND;
	bool success = false;
	
	//>> Nothing to do if empty
	if (this->_ms_choices.empty()) {
		return NOTHING_FOUND;
	}
	
	//>> Backtrack untill successful or the tree is empty
	while (not(success)) {
		//>> Try to switch off the current choice
		if (this->_ms_choices.back()->is_selected) {
			this->_ms_choices.back()->is_selected = false; // switch off!
			u = this->_ms_choices.back()->id;
			success = true;
		} else {
			u = this->ms_backward();
			if (not(u HAS_BEEN_FOUND)) { // tree is empty
				success = true;
			}
		}
	}
	
	//>> Return the switched node or NOTHING_FOUND
	return u;
}

/// //////////////////////////////////////////////////////////////////////// ///

void ms_graph::ms_restore_to_zero() {  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	//>> Nothing to do if already empty
	if (this->_ms_choices.empty()) {
		return;
	}
	
	//>> Restore nodes to its origin point
	this->nodes.restore(this->_ms_choices.front()->node_restore_point);
	
	//>> Empty the set and the tree
	this->_ms_current_set.clear();
	this->_ms_choices.clear();
}

/// //////////////////////////////////////////////////////////////////////// ///

std::vector<int> ms_graph::ms_find_set(const std::vector<float> price) { //<<<//
	//>> Init
	int UB;
	int ub;
	int node;
	std::vector<int> out, fout;

	//>> Set prices
	this->_ms_prices = price;

	//>> First forward
	this->ms_deep_forward();

	//>> Retrieve score
	ub = this->ms_score();
	UB = ub;

	while(UB > this->_ms_obj+RC_EPS) {
		//>> Switch to the nearest branch
		node = this->ms_backtrack();
		if (node HAS_BEEN_FOUND) {
			//>> Go to the leaf
			this->ms_deep_forward();
			//>> Retrieve score
			ub = this->ms_score();
			UB = std::max(ub, UB);

		} else {
			break;
		}
	}

	
	//>> Save result 
	out = this->_ms_current_set;
	
	//>> Transform the current "out" to fit the appropriate format
	for (int n=0 ; n<int(this->matrix.size()) ; n++) {
		//>> Find the representative node
		int p = this->parent[n];
		while ( p != this->parent[p] ) {
			p = this->parent[p];
		}
		//>> In or out ?
		if (std::find(out.begin(), out.end(), p) != out.end()) {
			fout.push_back(1);
		} else {
			fout.push_back(0); 
		}
	}

	//>> Restore the ->nodes
	this->ms_restore_to_zero();

	//>> Output
	return fout;
	
} 

/// //////////////////////////////////////////////////////////////////////// ///

float ms_graph::ms_score() { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<// 
	//>> Init
	float score = 0;

	//>> Known score for the selected nodes
	for (int i : this->_ms_current_set) {
		score+=this->_ms_prices[i];
	}

	//>> Supposed score for the rest of the nodes
	score+=this->ms_cplt();

	//>> Output
	return score;

}

/// //////////////////////////////////////////////////////////////////////// ///

float ms_graph::ms_cplt() { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	//>> Init
	float cplt=0;
	std::vector<int> buffer;
	std::vector<int> C, K; //>> "Candidat, 'K'lique"
	std::vector<float> current_price = this->_ms_prices;	
	int u  = NOTHING_FOUND;
	int v  = NOTHING_FOUND;
	int mu = NOTHING_FOUND;
	int mv = NOTHING_FOUND;

	//>> Save this->nodes state
	size_t restore_point = this->nodes.size();

	//>> Do the thing
	while (not(this->nodes.empty())) {
		//>> Cleaning
		K.clear();
		C.clear();

		//>> Select a node		
		v = this->nodes[0];
		K.push_back(v);

		//>> Select neighbors
		buffer = this->ms_find_neighbors(v);

		//>> Reduce to valid candidates
		for (int b : buffer) {
			if (this->nodes.contain(b)) {
				C.push_back(b);
			}
		}
	
		//>> Create a 'K'lique until there's no more Candidat
		while (not(C.empty())) {
			//>> Select a node			
			u = C[0];

			//>> Find its neighbors in C
			C = inter(C, this->ms_find_neighbors(u));

			//>> Add u to K
			K.push_back(u);

			//>> Update mv and mu
			if (int(K.size()) == 2) {
				//>> Init mu and mv
				if (-current_price[u] < -current_price[v]) {
					mu = u;
					mv = v;
				} else {
					mu = v;
					mv = u;
				}	
			} else if (-current_price[u] > -current_price[mv]) {
				mu = mv;
				mv = u;
			} else if (-current_price[u] > -current_price[mu]) {
				mu = u;
			}		
		}

		//>> Use K to update the cplt and current_price
		if (int(K.size()) == 1) {
			//>> Increase cplt
			cplt += current_price[K[0]];
		
			//>> Remove K[0]
			this->ms_remove_node(K[0]);

		} else {
			//>> Increase cplt
			cplt += current_price[mu];

			//>> Remove all element of K but v from this->nodes
			for (int k : K) {
				if (k != mv) {
					this->ms_remove_node(k);
				}
			}
			
			//>> Update prices of mv
			current_price[mv] -= current_price[mu];
		}		
		
	}

	//>> Restore nodes
	this->nodes.restore(restore_point);
	
	//>> Output
	return cplt;
}

/// //////////////////////////////////////////////////////////////////////// ///

void ms_graph::_ms_print_current_set() { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	std::cout << "[ ";
	for (int a : this->_ms_current_set) {
		std::cout << a << " ";
	}
	std::cout << "]" << std::endl;
}

/// //////////////////////////////////////////////////////////////////////// ///

void ms_graph::ms_set_node_sorter(const NodeSorter ns) { //<<<<<<<<<<<<<<<<<<<// 
	this->_ms_ns = ns;
}  

/// //////////////////////////////////////////////////////////////////////// ///

void ms_graph::ms_set_sn_mode(const SN_MODE sn_mode) { //<<<<<<<<<<<<<<<<<<<<<//
	this->_ms_sn_mode = sn_mode;
}

/// //////////////////////////////////////////////////////////////////////// ///

void ms_graph::ms_set_obj(const int obj) { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	this->_ms_obj = obj;
}

//============================================================================//
//====// Private Methods //===================================================//

int ms_graph::_ms_select_min_swap() const {
	if (this->nodes.empty()) {
		return NOTHING_FOUND;
	} else {
		return this->nodes[int(this->nodes.size())-1];
	}
}

/// //////////////////////////////////////////////////////////////////////// ///

int ms_graph::_ms_select_by_sorting() const { //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
	//>> Construct the nodes
	std::vector<ms_Node> vec;
	for (int i=0 ; i<int(this->nodes.size()) ; i++) {
		ms_Node n;
		n.id = this->nodes[i];
		n.weight = -this->_ms_prices[n.id];
		n.degree = int(this->ms_find_neighbors(n.id).size());
		vec.push_back(n);		
	}

	//>> Call the sorter
	std::sort(vec.begin(), vec.end(), this->_ms_ns);

	//>> Return the max
	if (vec.empty()) {
		return NOTHING_FOUND;
	} else { 
		return vec.back().id;
	}
}

/// //////////////////////////////////////////////////////////////////////// ///

int ms_graph::_ms_select_nearest_improvement() const { //<<<<<<<<<<<<<<<<<<<<<//
	//>> pick the first node with a negative price.	
	for (int i=0 ; i<int(this->nodes.size()) ; i++) {
		if (this->_ms_prices[this->nodes[i]] < 0) {
			return this->nodes[i]; 
		}
	}
	//>> If no node seems to be good, let's something else pick a node.
	return this->_ms_select_min_swap();
}

//============================================================================//
//====// Operators //=========================================================//

std::vector<int> gc::inter(const std::vector<int> &a, const std::vector<int> &b) { //<<//
	std::vector<int> out;
	for (int alpha : a) {
		if (std::find(b.begin(), b.end(), alpha) != b.end()) {
			out.push_back(alpha);
		}
	}
	return out;
} 










