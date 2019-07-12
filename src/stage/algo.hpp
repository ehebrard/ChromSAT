#ifndef __BNP_ALGO_HPP
#define __BNP_ALGO_HPP

//============================================================================//
//====// Includes //==========================================================//

#include <limits>
#include <vector>

#include "bi_graph.hpp"
#include "cliquer.hpp"
#include "dsatur.hpp"
#include "intstack.hpp"
#include "statistics.hpp"

//============================================================================//
//====// Defines //===========================================================//

//============================================================================//
//====// Namespace : begin //=================================================//

namespace gc {

//============================================================================//
//====// Enumerations //======================================================//

//============================================================================//
//====// Structures //========================================================//

//============================================================================//
//====// Typedef //===========================================================//

//============================================================================//
//====// Primitives //========================================================//

template <class graph_struct> class coloring_algorithm;

//============================================================================//
//====// Classes //===========================================================//

//############################################################################//
template <class graph_struct> class heuristic { //############################//

///////-/// ATTRIBUTES ///-/////////////////////////////////////////////////////
	
	private:
	coloring_algorithm<graph_struct>& env;
	degeneracy_finder<graph_struct>   df;
	std::vector<int>                  degeneracy;
	std::vector<int>                  dg_rank;

///////-/// CONSTRUCTORS ///-///////////////////////////////////////////////////
	
	public:
	heuristic(coloring_algorithm<graph_struct>& env);

///////-/// METHODS ///-////////////////////////////////////////////////////////
	
	public:
	int degeneracy_coloring(std::vector<int>& coloring);
	int dsatur_coloring(const int ub, std::vector<int>& coloring);
};

//############################################################################//
template <class graph_struct> class lower_bound { //##########################//

///////-/// ATTRIBUTES ///-/////////////////////////////////////////////////////
	
	private :
	coloring_algorithm<graph_struct>& env;
	std::vector<int>                  largest_cliques;
	gc::cliquer<graph_struct>         cq;
	gc::bi_graph                      B;

	public:
	int clique_sz;
	gc::bitset max_clique;

///////-/// CONSTRUCTORS ///-///////////////////////////////////////////////////

	public:
	lower_bound(coloring_algorithm<graph_struct>& env);

///////-/// METHODS ///-////////////////////////////////////////////////////////

	public:
	template <class random_it> int clique(random_it beg, random_it end);
	void                           get_largest_clique();
	int                            matching();
};

//############################################################################//
template <class graph_struct> class selector { //#############################//

///////-/// ATTRIBUTES ///-/////////////////////////////////////////////////////

	private:
	coloring_algorithm<graph_struct>& env;

///////-/// CONSTRUCTORS ///-///////////////////////////////////////////////////

	public:
	selector(coloring_algorithm<graph_struct>& env);

///////-/// METHODS ///-////////////////////////////////////////////////////////

	public:
	arc select();
	void make_choice();
};

//############################################################################//
template <class graph_struct> class coloring_algorithm { //###################//

///////-/// ATTRIBUTES ///-/////////////////////////////////////////////////////

	public:
	graph_struct&   G;
	gc::statistics& stats;
	gc::options&    options;

	int    UB;
	int    LB;
	int    depth;
	int    N;
	dsatur brelaz;

	heuristic<graph_struct>   h;
	lower_bound<graph_struct> l;
	selector<graph_struct>    s;

	std::vector<int> coloring;

///////-/// CONSTRUCTORS ///-///////////////////////////////////////////////////

	public:
	coloring_algorithm(graph_struct& g, statistics& stats, gc::options& options);

///////-/// METHODS ///-////////////////////////////////////////////////////////

	public:
	void update_lb(const int lb);
	void update_ub(const int ub);
	void backtrack();
	void find_coloring();

};

//============================================================================//
//====// Namespace : end //===================================================//

} // namespace gc

#endif // ALGORITHM
