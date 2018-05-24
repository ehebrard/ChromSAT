
#include "sparse_graph.hpp"
#include "utils.hpp"
#include <minicsp/mtl/Heap.h>
#include <algorithm>
#include <iomanip>

namespace gc
{

/****************************************************************/
//							SPARSE								//
/****************************************************************/
sparse_clique_finder::sparse_clique_finder(const sparse_graph& g)
    : g(g)
    , num_cliques(1)
{
	last_clique.resize(g.capacity()); // OK
    cliques.resize(g.capacity()); // Could we divide it by 2 ?
    clique_sz.resize(g.capacity());
    candidates.resize(g.capacity());
	// Initialization of the candidates of the first clique (all vertices from g.nodes)
	candidates[0].reserve(g.capacity()); // Allocate the necessary space
	for(std::vector<int>::const_iterator it = g.nodes.begin(); it < g.nodes.end(); ++it ) {
		candidates[0].push_back(*it);
	}
}

void sparse_clique_finder::clear() { num_cliques = 0; }

// Declare a new clique (it should be associated with a vertex to put in)
void sparse_clique_finder::new_clique()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].clear();
	// Why set its size to 0 ?
    clique_sz[num_cliques] = 0;
	// All vertices candidates
	for(std::vector<int>::const_iterator it = g.nodes.begin(); it < g.nodes.end(); ++it ) {
		candidates[num_cliques].push_back(*it);
	}    
    ++num_cliques;
}

void sparse_clique_finder::new_color()
{
    assert(num_cliques < g.capacity());
    cliques[num_cliques].clear();
    clique_sz[num_cliques] = 0;
    candidates[num_cliques].clear();
    ++num_cliques;
}

void sparse_clique_finder::insert(int v, int clq)
{
    cliques[clq].push_back(v);
    ++clique_sz[clq];
	// HERE we want candidates[clq] to be the intersected with N(v)
	// Should be sorted first (OK by construction)
	vertices_vec::iterator it;
	it = set_intersection(candidates[clq].begin(), candidates[clq].end(), g.adjacency[v].begin(), g.adjacency[v].end() , candidates[clq].begin());
	candidates[clq].resize(it-candidates[clq].begin());
    last_clique[v] = clq;
}


/*
void sparse_clique_finder::insert_color(int v, int clq, bitset& diff)
{
    cliques[clq].fast_add(v);
    ++clique_sz[clq];
    diff.copy(g.matrix[v]);
    diff.setminus_with(candidates[clq]);
    candidates[clq].union_with(g.matrix[v]);
}
*/
// DISPLAY
void sparse_graph::describe(std::ostream& os) const
{
    os << "# vertices = " << capacity() << std::endl;
}

std::ostream& sparse_graph::display_adjacency(std::ostream& os) const
{
	os << "Displaying adjacency\n";
	int i(0);
    for (std::vector<vertices_vec>::const_iterator it = adjacency.begin(); it < adjacency.end(); ++it) {
        os << "Vertex " << i << "   Neighbors : ( ";
		i++;
		for (std::vector<int>::const_iterator itN = (*it).begin(); itN < (*it).end(); ++itN) {
			os << *itN << " ";
		}
	    os << ")";
		os << " degree =" << (*it).size() << "\n";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const sparse_graph& x)
{
    return x.display_adjacency(os);
}

void sparse_clique_finder::display()
{
	std::cout << " Maximal cliques " << std::endl;
	int i = 0;
	for(std::vector<vertices_vec>::const_iterator it = cliques.begin(); it != cliques.end(); ++it) {
		if (!(*it).empty()){
			std::cout << " Clique " << i << " { ";
			i++;
			for (std::vector<int>::const_iterator itN = (*it).begin(); itN < (*it).end(); ++itN) {
				std::cout << *itN << " ";
			}
			std::cout << " } size " << (*it).size() <<std::endl;
		}
	}
}

/****************************************************************/
//						BRONKERBOSCH							//
/****************************************************************/

BronKerbosch::BronKerbosch(const sparse_graph& g)
	:g(g)
	, num_cliques(0)
{
	maximal_cliques.resize(g.capacity()); // Look for an upper bound to #maximal_cliques
										  // to prevent memory issues
    clique_sz.reserve(g.capacity());
    //actual_clique.resize(g.capacity());
	candidates.reserve(g.capacity()); // 
	banned.reserve(g.capacity());	// Too large

	// Set the candidates as all nodes
	for(std::vector<int>::const_iterator it = g.nodes.begin(); it < g.nodes.end(); ++it ) {
		candidates.push_back(*it);
	}
}

void BronKerbosch::clear() { num_cliques = 0; }

// Add a clique to maximal_cliques
void BronKerbosch::add_max_clique(vertices_vec& clique)
{
	std::cout << "ADDING CLIQUE " << std::endl;
	assert(num_cliques < g.capacity());	
    maximal_cliques[num_cliques].clear();
	for(std::vector<int>::const_iterator it = clique.begin(); it < clique.end(); ++it ) {
		maximal_cliques[num_cliques].push_back(*it);
	}
	clique_sz[num_cliques] = maximal_cliques[num_cliques].size();    
    ++num_cliques;
}

vertices_vec BronKerbosch::intersect(vertices_vec const& v1,vertices_vec const& v2) // v1 & v2 are assumed to be sorted
{
	vertices_vec v_intersection;
	std::set_intersection(v1.begin(), v1.end(),
						  v2.begin(), v2.end(),
						  std::back_inserter(v_intersection));	
	return v_intersection;
}

vertices_vec BronKerbosch::unite_vector_element(vertices_vec& v,const int i) // v is assumed to be sorted
{
	// Increase size by one
	v.resize(v.size() + 1);
	
	// Case i > last of v (empty vector works too)
	if (v[v.size()-2] < i ){
		v[v.size()-1] = i;
		return v;
	}

	// Find where to place the new value
	auto upper = std::upper_bound(v.begin(), v.end(), i);

	// Shift the elements above i to the right
	for(std::vector<int>::reverse_iterator rit = v.rbegin(); &*rit != &*upper ; ++rit)	
		*rit = *(rit+1);

	// Insert new value
	*upper = i;
	
	return v;
}

// Display calls
void BronKerbosch::bronkerbosch_calls_display(vertices_vec clique, vertices_vec candidates, vertices_vec banned)
{
	std::cout << "BronKerbosch( {" ;
	for(std::vector<int>::const_iterator it1 = clique.begin(); it1 != clique.end(); ++it1){
		std::cout << " " << *it1;}
	std::cout << " }, {";
	for(std::vector<int>::const_iterator it2 = candidates.begin(); it2 != candidates.end(); ++it2){
		std::cout << " " << *it2;}
	std::cout << " }, {";
	for(std::vector<int>::const_iterator it3 = banned.begin(); it3 != banned.end(); ++it3){
		std::cout << " " << *it3;}
	std::cout << " })"<< std::endl;
}

void BronKerbosch::find_cliques_withoutPivot(vertices_vec clique, vertices_vec candidates, vertices_vec banned)
{
	bronkerbosch_calls_display(clique, candidates, banned);	

	if(candidates.empty() && banned.empty())
		add_max_clique(clique);
	
	for(std::vector<int>::const_reverse_iterator rit = candidates.crbegin(); rit != candidates.crend(); ++rit) {
		// Recursive call
		find_cliques_withoutPivot(unite_vector_element(clique, *rit), intersect(candidates, g.adjacency[*rit]), intersect(banned, g.adjacency[*rit]));
		// P := \ {v}
		candidates.pop_back(); // Vector -> more efficient to treat candidates backward

		// X := X U {v} // Use unite_vector_element
		banned.push_back(*rit); 
		std::sort (banned.begin(), banned.end());

/* Use print_container instead ?
		std::cout << "candidates left: ";
		for(std::vector<int>::const_iterator it2 = candidates.begin(); it2 != candidates.end(); ++it2){
			std::cout << " " << *it2;}
		std::cout << std::endl;

		std::cout << "banned : ";
		for(std::vector<int>::const_iterator it3 = banned.begin(); it3 != banned.end(); ++it3){
		std::cout << " " << *it3;}
		std::cout << std::endl;
*/
	}		
}

void BronKerbosch::order_by_degree()
{	
	by_degree.resize(g.capacity());
	degree.resize(g.capacity());
	for (auto v : g.nodes) {
        degree[v] = g.adjacency[v].size();
        //by_degree[degree[v]].push_back(v); FALSE
    }

	for(std::vector<int>::const_iterator it = by_degree.begin(); it != by_degree.end(); ++it){
	std::cout << " " << *it;}
	std::cout << std::endl;
}

// Find the vertice with max degree within two vectors of vertices
int pivot(vertices_vec v1, vertices_vec v2)
{
	
}

/*										//               R                      P                      X
void BronKerbosch::find_cliques_withPivot(vertices_vec clique, vertices_vec candidates, vertices_vec banned)
{
	bronkerbosch_calls_display(clique, candidates, banned);	

	if(candidates.empty() && banned.empty())
		add_max_clique(clique);

	// Choose a pivot among P U X with max_degree
	pivot(candidates, banned);
	restricted_candidates = intersect(candidates, g.adjacency[pivot]);
	// TREAT CANDIDATES BACKWARD
	for(std::vector<int>::const_reverse_iterator rit = restricted_candidates.crbegin(); rit != restricted_candidates.crend(); ++rit) {
		// RECURSIVE CALL
		find_cliques_withoutPivot(unite(clique, *rit), intersect(candidates, g.adjacency[*rit]), intersect(banned, g.adjacency[*rit]));
		// P := \ {v}
		candidates.pop_back(); // MUCH BETTER TO TREAT CANDIDATES BACKWARDS AND USE POP_BACK !

		// X := X U {v} // Need to be sorted -> Consuming 
		banned.push_back(*rit); 
		std::sort (banned.begin(), banned.end());

		std::cout << "candidates left: ";
		for(std::vector<int>::const_iterator it2 = candidates.begin(); it2 != candidates.end(); ++it2){
			std::cout << " " << *it2;}
		std::cout << std::endl;

		std::cout << "banned : ";
		for(std::vector<int>::const_iterator it3 = banned.begin(); it3 != banned.end(); ++it3){
		std::cout << " " << *it3;}
		std::cout << std::endl;
	}		
}
/*
		// R U {v}		
		actual_clique.push_back(*rit); 
		// P inter N(v)		
		vertices_vec::iterator itr;
		itr = set_intersection(candidates.begin(), candidates.end(), g.adjacency[*rit].begin(), g.adjacency[*rit].end() , candidates.begin());
		candidates.resize(itr-candidates.begin());
		// X inter N(v)
		itr = set_intersection(banned.begin(), banned.end(), g.adjacency[*rit].begin(), g.adjacency[*rit].end() , banned.begin());
		banned.resize(itr-banned.begin());

*/

void BronKerbosch::display()
{
	std::cout << " Maximal cliques " << std::endl;
	int i = 0;
	for(std::vector<vertices_vec>::const_iterator it = maximal_cliques.begin(); it != maximal_cliques.end(); ++it) {
		if (!(*it).empty()){
			std::cout << " Clique " << i << " { ";
			i++;
			for (std::vector<int>::const_iterator itN = (*it).begin(); itN < (*it).end(); ++itN) {
				std::cout << *itN << " ";
			}
			std::cout << " } size " << (*it).size() <<std::endl;
		}
	}
}



// degeneracy_finder::degeneracy_finder(const graph& g)
//     : g(g)
// {
//
//      degrees.resize(g.capacity());
//      iterators.resize(g.capacity());
//      ordered.initialise(0, g.capacity()-1, bitset::empt);
// }
//
//
//
// void degeneracy_finder::get_degeneracy_order( std::vector< int >& order  ) {
//              for (auto v : g.nodes) {
//                  // auto vd = g.neighbor(v).size();
//                  if (vd >= buckets.size())
//                      buckets.resize(vd + 1);
//                  buckets[vd].push_front(v);
//                  degrees[v] = vd;
//                  iterators[v] = buckets[vd].begin();
//              }
//
//              ordered.clear();
//     while (true) {
//         size_t i{0};
//         for (; i != buckets.size(); ++i)
//             if (!buckets[i].empty())
//                 break;
//         if (i == buckets.size())
//             break;
//         auto v = buckets[i].back();
//         order.push_back(v);
//         buckets[i].pop_back();
//         ordered.fast_add(v);
//         for (auto u : neighbor(v)) {
//             if (ordered.fast_contain(u))
//                 continue;
//             auto& ud = degrees[u];
//             buckets[ud].erase(iterators[u]);
//             --ud;
//             buckets[ud].push_front(u);
//             iterators[u] = buckets[ud].begin();
//         }
//     }
// }


} // namespace gc
