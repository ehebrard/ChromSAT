#ifndef __CG_MS_GRAPH_HH
#define __CG_MS_GRAPH_HH

//============================================================================//
//====// Includes //==========================================================//

#include <set>
#include <utility>

#include "../ca_graph.hpp"

//============================================================================//
//====// Namespace : begin //=================================================//
namespace gc
{

//============================================================================//
//====// Enumerations //======================================================//

enum SN_MODE {
    SN_MODE_MIN_SWAP = 0,
    SN_MODE_SORT = 1,
    SN_MODE_NEAR = 2
};

//============================================================================//
//====// Structures //========================================================//

struct ms_Node {
	int   id;
	float weight; //>> This value is positive. The higher, the better for our solution
	int   degree;
};

struct ms_Choice {
    int    id;
    bool   is_selected;
    size_t node_restore_point;      //>> Size of ms_graph::nodes before removing .id 
    size_t neighbors_restore_point; //>> Size of ms_graph::nodes before removing
                                    //>> the neighbors of .id
};

//============================================================================//
//====// Typedef //===========================================================//

//>> A simple shared pointer to make the memory-management easier.
typedef std::shared_ptr<ms_Choice> ms_pChoice;

//>> Function pointer towards the ordering relation.
typedef bool (* NodeSorter)(ms_Node, ms_Node);

//============================================================================//
//====// Class //=============================================================//

class ms_graph : public ca_graph {
    
///////-/// ATTRIBUTES ///-/////////////////////////////////////////////////////

    private:
    std::vector<ms_pChoice> _ms_choices;
    std::vector<int>        _ms_current_set;
    NodeSorter              _ms_ns;
    std::vector<float>      _ms_prices;
    SN_MODE                 _ms_sn_mode;
    int                     _ms_obj;
    
///////-/// CONSTRUCTORS ///-///////////////////////////////////////////////////

    public:
    ms_graph(): ca_graph() {}
    ms_graph(int nv): ca_graph(nv) {}

///////-/// METHODS ///-////////////////////////////////////////////////////////

    public:
    int ms_select_node() const; /*
    * >> Wrapper method <<
    * Call the appropriate methods to select a vertex according to the
    * attributes _ms_sn_mode.
    */
    
    std::vector<int> ms_find_neighbors(const int u) const; /*
    * Find all the current neighbors of u and return them as a vector.
    */
    
    void ms_remove_node(const int u); /*
    * Remove the vertex u from the .nodes intstack.
    */
    
    void ms_remove_nodes(const std::vector<int> su); /*
    * Call .ms_remove_node(u) for every u in su.
    */
    
    int ms_forward(); /*
    * Branch on the current point.
    * Return the choosen vertex.
    */
    
    int ms_backward(); /*
    * Restore the search to the previous point.
    * And
    * Return the current vertex.
    */
    
    void ms_deep_forward(); /*
    * Branch until a max set is found.
    */
    
    int ms_backtrack(); /*
    * Call .ms_backtrack() until it is possible to call ->ms_forward(sn_mode).
    * On finding the node where we can call ->ms_forward(), it switch the
    * selection-value (ms_Choice.is_selected) to false. It means the choosen vertex
    * will not be part of any set created in this branch of the tree.
    * Return the choosen vertex or NOTHING_FOUND.
    */  
    
    void ms_restore_to_zero(); /*
    * Restore .nodes to its states before the maxSet search i.e. before the first
    * .ms_forward().
    */
    
    float ms_score(); /*
    * Give the score of the current set.
    * True score if it is a max set.
    * Estimated score in any other case.
    */
        
    std::vector<int> ms_find_set(const std::vector<float> price); /*
    * Find a set that with the given prices beats the obj,
    * i.e. its global price is under obj.    
    */

    float ms_cplt(); /*
    * Search a lb for the remaining node (price-wise).
    */
    
    void ms_set_node_sorter(const NodeSorter ns);   
    void ms_set_sn_mode(const SN_MODE sn_mode);
    void ms_set_obj(const int obj);

    private:    
    void _ms_print_current_set();

    int _ms_select_min_swap() const; /*
    * Find the node which make the work of our intstack
    * easier.
    */   
    
    int _ms_select_by_sorting() const; /*
    * Find the maximum node by the ordering relation
    * ->_ms_node_sorter.
    */
    
    int _ms_select_nearest_improvement() const; /*
    * Find the first node with a negative price
    */
    
       
};
//============================================================================//
//====// Operator //==========================================================//

std::vector<int> inter(const std::vector<int> &a, const std::vector<int> &b) ;
/* Intersection "operator" 
 */


//============================================================================//
//====// Namespace : end //===================================================//
}

#endif
