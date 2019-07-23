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
    SN_MODE_MIN_SWAP = 0
};

//============================================================================//
//====// Structures //========================================================//

struct ms_StateSaver {
    int      u;             //>> Selected node.
    std::set<int> ps;       //>> The set of Previous Siblings.
    size_t   previous_size; //>> The size of .nodes before the transformation was applied.
};

//============================================================================//
//====// Typedef //===========================================================//

typedef std::shared_ptr<ms_StateSaver> ms_pStateSaver;

//============================================================================//
//====// Class //=============================================================//

class ms_graph : public ca_graph {
    
///////-/// ATTRIBUTES ///-/////////////////////////////////////////////////////

    private:
    std::vector<ms_pStateSaver> _ms_states;
    
///////-/// CONSTRUCTORS ///-///////////////////////////////////////////////////

    public:
    ms_graph(): ca_graph() {}
    ms_graph(int nv): ca_graph(nv) {}

///////-/// METHODS ///-////////////////////////////////////////////////////////

    public:
    int ms_select_node(const SN_MODE sn_mode) const; /*
    * >> Wrapper method <<
    * Call the appropriate methods to select a vertex according to the
    * choosen SN_MODE.
    */
    
    std::vector<int> ms_find_neighbors(const int u) const; /*
    * Find all the current neighbors of u and return them as a set.
    */
    
    bool ms_remove_node(const int u, const bool save = false); /*
    * Remove the vertex u from the .nodes intstack. If save is true then it saves
    * this modification to be able to undo it later.
    * Return true if the node was effectively removed.
    * Return false if there was no node to remove.
    */
    
    bool ms_remove_nodes(const std::vector<int> su, const bool save = false); /*
    * Call .ms_remove_node(u, save) for every u in su.
    */
    
    int ms_forward(const SN_MODE sn_mode); /*
    * Branch on the current point. Use the choosen SN_MODE for the node
    * selection.
    * Return the choosen vertex.
    */
    
    int ms_backward(); /*
    * Restore the search to the previous point.
    * Return the current vertex.
    */
    
    int ms_deep_forward(const SN_MODE sn_mode); /*
    * Branch until a max set is found.
    */
    
    int ms_backtrack(const SN_MODE sn_mode); /*
    * Call .ms_backtrack() until it is possible to call .ms_forward(sn_mode).
    * Return the choosen vertex.
    */  
    
    void ms_restore(); /*
    * Restore .nodes to its states before the maxSet search i.e. before the first
    * .ms_forward().
    */
    
    float ms_score(std::vector<float> price); /*
    * Give the score of the current set
    */
    
    std::vector<int> ms_retrieve_set(); /*
    * Retrieve the current set
    */
    
    std::vector<int> ms_find_set(const std::vector<float> price, const float obj, const SN_MODE ns_mode); /*
    * Find a set that with the given prices beats the obj.
    * if mode>0  : the set has a score higher than obj.
    * if mode<=0 : the set has a score lower than obj
    */
    
    private:
    int _ms_select_min_swap() const; /*
    * Try to find the vertex which is the closest to the limit of the intstack.
    */    
};
//============================================================================//
//====// Namespace : end //===================================================//
}

#endif
