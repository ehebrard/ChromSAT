//============================================================================//
//====// Includes //==========================================================//

#include "ms_graph.hpp"

//============================================================================//
//====// Namespace //=========================================================//
using namespace gc;

//============================================================================//
//====// Public Methods //====================================================//

int ms_graph::ms_select_node(const SN_MODE sn_mode) const { /// //////////// ///
    if (sn_mode == SN_MODE_MIN_SWAP) {
        return this->_ms_select_min_swap();   
    } else {
        return -1;
    }
}

std::vector<int> ms_graph::ms_find_neighbors(const int u) const { /// ////// ///
    return this->matrix[u];
}

bool ms_graph::ms_remove_node(const int u, const bool save) { /// ////////// ///    
    //>> Save
    if (save) {
        ms_pStateSaver ms_ss = std::make_shared<ms_StateSaver>();
        ms_ss->u = u;
        ms_ss->previous_size = this->nodes.size();    
        this->_ms_states.push_back(ms_ss);    
    }
    
    //>> Remove
    if (this->nodes[u] < int(this->nodes.size())) {
        this->nodes.remove(u);       
        return true;
    }
    return false;
}

bool ms_graph::ms_remove_nodes(const std::vector<int> su, const bool save) { ////////
    //>> Init
    bool out = false;
    
    //>> Remove
    for(int u : su) {
        out = out or this->ms_remove_node(u, save);    
    }
    
    //>> Return
    return out;
}

int ms_graph::ms_forward(const SN_MODE sn_mode) { /// ////////////////////// ///
    //>> Select the node
    int u = this->ms_select_node(sn_mode);
    
    if (u != -1) {
        //>> Find its neighborhood
        std::vector<int> su = this->ms_find_neighbors(u);    
        
        //>> Remove the node
        this->ms_remove_node(u, true);
        
        //>> Remove its neighborhood
        this->ms_remove_nodes(su, false);
    }
    
    //>> Return
    return u;    
}

int ms_graph::ms_deep_forward(const SN_MODE sn_mode) {
    while (this->ms_forward(sn_mode) != -1);
}

int ms_graph::ms_backward() { /// ////////////////////////////////////////// ///
    //>> Restore to previous point
    this->ms_restore();
           
    //>> Return the current node or -1 if root
    if (0 == int(this->_ms_states.size())) {
        return -1;
    } else {
        return this->_ms_states.back()->u;
    }
}

int ms_graph::ms_backtrack(const SN_MODE sn_mode) { /// ////////////////////////////// ///
    int  buffer;
    int  select;
    while ((select==-1)and(buffer!=-1)) {
        buffer = this->ms_backward();
        select = this->ms_forward(sn_mode);
    }
    return select;    
}


void ms_graph::ms_restore() { /// ////////////////////////////////////////// ///
    if (int(this->_ms_states.size()) != 0) {
        //>> Call restore
        this->nodes.restore(this->_ms_states.back()->previous_size);
                
        //>> Delete saved state
        this->_ms_states.pop_back();
    }
}

std::vector<int> ms_graph::ms_find_set(const std::vector<float> price, const float obj, const SN_MODE sn_mode) {
    //>> Init
    float UB;
    float ub;
    
    //>> First set
    this->ms_deep_forward(sn_mode);
    UB = this->ms_score(price);
    ub = UB;

    //>> Searching for the right set
    while (UB >= obj) {
        if (this->ms_backtrack(sn_mode) != -1) {
            this->ms_deep_forward(sn_mode);
            ub = this->ms_score(price);
            UB = std::min(ub,UB);
        } else {
            break;
        }             
    }

    //>> Return the constructed set
    return this->ms_retrieve_set();

}

float ms_graph::ms_score(std::vector<float> price) {
    float out = 0;
    for (auto s : this->_ms_states) {
        out += price[s->u];
    }
    return out;
}

std::vector<int> ms_graph::ms_retrieve_set() {
    std::vector<int> out;
    for (auto s : this->_ms_states) {
        out.push_back(s->u);
    }
    return out;
}

//============================================================================//
//====// Private Methods //===================================================//

int ms_graph::_ms_select_min_swap() const {
    int out = -1;
    for (int i=int(this->nodes.size())-1 ; i<=0 ; i--) {
        int elmnt = this->nodes[i];
        if (this->_ms_states.back()->ps.find(elmnt) != this->_ms_states.back()->ps.end()) {
            out = elmnt;
            break;            
        }    
    }
    return out;
}















