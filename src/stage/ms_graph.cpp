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

//============================================================================//
//====// Public Methods //====================================================//

int ms_graph::ms_select_node(const SN_MODE sn_mode, std::vector<float> w) const { ///
    if (sn_mode == SN_MODE_MIN_SWAP) {
        return this->_ms_select_min_swap();   
    } else if (sn_mode == SN_MODE_SORT) {
        return this->_ms_select_by_sorting(w);
    } else if (sn_mode == SN_MODE_NEAR) {
        return this->_ms_select_nearest_improvement(w);
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
    this->nodes.remove(u);       
   
   
    return false;
}

bool ms_graph::ms_remove_nodes(const std::vector<int> su, const bool save) { ///
    //>> Init
    bool out = false;
    
    //>> Remove
    for(int u : su) {
        out = out or this->ms_remove_node(u, save);    
    }
    
    //>> Return
    return out;
}

int ms_graph::ms_forward(const SN_MODE sn_mode, std::vector<float> w) { /// ////////////////////// ///
    //>> Select the node
    int u = this->ms_select_node(sn_mode, w);
    
    if (u != -1) {
        //>> Add to previous siblings
        if (not(this->_ms_states.empty())) {
            this->_ms_states.back()->ps.insert(u);
        } else {
            this->_ms_root->ps.insert(u);
        }
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

int ms_graph::ms_deep_forward(const SN_MODE sn_mode, std::vector<float> price) {

    std::vector<float> w;
    for (auto a : price) {
        w.push_back(-a);
    }

    int i = this->ms_forward(sn_mode,w);
    
    while (i != -1) {
        i = this->ms_forward(sn_mode, w);
        //this->_ms_print_current_set();
        if (this->ms_score(price)  > -1-RC_EPS) {
        	i = -1;
        }
    }
    return i;
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

int ms_graph::ms_backtrack(const SN_MODE sn_mode, std::vector<float> w) { /// //////////////////// ///
    int  buffer=1;
    int  select=-1;
    while ((select==-1)and(buffer!=-1)) {
        buffer = this->ms_backward();
        select = this->ms_forward(sn_mode, w);
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
    std::vector<int> out;
    int revert = 0;    

    std::vector<float> w;
    for (auto a : price) {
        w.push_back(-a);
    }

    //>> First set
    this->ms_deep_forward(sn_mode, price);
    UB = this->ms_score(price);
    ub = UB;
    //>> Searching for the right set
    while (UB >= -1) {
            //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            revert = this->ms_backtrack(sn_mode, w);
        if (revert != -1) {
            this->ms_deep_forward(sn_mode, price);
            ub = this->ms_score(price);
            UB = std::min(ub,UB);
        } else {
            break;
        }          
    }

    //>> Return the constructed set
    out = this->ms_retrieve_set();
    while (!this->_ms_states.empty()) {
        this->ms_restore();
    }
    return out;
}

float ms_graph::ms_score(std::vector<float> price) { /// /////////////////// ///
    float out = 0;
    for (auto s : this->_ms_states) {
        out += price[s->u];
    }
    out -= this->ms_cplt(price);
    return out;
}

std::vector<int> ms_graph::ms_retrieve_set() { /// ///////////////////////// ///
    std::vector<int> out(int(this->matrix.size()), 0);
    for (auto s : this->_ms_states) {
        out[s->u] = 1;
    }
    return out;
}

float ms_graph::ms_cplt(const std::vector<float> price) {
    
    //>> Init
    std::vector<float> w;
    for (auto a : price) {
        w.push_back(-a);
    }
    float            nbreCol   = 0;
    size_t           restorePt = this->nodes.size();
    std::vector<int> C, K;
    std::vector<int> N ;
    int              u         = -1;
    int              v         = -1;
    int              mu        = -1;
    int              mv        = -1;
	
    //>> Do the thing
    while (int(this->nodes.size())>0) {
        //>> Empty some vectors
        K.clear();
        C.clear();

        //>> Select a starting node
        v = this->nodes[0];
        K.push_back(v);        
        
        //>> Restrain its neighborhood to valid nodes
        C = this->ms_find_neighbors(v);
        for (int i=0 ; i<int(C.size()) ; i++) {
            if (not(this->nodes.contain(C[i]))) {
                C.erase(std::remove(C.begin(), C.end(), C[i]), C.end());
                i--;
            }
        }
        
        //>> Initialize the max indices.
        mv = v;
        mu = mv;

        //>> Complete the set
        
        while (not(C.empty())) {
            //>> Select 
            u = C[0];
            K.push_back(u);
            C.erase(std::remove(C.begin(), C.end(), u), C.end());

            //>> Intersection          
            N = this->ms_find_neighbors(u);   
            for (int i=0 ; i<int(C.size()) ; i++) {
                if (std::find(N.begin(), N.end(), C[i]) == N.end()) {
                    C.erase(std::remove(C.begin(), C.end(), C[i]), C.end());
                    i--;
                }
            }
            
            //>> Update max
            if (w[u] >= w[v]) {
                mu = mv;
                mv = u;
            } else if (int(K.size()) == 2) {
                mu = u;
            }
        }
        
        //>> Update
        nbreCol += w[mu];
        w[mv]   -= w[mu];
        for (int k : K) {
            if (k != mv) {
                this->nodes.remove(k);
            }
        }
        if (int(K.size())==1) {
            this->nodes.remove(K[0]);
        }
    }

    //>> Restore & Return
    this->nodes.restore(restorePt);
    return nbreCol;
}

void ms_graph::_ms_print_current_set() { /// /////////////////////////////// ///
    std::cout << "[ ";
    for (auto s : this->_ms_states) {
        std::cout << s->u << " ";
    }
    std::cout << "]" << std::endl;
}


void ms_graph::ms_set_node_sorter(NodeSorter ns) {
    this->_ms_ns = ns;
}  


//============================================================================//
//====// Private Methods //===================================================//

int ms_graph::_ms_select_min_swap() const {
    int out = -1;
    for (int i=int(this->nodes.size())-1 ; i>=0 ; i--) {
        int elmnt = this->nodes[i];
        if (not(this->_ms_states.empty())) {
            if (this->_ms_states.back()->ps.find(elmnt) == this->_ms_states.back()->ps.end()) {
                out = elmnt;
                break;  
            }  
        } else {
            if (this->_ms_root->ps.find(elmnt) == this->_ms_root->ps.end()) {
                out = elmnt;
                break;  
            } 
        }  
    }
    return out;
}


int ms_graph::_ms_select_by_sorting(std::vector<float> w) const { /// ///////////// ///
    //>> Init
    int out = -1;
  
    //>> Make the vector that will be sort
    std::vector<ms_Node> vec;
    for (int i=0 ; i<int(this->nodes.size()) ; i++) {
        int              id     = this->nodes[i];
        float            weight = w[id];
        int              deg    = 0;
        std::vector<int> N      = this->ms_find_neighbors(id);
        for (auto n : N) {
            if (this->nodes.contain(n)) {
                deg++;
            }
        }    
        ms_Node add = {id, weight, deg};
        vec.push_back(add);
    }
    
    //>> Sort it
    std::sort(vec.begin(), vec.end(),  this->_ms_ns);
    
    //>> Take the maximum
    if (not(vec.empty())) {
        for (int i=int(vec.size())-1 ; i>=0 ; i--) {
            if (not(this->_ms_states.empty())) {
                if (this->_ms_states.back()->ps.find(vec[i].id) == this->_ms_states.back()->ps.end()) {
                    out = vec[i].id;
                    break;  
                }  
            } else {
                if (this->_ms_root->ps.find(vec[i].id) == this->_ms_root->ps.end()) {
                    out = vec[i].id;
                    break;  
                } 
            }
        }
        
        //>> If it can't do the thing call someone else
        if (out == -1) {
            out = this->_ms_select_min_swap();
        }
    } 
    
    return out;

}


int ms_graph::_ms_select_nearest_improvement(std::vector<float> w) const {
    //>> Init
    int out = -1;
    
    //>> Do the thing
    for (int i=0 ; i<int(this->nodes.size()) ; i++) {
        if (w[this->nodes[i]] > 0) {
            if (not(this->_ms_states.empty())) {
                if (this->_ms_states.back()->ps.find(this->nodes[i]) == this->_ms_states.back()->ps.end()) {
                    out = this->nodes[i];
                    break;  
                }  
            } else {
                if (this->_ms_root->ps.find(this->nodes[i]) == this->_ms_root->ps.end()) {
                    out = this->nodes[i];
                    break;  
                } 
            }
        }
    } 
    
    //>> If it can't do the thing call someone else
    if (out == -1) {
        out = this->_ms_select_min_swap();
    }
    return out;
}












