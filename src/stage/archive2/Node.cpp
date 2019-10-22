//============================================================================//
//====// Includes //==========================================================//

#include <vector>
#include <limits>

#include "Generique.hpp"
#include "Node.hpp"

using namespace std;

//============================================================================//
//====// Constructors //======================================================//

Node::Node() { //=============================================================//
	this->_status    = bnp::NS_CREATED;
	this->_nullified = vector<int>();
	this->_t         = bnp::T_NONE;
	this->_u         = -1;
	this->_v         = -1;
	this->_parent    = NULL;
	this->_children  = vector<Node*>();
	this->_ub        = numeric_limits<float>::max();
	this->_lb        = 0;
	this->_w         = vector<float>();
	this->_c         = vector<int>();
}

Node::Node(Node n, bnp::Transform t, int u, int v) { //=======================//
	this->_status    = bnp::NS_CREATED;
	this->_nullified = vector<int>(); 
	this->_t         = t;
	this->_u         = u;
	this->_v         = v;
	this->_parent    = &n;
	this->_children  = vector<Node*>();
	this->_ub        = numeric_limits<float>::max();
	this->_lb        = 0;
	this->_w         = vector<float>();
	this->_c         = vector<int>();
}

//============================================================================//
//====// Operators //=========================================================//

bool operator<(Node const& one, Node const& another) {
	return one.getUB() < another.getLB();
}

bool operator<=(Node const& one, Node const& another) {
	return one.getUB() <= another.getLB();
}

bool operator==(Node const& one, Node const& another) {
	return (one.getUB()==another.getUB())and(one.getLB()==another.getLB());
}

bool operator>=(Node const& one, Node const& another) {
	return one.getLB() >= another.getUB();
}

bool operator>(Node const& one, Node const& another) {
	return one.getLB() > another.getUB();
}

//============================================================================//
//====// Getters & Setters //=================================================//

float Node::getLB() const { //================================================//
	return this->_lb;
}

float Node::getUB() const { //================================================//
	return this->_ub;
}

vector<float> Node::getWeights() const { //===================================//
	return this->_w;
}

vector<int> Node::getColumns() const { //=====================================//
	return this->_c;
}

vector<int> Node::getNullified() const { //===================================//
	return this->_nullified;
}

Node* Node::getParent() const { //============================================//
	return this->_parent;
}

Node* Node::getChild(int i) const { //========================================//
	return this->_children[i];
}

int Node::getNChildren() const { //===========================================//
	return this->_children.size();
}

bnp::NodeStatus Node::getStatus() const { //==================================//
	return this->_status;
}

//>> Setters
void Node::setUB(float ub) { //===============================================//
	this->_ub = ub;
}

void Node::setLB(float lb) { //===============================================//
	this->_lb = lb;
}

void Node::setStatus(bnp::NodeStatus ns) { //=================================//
	this->_status = ns;
}

//>> Adders
void Node::addNullified(int i) { //===========================================//
	this->_nullified.push_back(i);
}

void Node::addChild(Node* n) { //=============================================//
	this->_children.push_back(n);
}

void Node::addWeight(float w) { //============================================//
	this->_w.push_back(w);
}

void Node::addColumn(int c) { //==============================================//
	this->_c.push_back(c);
}

//>> Clearers
void Node::clearChild() { //==================================================//
	this->_children.clear();
}

//>> Updaters
void Node::updateStatus() { //================================================//
	bnp::NodeStatus save;
	if (this->_status == bnp::NS_SOLVED) {
		if (this->_lb == this->_ub) {
			this->_status = bnp::NS_PERFECT;
		} else {
			save = this->_status;
			this->_status = bnp::NS_COMPLETED;
			for (int i=0 ; i<int(this->_children.size()) ; i++) {
				if (   (this->_children[i]->getStatus() != bnp::NS_COMPLETED)
				    and(this->_children[i]->getStatus() != bnp::NS_COMPLETED)) { 
					this->_status = save;
					break;
				}			
			}
		}
	}
}

void Node::updateNullified(vector<vector<int>> columns) { //==================//
	//>> Copy the nullified of the father node if not root
	if (this->_parent == NULL) {
		return;	
	}
	vector<int> pn = this->_parent->getNullified();
	this->_nullified.clear();
	for(auto c : pn) {
		this->_nullified.push_back(c);
	}	

	//>> Remove the bad columns for LINK
	if (this->_t == bnp::T_LINK) {
		for (int i=0 ; i < int(columns.size()) ; i++) {
			if (   (find(columns[i].begin(), columns[i].end(), this->_u) != columns[i].end() )
			    and(find(columns[i].begin(), columns[i].end(), this->_v) != columns[i].end() ) ) {
				this->_nullified.push_back(i);
			}
		}
	}

	//>> Remove the bad columns for MERGE
	if (this->_t == bnp::T_MERGE) {
		for (int i=0 ; i < int(columns.size()) ; i++) {
			if (   (find(columns[i].begin(), columns[i].end(), this->_u) != columns[i].end())
			    and(find(columns[i].begin(), columns[i].end(), this->_v) == columns[i].end()) ) {
				this->_nullified.push_back(i);
			} else if (   (find(columns[i].begin(), columns[i].end(), this->_u) == columns[i].end())
			           and(find(columns[i].begin(), columns[i].end(), this->_v) != columns[i].end()) ) {
				this->_nullified.push_back(i);
			}
		}
	}
}

