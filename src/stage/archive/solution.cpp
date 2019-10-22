//============================================================================//
//====// Includes //==========================================================//

#include "solution.hpp"

#include <iostream>
#include <string>
#include <limits>

using namespace std;

//============================================================================//
//====// Constructors //======================================================//

Solution::Solution() { //=====================================================//
	this->_name = "unknown";
	this->_lb   = 0;
	this->_ub   = numeric_limits<float>::max();
	this->_vWeights = vector<float>();
	this->_vStables = vector<vector<int>>();
}

Solution::Solution(BP_node* n) { //============================================//
	//#URGENT #TBD
	this->_name     = n->getGraph()->getName()+".n"+to_string(n.getID());
	this->_lb       = n->getLB2(); //#TBD
	this->_ub       = n->getUB2(); //#TBD
	this->_vWeights = n->getVec(); //#TBD
	this->_vStables = n->getCol(); //#TBD
}

//============================================================================//
//====// Methods //===========================================================//

void Solution::print(SolutionOutputMode som /*= SOM_DOT*/) { //===================//
	if        (som == SOM_DOT) {
		this->_printDOT();
	} else if (som == SOM_SOL) {
		this->_printSOL();
	} else if (som == SOM_STD) {
		this->_printSTD();
	}
}

void Solution::_printSTD() { //===============================================//
	//>> Declaration	
	int i;

	//>> Print
	cout << this->_name << endl;
	cout <<    "Lower bound: " << this->_lb 
	     << " | Upper bound: " << this->_ub << endl;
	cout << "Sets:" << endl;
	for(i=0 ; i<this->_vStables.size() ; i++) {
		cout << this->_vWeights[i] << " x [";
		for(auto j : this->_vStables[i]) {
			cout << " " << j;
		}
		cout << "]" << endl;
	}
}

void Solution::_printDOT() { //===============================================//
	//#TBD
}

void Solution::_printSOL() { //===============================================//
	//#TBD
}

//============================================================================//
//====// Other functions //===================================================//



//============================================================================//
//====// Getters & Setters //=================================================//

float Solution::getLB() const {
	return this->_lb;
}

float Solution::getUB() const {
	return this->_ub;
}
