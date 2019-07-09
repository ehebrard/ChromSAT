//============================================================================//
//====// Includes //==========================================================//

#include "transform.hpp"

//============================================================================//
//====// Constructors //======================================================//

TransformationTable::TransformationTable() {}

TransformationTable::TransformationTable(bp_graph* g) { //====================//
	//>> Attribute declaration and initialization
	this->_transformations = list<bp::Transformation>();
	this->_ref             = g;
	this->_currentSize     = this->_ref->getVerticesNb();
}

TransformationTable::TransformationTable(TransformationTable tt, //===========//
					 bp::Transformation t) {
	//>> Get the graph	
	this->_ref = tt.getRef();

	//>> Copy of tt
	this->_transformations = list<bp::Transformation>(tt);

	//>> Add t
	this->add(t);
}

//============================================================================//
//====// Methods //===========================================================//

void TransformationTable::add(bp::Transformation t) { //======================//
	//>> Index validity
	if ((t.node1 >= this->_currentSize)or(t.node2 >= this->_currentSize)) {
		cerr << "Error: Tansformation index out of range.";
		exit(EXIT_FAILURE);
	}
	
	//>> Decrease currentSize if necessary
	if (t.type == bp::TT_MERGE) {
		this->_currentSize--;
	}

	//>> Append t
	this->_transformations.push_back(t);
}

list<int> TransformationTable::regenerate(list<int> col) {
	
}

//============================================================================//
//====// Getters & Setters //=================================================//

BP_graph* TransformationTable::getRef() {
	return this->_ref;
}

