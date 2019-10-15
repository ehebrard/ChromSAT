//============================================================================//
//====// Includes //==========================================================//

#include "bp_graph.hpp"

//>> Debug colors
#define xX "\033[1;41m"<<
#define Xx <<"\033[0m"


//============================================================================//
//====// Constructors //======================================================//

BP_graph::BP_graph() {
	this->_name        = "unknown";
	this->_nodes       = vector<BP_node>();
}

//============================================================================//
//====// Methods //===========================================================//


void BP_graph::load(string filename, bp::Format ext) { //=====================//
	//>> Select the proper loader
	if (ext == bp::F_TGF) {
		this->_loadTGF(filename);
	} else if (ext == bp::F_CLQ) {
		this->_loadCLQ(filename);
	} else if (ext == bp::F_DOT) {
		this->_loadDOT(filename);
	}	
	
	//>> Init/Reinit the _nodes
	this->_nodes = vector<BP_node>();
	this->_nodes.push_back(BP_node(this));

	//>> Init the _root and _currentNode
	this->_root = &this->_nodes.back();
	this->_currentNode = this->_root;

	//>> Init the _incumbent
	//#TBD
}
//>> Due to the length of these methods the _loadEXT methods have been placed
//>> at the end of this file.

void BP_graph::run(bp::GenMode gm /*= GM_SELF*/, //===========================//
                   bp::Generator g /*= NULL*/) { 
	//#TBD
}

bp::BacktrackStatus BP_graph::backtrack() { //====================================//
	//>> Perform a backtrack action depending on the _currentNode
	//>> see bp_graph.hpp for further details.
	if (this->_currentNode == NULL) {
		cerr << "ERROR: current node is NULL. Cannot backtrack!" << endl;
		return bp::BTS_FAILURE;

	} else if (this->_currentNode->getParent() == NULL) {
		return bp::BTS_ROOT;

	} else if (this->_currentNode->getParent() != NULL) {
		this->_currentNode = this->_currentNode->getParent();
		return bp::BTS_SUCCESS;

	} else {
		cerr << "ERROR: Unknown error occured while backtracking." 
		     << endl;
		return bp::BTS_FAILURE;
	}
}

bp::ForwardStatus BP_graph::forward() { //========================================//
	//>> Perform a forward action depending on the _currentNode
	//>> see bp_graph.hpp for further details.
	if ( this->_currentNode->getStatus() == bp::ND_VALIDATED ) {
		return bp::FWS_VALIDATED;

	} else if (     (this->_currentNode->getChildren().size() == 0)
	            and (this->_currentNode->getStatus() == bp::ND_TERMINATED) ) {
		this->_currentNode->setStatus(bp::ND_VALIDATED);
		return bp::FWS_VALIDATED;
	
	} else if (     (this->_currentNode->getChildren().size() == 2)
	            and (this->_currentNode->getChildren()[0]->getStatus() == bp::ND_VALIDATED)

	            and (this->_currentNode->getChildren()[1]->getStatus() == bp::ND_VALIDATED)
	            and (this->_currentNode->getStatus() == bp::ND_TERMINATED) ) {
		this->_currentNode->setStatus(bp::ND_VALIDATED);
		return bp::FWS_VALIDATED;

	} else if (     (this->_currentNode->getChildren().size() == 2)
	            and (this->_currentNode->getChildren()[0]->getStatus() != bp::ND_VALIDATED)
	            and (this->_currentNode->getChildren()[0]->getStatus() != bp::ND_TERMINATED) ) {
		this->_currentNode = this->_currentNode->getChildren()[0];
		return bp::FWS_SUCCESS;

	} else if (     (this->_currentNode->getChildren().size() == 2)
	            and (this->_currentNode->getChildren()[1]->getStatus() != bp::ND_VALIDATED)
	            and (this->_currentNode->getChildren()[1]->getStatus() != bp::ND_TERMINATED) ) {
		this->_currentNode = this->_currentNode->getChildren()[1];
		return bp::FWS_SUCCESS;

	} else {
		cerr << "ERROR: an error as occured while forwarding." << endl;
		return bp::FWS_FAILURE;
	} 
}

void BP_graph::cut() { //=====================================================//
	this->_currentNode->setStatus(bp::ND_VALIDATED);
}

void BP_graph::evaluateCurrentNode(bp::GenMode gm /*= GM_SELF*/, //===========// 
                                   bp::Generator g /*= NULL*/) {

	cout << xX "Begin" Xx << endl;

	//>> Declaration
	Solution s;
	cout << xX "Declaration ok" Xx << endl;

	//>> Verify that _currentNode is not NULL
	if (this->_currentNode == NULL) {
		cerr << "ERROR: current node is NULL. Cannot be solved!" 
	             << endl;
		exit(EXIT_FAILURE);
	}
	cout << xX "Not NULL ok" Xx << endl;

	//>> Initialize the node
	this->_currentNode->init();
	this->_currentNode->begin();
	this->_currentNode->createMasterProblem();
	if (gm == bp::GM_SELF) {
		this->_currentNode->createGeneratorProblem();
	} else {
		/*this->_currentNode->setGenerator(g);*/ exit(-1);
	}
	cout << xX "Init node ok" Xx << endl;

	//>> Solve the node
	this->_currentNode->solve();
	cout << xX "Solve ok" Xx << endl;
	
	//>> Retrieve solution
	s = this->_currentNode->retrieveSolution(); // <= The segmentation fault is here
	                                            // <= 'cause it's not define
	cout << xX "Retrieve ok" Xx << endl;

	//>> End the node
	this->_currentNode->end();
	cout << xX "End ok" Xx << endl;

	//>> Generate children
	this->_currentNode->branch();
	cout << xX "Branch ok" Xx << endl;
}

BP_node* BP_graph::addNode(BP_node* n, const bp::Transformation t) { //=============//
	this->_nodes.push_back(BP_node(n, t));
	return &this->_nodes.back();
}

//============================================================================//
//====// Other Functions //===================================================//

//>> none!

//============================================================================//
//====// Getters & Setters //=================================================//

gc::ca_graph BP_graph::getAdjacencyGraph() const { //=========================//
	return this->_adjGraph;
}

vector<vector<int>> BP_graph::getAdjacencyMatrix() const { //=================//
	return this->_adjGraph.matrix;
}

string BP_graph::getName() const { //=========================================//
	return this->_name;
}

Solution BP_graph::getIncumbent() const { //==================================//
	return *(this->_incumbent);
}

int BP_graph::getVerticesNb() const { //======================================//
	return this->_adjGraph.size();
}

//============================================================================//
//====// LoadEXT //===========================================================//

void BP_graph::_loadTGF(string filename) { //=================================//
	//>> Declaration 
	int            V(0);     // Vertices
	string         line;     // Buffer de lecture
	vector<string> tks;      // Vecteur de token
	int            a, b;
	ifstream       f;        // Flux de lecture entrant

	//>> Open file
	f.open(filename);
	
	//>> Verify if the file is correctly opened
	if (!f.is_open()) {
		cout << "Impossible d'ouvrir le fichier en entrée." << endl;
		exit(EXIT_FAILURE);
	}

	//>> Read the veertices
	while (getline(f, line)) {
		if (line == "#") {
			break;
		} else {
			V++;
		}
	}

	//>> Create the ca_graph
	this->_adjGraph = gc::ca_graph(V);

	//>> Read the edges
	while (getline(f, line)) {
		stringstream   sepline;
		sepline << line;
		sepline >> a >> b;
		this->_adjGraph.add_edge(a,b);
	}

	// Fermer fichier
	f.close();
	
	// Séparer et sauvegarder le nom du fichier
	stringstream   sepline;
	sepline << filename;
	tks = vector<string>();
	while (getline(sepline, line, '/')) {
		tks.push_back(line);
	} 
	if (!tks.empty()){
		stringstream   sepline;
		sepline << tks.back();
		getline(sepline, this->_name, '.');
	} else {
		this->_name = "unknown";
	}
}

void BP_graph::_loadCLQ(string filename) { //=================================//
	//#TBD
}

void BP_graph::_loadDOT(string filename) { //=================================//
	//#TBD
}

//============================================================================//
//====// Main //==============================================================//

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cout << "Usage : node <file>" << endl;
		exit(EXIT_FAILURE);
	}

	BP_graph bpg = BP_graph();
	bpg.load(argv[1], bp::F_TGF);
	bpg.evaluateCurrentNode(); // <= There's a segmentation fault in there
	//Solution s = bpg.getCurrentSolution();
	//s.print(SOM_SOL);
	//s.print(SOM_DOT);
}

//============================================================================//















