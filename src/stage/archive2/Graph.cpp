//============================================================================//
//====// Includes //==========================================================//

#include <string>
#include <vector>

#include <ilcplex/ilocplex.h>

#include "Generique.hpp"
#include "Graph.hpp"

ILOSTLBEGIN;

//============================================================================//
//====// Constructors //======================================================//

Graph::Graph() {}

//============================================================================//
//====// Methods //===========================================================//

void Graph::load(string filename, bnp::InFormat in) { //======================//
	//>> Declarations
	IloInt i;
	IloInt nbVertices(this->_graph.matrix.size());

	//>> Select format
	if (in == bnp::I_TGF) {
		this->_loadTGF(filename);
	} else if (in == bnp::I_DOT) {
		cout << "ERROR: This format (.dot) is not supported yet." << endl;
		exit(EXIT_FAILURE);
	} else if (in == bnp::I_CLQ) {
		cout << "ERROR: This format (.clq) is not supported yet." << endl;
		exit(EXIT_FAILURE);
	} else {
		cout << "ERROR: Unknown format." << endl;
		exit(EXIT_FAILURE);
	}

	//>> Create a root node
	this->_nodes = vector<Node>();
	this->_nodes.push_back(Node());
	this->_rootNode    = &(this->_nodes.back());
	this->_currentNode = &(this->_nodes.back());
	this->_incumbent   = &(this->_nodes.back());

	//>> Create the column vector
	this->_columns = vector<vector<int>>();

	//>> Create the master problem
	//>>//>> Environment
	this->_env = IloEnv();

	//>>//>> Create a model
	this->_masterModel = IloModel(this->_env);
	
	//>>//>> Create the objective
	this->_masterObj = IloMinimize(this->_env);
	IloAdd(this->_masterModel, this->_masterObj);
	
	//>>//>> Create constraints : "At least one color by vertex."
	this->_masterXRange = IloRangeArray(this->_env, nbVertices, 1, IloInfinity);
	IloAdd(this->_masterModel,this->_masterXRange);

	//>>//>> Create constraints : "A column cannot be selected more than once."
	//>>//>> NB => We don't have columns right now so it's empty.
	this->_masterLRange = IloRangeArray(this->_env);

	//>>//>> Create the decision vector
	this->_masterLVector = IloNumVarArray(this->_env);
	
	//>>//>> Create trivial columns
	for (i=0 ; i<nbVertices ; i++) {
		this->_masterLRange.add(IloRange(this->_env, 0, 1));
		this->_masterLVector.add(IloNumVar(this->_masterObj(1)       +
		                                   this->_masterXRange[i](1) +
		                                   this->_masterLRange[i](1) ));
	}

	//>>//>> Create solver
	this->_masterSolver = IloCplex(this->_masterModel);
} 

void Graph::print(bnp::OutFormat out) { //====================================//
	// #TBD
} 

void Graph::solve() { //======================================================//
	//>> Declaration
	IloInt      i;
	float       result = 0;
	IloInt      nbVertices(this->_graph.matrix.size());
	IloNumArray price(this->_env, nbVertices);
	IloNumArray column(this->_env, nbVertices);
	

	//>> Solve to optimality
	for(;;) {
		//>> Solve master problem with current patterns
		this->_masterSolver.solve();

		//>> Retrieve costs / price
		for (i=0 ; i<nbVertices ; i++) {
			price[i] = -this->_masterSolver.getDual(this->_masterXRange[i]);
		}

		//>> Call the generator method/function
		column = this->_gen(price, this->_graph); // <= Carefull here

		//>> Evaluate the improvement of this pattern :
		for (i = 0 ; i<nbVertices ; i++) {
			result += column[i]*price[i];
		}
		
		//>> Test the value if the new pattern may improve the master
		//>> pattern.
		if (result+1 > -RC_EPS) {
			break; // The best pattern is not usefull => break
		}

		//>> Add the new pattern to the master problem
		this->_masterLRange.add(IloRange(this->_env, 0, 1));
		i = this->_masterLRange.getSize()-1;
		this->_masterLVector.add(IloNumVar( this->_masterObj(1)         + 
		                                    this->_masterXRange(column) +
		                                    this->_masterLRange[i](1)   ));

	}

	//>> Retrieve solution values
	//>>//>> lb
	this->_currentNode->setLB(this->_masterSolver.getObjValue());
	//>>//>> ub, c, w
	float ub(0);
	for (IloInt j = 0; j < this->_masterLVector.getSize(); j++) {
		if (this->_masterSolver.getValue(this->_masterLVector[j]) != 0) {
			ub += 1;
			this->_currentNode->addColumn(j);
			this->_currentNode->addWeight( this->_masterSolver.getValue(
			                               this->_masterLVector[j])
			                             );
		}
   	}
} 

void Graph::branch(int u, int v) { //=========================================//
	//>> Clear
	this->_currentNode->clearChild();
	//>> Add left
	this->_nodes.push_back(Node(*(this->_currentNode), bnp::T_LINK, u, v));
	this->_currentNode->addChild(&(this->_nodes.back()));
	//>> Add right
	this->_nodes.push_back(Node(*(this->_currentNode), bnp::T_MERGE, u, v));
	this->_currentNode->addChild(&(this->_nodes.back()));
} 

void Graph::forward() { //====================================================//
	this->_currentNode->updateStatus();
	if (this->_currentNode->getStatus() != bnp::NS_COMPLETED) {	
		for(int i=0 ; i<(this->_currentNode->getNChildren()) ; i++) {
			if (this->_currentNode->getChild(i)->getStatus() != bnp::NS_COMPLETED) {
				this->_currentNode = this->_currentNode->getChild(i);
				break;
			}
		}
	}
	this->_currentNode->updateNullified(this->_columns);
	this->updateColRange();
} 

void Graph::backtrack() { //==================================================//
	this->_currentNode->updateStatus();
	if (this->_currentNode->getParent() != NULL) {
		this->_currentNode = this->_currentNode->getParent();
	}
	this->_currentNode->updateStatus();
	this->updateColRange(); // #TBD
} 

void Graph::cut() { //========================================================//
	this->_currentNode->setStatus(bnp::NS_COMPLETED);
} 


//============================================================================//
//====// Getters & Setters //=================================================//

//>> GETTERS
Node* Graph::getRoot() const { //=============================================//
	return this->_rootNode;
} 

Node* Graph::getIncumbent() const { //========================================//
	return this->_incumbent;
} 

Node* Graph::getCurrentNode() const { //======================================//
	return this->_currentNode;
} 

string Graph::getName() const { //============================================//
	return this->_name;
} 

vector<int> Graph::getCol(int i) const { //===================================//
	return this->_columns[i];
}

int Graph::getNCol() const { //===============================================//
	return this->_columns.size();
}

//>> ADDERS
void Graph::addCol(IloInt i) { //=============================================//
	vector<int> col;
	col.push_back(i);
	this->_columns.push_back(col);
}

void Graph::addCol(IloNumArray col) { //======================================//
	vector<int> c;
	for(IloInt i = 0 ; i<col.getSize() ; i++) {
		if (col[i] == 1) {
			c.push_back(i);
		}
	}
	this->_columns.push_back(c);
}

//>> SETTERS
void Graph::setMethod(bnp::Generator gen) { //================================//
	this->_gen = gen;
} 

//>> UPDATERS
void Graph::updateColRange() { //=============================================//
	//>> Reinitialize LRange
	for(int i=0 ; i<(this->_masterLRange.getSize()) ; i++) {
		this->_masterLRange[i].setUB(1);
	}

	//>> Load 0 in LRange
	for(int i : this->_currentNode->getNullified()) {
		this->_masterLRange[i].setUB(0);
	}
}


//============================================================================//
//====// LoadEXT //===========================================================//

void Graph::_loadTGF(string filename) { //====================================//
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

	//>> Read the vertices
	while (getline(f, line)) {
		if (line == "#") {
			break;
		} else {
			V++;
		}
	}

	//>> Create the ca_graph
	this->_graph = gc::ca_graph(V);

	//>> Read the edges
	while (getline(f, line)) {
		stringstream   sepline;
		sepline << line;
		sepline >> a >> b;
		this->_graph.add_edge(a,b);
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

void Graph::_loadCLQ(string filename) { //====================================//
	//#TBD
}

void Graph::_loadDOT(string filename) { //====================================//
	//#TBD
}

//============================================================================//
//====// PrintEXT //==========================================================//

void Graph::_printSTD() { //==================================================//
	// #TBD
} 

void Graph::_printDOT() { //==================================================//
	// #TBD
} 

void Graph::_printSOL() { //==================================================//
	// #TBD
} 
