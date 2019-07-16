//============================================================================//
//====// Includes //==========================================================//

#include <limits>
#include <iostream>

#include "BnP.hpp"

//============================================================================//
//====// Define //============================================================//

#define wW "\033[33;1m"<<
#define Ww <<"\033[0m"
#define eE "\033[31;1m"<<
#define Ee <<"\033[0m"
#define uU "\033[1;4m"<<
#define Uu <<"\033[0m"
#define mM "\033[35;1m"<<
#define Mm <<"\033[0m"
#define gG "\033[32;1m"<<
#define Gg <<"\033[0m"

//============================================================================//
//====// Namespace //=========================================================//

ILOSTLBEGIN;
using namespace bnp; // We use every element of this namespace here.

//============================================================================//
//====// Constructors //======================================================//

BnP::BnP() {}

//============================================================================//
//====// Public Methods //====================================================//

void BnP::load(string filename, InFormat in) { // /////////////////////////// //
	//>> Declarations
	IloInt i;
	IloInt nbVertices;

	//>> Select format
	if (in == bnp::I_TGF) {
		this->_loadTGF(filename);
	} else if (in == bnp::I_DOT) {
		cout << eE "ERROR: This input format (.dot) is not supported yet." Ee << endl;
		exit(EXIT_FAILURE);
	} else if (in == bnp::I_CLQ) {
		cout << eE "ERROR: This input format (.clq) is not supported yet." Ee << endl;
		exit(EXIT_FAILURE);
	} else {
		cout << eE "ERROR: Unknown input format." Ee << endl;
		exit(EXIT_FAILURE);
	}

	//>> Retrieve the number of Vertices
	nbVertices = this->_graph.matrix.size();

	//>> Create a root node
	this->_rootNode         = make_shared<Node>();
	this->_rootNode->status = NS_CREATED;
	this->_rootNode->t      = T_NONE;
	this->_rootNode->u      = -1;
	this->_rootNode->v      = -1;
	this->_rootNode->su     = -1;
	this->_rootNode->sv     = -1;
	this->_rootNode->lb     = 0;
	this->_rootNode->ub     = float(this->_graph.size());

	//>> Make it _currentNode and add it to _nodes
	this->_currentNode = this->_rootNode;
	this->_nodes.push_back(this->_rootNode);
	this->_nodes.shrink_to_fit();

	//>> Create the column vector
	this->_columns = vector<vector<int>>();

	//>> Create the master problem
	//>>//>> Environment
	this->_env = IloEnv();

	//>>//>> 
	this->price = IloNumArray(this->_env,IloInt(int(this->_graph.size())));
	this->column =  IloNumArray(this->_env,IloInt(int(this->_graph.size())));

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
		this->addCol(i);
		this->_masterLRange.add(IloRange(this->_env, 0, 1));
		this->_masterLVector.add(IloNumVar(this->_masterObj(1)       +
		                                   this->_masterXRange[i](1) +
		                                   this->_masterLRange[i](1) ));
	}

	//>>//>> Create solver
	this->_masterSolver = IloCplex(this->_masterModel);
	
}

void BnP::print(OutFormat out) { // ///////////////////////////////////////// //
	//>> Check if there's something to print:
	if (this->_currentNode == nullptr) {
		cout << wW"WARNING: There's nothing to print!"Ww << endl;
		return;
	}

	//>> Redirecting to private function by format
	if (out == O_STD) {
		this->_printSTD();
	} else if (out == O_DOT) {
		cout << eE "ERROR: This output format (.dot) is not supported yet." Ee << endl;
		exit(EXIT_FAILURE);
	} else if (out == O_SOL) {
		cout << eE "ERROR: This output format (.sol) is not supported yet." Ee << endl;
		exit(EXIT_FAILURE);
	} else {
		cout << eE "ERROR: Unknown output format." Ee << endl;
		exit(EXIT_FAILURE);
	}
}

void BnP::solve() { // ////////////////////////////////////////////////////// //

	//>> Declaration
	bool        success;
	int         i;
	float       result = 0;
	IloInt      nbVertices(this->_graph.matrix.size());
	
	//>> Update the range thanks to the nullified list
	//>>//>> Clear nullifying constraints
	for( i=0 ; i < this->_masterLRange.getSize() ; i++ ) {
		this->_masterLVector[i].setUB(1);
	}

	//>>//>> Add the constraint
	for(int j : this->_currentNode->nullified) {
		this->_masterLVector[j].setUB(0);
	}
	cout << wW "Zero-Col: " << int(this->_currentNode->nullified.size()) Ww << endl;
	
	//>> Solve to optimality
	for(;;) {
		//>> Solve master problem with current patterns
		success = this->_masterSolver.solve();
		if (!success) {
			cout << wW "Solve() failed" Ww << endl;			
			for(IloInt l=0 ; l<this->_masterLVector.getSize() ; l++) {
				cout << "Col" << l << " in ";
				cout << "[" << this->_masterLVector[l].getLB();
				cout << ";" << this->_masterLVector[l].getUB() << "]";
				cout << " is valued " << this->_masterSolver.getValue(this->_masterLVector[l]) << endl;					
			}
			exit(EXIT_FAILURE);
		} else {
			cout << wW "Solve() succeed" Ww << endl;
		}

		//>> Retrieve costs / price
		for (i=0 ; i<nbVertices ; i++) {
			price[i] = -this->_masterSolver.getDual(this->_masterXRange[i]);
		}

		//>> Call the generator method/function
		column = this->_gen(price, this->_graph); // <= Carefull here

		//>> Evaluate the improvement of this pattern :
		result = 0;
		for (i = 0 ; i<nbVertices ; i++) {
			result += column[i]*price[i];
		}
		
		//>> Test the value if the new pattern may improve the master
		//>> pattern.
		if (result+1 > -RC_EPS) {
			break; // The best pattern is not usefull => break
		}

		//>> Add the new pattern to the master problem
		this->addCol(column);
		this->_masterLRange.add(IloRange(this->_env, 0, 1));
		i = this->_masterLRange.getSize()-1;
		this->_masterLVector.add(IloNumVar(   this->_masterObj(1)         
		                                    + this->_masterXRange(column)
		                                    + this->_masterLRange[i](1)   ));

	}

	//>> Retrieve solution values
	//>>//>> lb
	this->_currentNode->lb = this->_masterSolver.getObjValue();
	//>>//>> ub, columns, weights
	float ub(0);
	for (IloInt j=0 ; j<this->_masterLVector.getSize() ; j++) {
		if (this->_masterSolver.getValue(this->_masterLVector[j]) != 0) {
			ub += 1;
			this->_currentNode->columns.push_back(j);
			this->_currentNode->weights.push_back(this->_masterSolver.getValue(this->_masterLVector[j]));
		}
	}
	this->_currentNode->ub = ub;
	//>>//>> Update the status of _currentNode
	this->_currentNode->status = NS_SOLVED;

	//>> Check if it's the new incumbent and replace it if necessary
	if (this->_currentNode->lb == this->_currentNode->ub) {
		if (   (this->_incumbent   == nullptr         )
		    or (this->_currentNode <= this->_incumbent) ) {
			this->_incumbent = this->_currentNode;
		}
	}
}

void BnP::forward(int u , int v) { // //////////////////////////////////////// //

	//>> Declaration
	Transform t;

	//>> Check if a child cannot be generated (<=> exit here)
	if (this->_currentNode->status == NS_2C) {
		cout << wW "WARNING: The current node cannot generate another child." Ww << endl;
		return;	
	}

	//>> Select the transformation and update the status
	if (this->_currentNode->status == NS_1C) { // Create second child
		t = T_LINK;
		this->_currentNode->status = NS_2C;
	} else { // Create first child
		t = T_MERGE;
		this->_currentNode->status = NS_1C;
	}

	//>> Find appropriate value for u and v if none was found before.
	if ((u==-1)or(v==-1)) {
		pair p(this->_choice(this->_graph));
		u = get<0>(p);
		v = get<1>(p);
	}

	//>> Load the vertices already used by this node to generate children
	//>> (if any)
	if ((this->_currentNode->su != -1) and (this->_currentNode->sv != -1)) {
		u = this->_currentNode->su;
		v = this->_currentNode->sv;
	} else {
		this->_currentNode->su = u;
		this->_currentNode->sv = v;
	}
	
	//>> Create the child 
	pNode child      = make_shared<Node>();
	child->status    = NS_CREATED;
	child->nullified = this->_currentNode->nullified;
	child->depth     = this->_currentNode->depth;
	child->lb        = 0;
	child->ub        = numeric_limits<float>::max();
	child->t         = t;
	child->su        = -1;
	child->sv        = -1;
	child->u = this->_currentNode->su;
	child->v = this->_currentNode->sv;
	
	//>> Add the child to _nodes
	this->_nodes.push_back(child);
	this->_nodes.shrink_to_fit();

	//>> Update _currentNode
	this->_currentNode = child;

	//>> Update the graph and the depth
	if (t == T_LINK) {
		cout << eE "FORWARD: LINK  " << child->u << " " << child->v Ee << endl;
		this->_graph.addition(child->u, child->v);		
	} else {
		cout << eE "FORWARD: MERGE " << child->u << " " << child->v Ee << endl;
		this->_graph.contract(child->u, child->v);
		this->_currentNode->depth++;
	}

	//>> Update the nullified list
	//>> See the publication of ??? for further details.
	if (t == T_LINK) {
		for (int i=this->_nbVertices ; i < int(this->_columns.size()) ; i++) {
			if (   (find(this->_columns[i].begin(), this->_columns[i].end(), child->u) != this->_columns[i].end() )
			    and(find(this->_columns[i].begin(), this->_columns[i].end(), child->v) != this->_columns[i].end() ) ) {
				child->nullified.insert(i);
			}
		}
	} else {
		for (int i=this->_nbVertices ; i < int(this->_columns.size()) ; i++) {
			if (   (find(this->_columns[i].begin(), this->_columns[i].end(), child->u) != this->_columns[i].end())
			    and(find(this->_columns[i].begin(), this->_columns[i].end(), child->v) == this->_columns[i].end()) ) {
				child->nullified.insert(i);
			} else if (   (find(this->_columns[i].begin(), this->_columns[i].end(), child->u) == this->_columns[i].end())
			           and(find(this->_columns[i].begin(), this->_columns[i].end(), child->v) != this->_columns[i].end()) ) {
				child->nullified.insert(i);
			}
		}
	}
}

void BnP::backward() { /// ///////////////////////////////////////////////// ///
	//>>
	cout << eE "BACKWARD" Ee << endl;

	//>> Stop if this node is root	
	if (this->_currentNode->t == T_NONE) {
		return;
	}

	//>> Remove (and delete if required) the current node
	this->_nodes.pop_back();
	
	//>> Select the parent node
	this->_currentNode = this->_nodes.back();
	
	

	//>> Update the graph
	this->_graph.undo();	
}


void BnP::selectIncumbent() { /// ////////////////////////////////////////// ///
	if (this->_incumbent != nullptr) {
		this->_currentNode == this->_incumbent;
	}
}

void BnP::run() { /// ////////////////////////////////////////////////////// ///
	//>> Declaration
 	float UB = float(this->_graph.size());
	float LB = 0;
	float ub;
	float lb;
	
	//>> Branch & Price
	while (LB < UB) {

		//>> Solve current node
		this->solve();
		cout << wW this->_masterSolver.getStatus() Ww << endl;
		this->print(O_STD);

		//>> Retrieve values
		lb = this->getCurrentLB();
		ub = this->getCurrentUB();

		//>> Update values
		if ((this->getCurrentDepth() == 0)and(lb > LB)) {
			LB = lb;
		}
	
		if (ub < UB) {
			UB = ub;
		}

		cout << "LB: " << LB << "  || UB: " << UB << endl;

		//>> Save if LB == UB
		if (LB >= UB) {
			this->_incumbent = this->_currentNode;
			break;
		}		

		//>> Find next edge
		pair<int,int> c = this->_choice(this->_graph);
		cout << get<0>(c) << " " <<  get<1>(c) << endl;

		//>> Select next step
		if ((UB > lb)and
		    (     (get<0>(c)!=0)
		        or(get<1>(c)!=0)
		   )) {
			this->forward(get<0>(c),get<1>(c));
		} else if (this->getCurrentDepth() != 0) {
			this->backward();
			while ((this->getCurrentStatus() == NS_2C)and(this->getCurrentDepth() != 0)) {
				this->backward();
			}
			this->forward();
		} else {
			cout << eE uU "BREAK" Uu Ee << endl;
			break;
		}
	}
}

//============================================================================//
//====// Getters, Setter, etc //==============================================//

void BnP::setGenerator(Generator gen) { /// //////////////////////////////// ///
	this->_gen = gen;
}

void BnP::setChoice(Choice choice) { /// /////////////////////////////////// ///
	this->_choice = choice;
}

void BnP::setDiscreetMode() { /// ////////////////////////////////////////// ///
	this->_masterSolver.setOut(this->_env.getNullStream());
}

void BnP::setNoisyMode() { /// ///////////////////////////////////////////// ///
	this->_masterSolver.setOut(cout);
}

NodeStatus BnP::getCurrentStatus() const { /// ///////////////////////////// ///
	return this->_currentNode->status;
}

gc::ca_graph& BnP::getRefGraph() { /// ///////////////////////////////////// ///
	return this->_graph;
}

float BnP::getCurrentLB() const { /// ////////////////////////////////////// ///
	return this->_currentNode->lb;
}

float BnP::getCurrentUB() const { /// ////////////////////////////////////// ///
	return this->_currentNode->ub;
}

int BnP::getCurrentDepth() const { /// ///////////////////////////////////// ///
	return this->_currentNode->depth;
}

void BnP::addCol(IloInt i) { /// /////////////////////////////////////////// ///
	vector<int> trivcol;
	trivcol.push_back(i);
	trivcol.shrink_to_fit();
	this->_columns.push_back(trivcol);
	this->_columns.shrink_to_fit();
}

void BnP::addCol(IloNumArray ilocol) { /// ///////////////////////////////// ///
	vector<int> col;
	for (int i=0 ; i<ilocol.getSize() ; i++) {
		if (ilocol[i] == 1) {
			col.push_back(i);
			col.shrink_to_fit();
		}
	}
	this->_columns.push_back(col);
	this->_columns.shrink_to_fit();
}



//============================================================================//
//====// Private Methods //===================================================//

void BnP::_loadTGF(string filename) { /// //////////////////////////////// ///
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
	this->_nbVertices = V;

	//>> Create the ca_graph
	this->_graph = gc::ca_graph(V);

	//>> Read the edges
	while (getline(f, line)) {
		stringstream   sepline;
		sepline << line;
		sepline >> a >> b;
		this->_graph.add_edge(a,b);
	}

	// Close file
	f.close();
	
	// Retrieve the file name.
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

void BnP::_loadDOT(string filename) { /// ////////////////////////////////// ///
}

void BnP::_loadCLQ(string filename) { /// ////////////////////////////////// ///
}

void BnP::_printSTD() { /// //////////////////////////////////////////////// ///
	//>> Shortening names
	float lb           = this->_currentNode->lb;
	float ub           = this->_currentNode->ub;
	vector<float> wgt  = this->_currentNode->weights; 
	vector<int>   colI = this->_currentNode->columns;	

	//>> Printing!
	cout << endl;
	cout << uU "Current node:" Uu << endl;
	cout << "Upper Bound: " << ub << endl;
	cout << "Lower Bound: " << lb << endl;
	cout << uU "Selected columns:" Uu << endl;
	for(int i=0 ; i<int(colI.size()) ; i++) {
		cout << uU "Column " << colI[i] Uu << "(" << wgt[i] << "x):" << endl;
		cout << "[" << this->_columns[colI[i]][0];
		for(int j=1 ; j<int(this->_columns[colI[i]].size()) ; j++) {
			cout << " " << this->_columns[colI[i]][j];
		}
		cout << "]" << endl;
	} 
	cout << endl << "===================================" << endl  << endl;
}

void BnP::_printDOT() { /// //////////////////////////////////////////////// ///
}

void BnP::_printSOL() { /// //////////////////////////////////////////////// ///
}

//============================================================================//
//====// Operators //=========================================================//

bool operator<(Node const& a, Node const& b) { /// ///////////////////////// ///
	return a.ub < b.lb;
}

bool operator<=(Node const& a, Node const& b) { /// //////////////////////// ///
	return a.ub <= b.lb;
}

bool operator==(Node const& a, Node const& b) { /// //////////////////////// ///
	return (a.ub==b.ub)and(a.lb==b.lb); 
}

bool operator>=(Node const& a, Node const& b) { /// //////////////////////// ///
	return a.lb >= b.ub;
}

bool operator>(Node const& a, Node const& b) { /// ///////////////////////// ///
	return a.lb > b.ub;
}

