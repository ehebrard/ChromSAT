//============================================================================//
//====// Includes //==========================================================//

#include <limits>
#include <iostream>
#include <ctime>

#include <fstream>
#include <sstream>

#include "BnP.hpp"
#include "../dimacs.hpp"

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
//====// Global //============================================================//

string color[] = {"blue","red","yellow","green","purple","chocolate",
                  "lightblue2","magenta3","navy","yellow4","lightsteelblue",
                  "antiquewhite4","cyan","deeppink"};
int colorcount = 14;

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
	if (in == I_TGF) {
		this->_loadTGF(filename);
	} else if (in == I_DOT) {
		cout << eE "ERROR: This input format (.dot) is not supported yet." Ee << endl;
		exit(EXIT_FAILURE);
	} else if (in == I_CLQ) {
		cout << eE "ERROR: This input format (.clq) is not supported yet." Ee << endl;
		exit(EXIT_FAILURE);
	} else if (in == I_DIMACS){
		this->_loadDIMACS(filename);
	} else {
		cout << eE "ERROR: Unknown input format." Ee << endl;
		exit(EXIT_FAILURE);
	}

	//>> Retrieve the number of Vertices
	nbVertices = this->_graph.matrix.size();
	this->_nbVertices = nbVertices;

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
	
	//>> and make it the _incumbent
	this->_incumbent = this->_currentNode;

	//>> Create the column vector
	this->_columns = vector<set<int>>();

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

	//>> Print the graph
	/*int w = 0;
	for (auto i : this->_graph.matrix) {
		cout << w << ": [ ";
		for(int j : i) {
			cout << j << " ";
		}
		cout << "]" << endl;
		w++;
	}*/
	
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
		this->_printDOT();
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

	//>> Solve to optimality
	for(;;) {
		//>> Solve master problem with current patterns
		success = this->_masterSolver.solve();
		if (!success) {
			cout << eE "Solve() failed" Ee << endl;			
			for(IloInt l=0 ; l<this->_masterLVector.getSize() ; l++) {
				cout << "Col" << l << " in ";
				cout << "[" << this->_masterLVector[l].getLB();
				cout << ";" << this->_masterLVector[l].getUB() << "]";
				cout << " is valued " << this->_masterSolver.getValue(this->_masterLVector[l]) << endl;					
			}
			exit(EXIT_FAILURE);
		}

		//>> Retrieve costs / price
		for (i=0 ; i<nbVertices ; i++) {
			price[i] = -this->_masterSolver.getDual(this->_masterXRange[i]);
		}
		
		//>> Call the generator method/function
		vector<int> storage = this->_gen(price, this->_graph);
		for(int k=0 ; k<int(storage.size()) ; k++) {
			column[k] = storage[k];
		}

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
	/*
	//>> Check if it's the new incumbent and replace it if necessary
	if (this->_currentNode->lb == this->_currentNode->ub) {
		if (   (this->_incumbent   == nullptr         )
		    or (this->_currentNode <= this->_incumbent) ) {
			this->_incumbent = this->_currentNode;
		}
	}*/
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
		pair p(this->_choice(this->_graph, this->_columns, *(this->_currentNode)));
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

	//>> Update _currentNode
	this->_currentNode = child;

	//>> Update the graph and the depth
	if (t == T_LINK) {
		this->_graph.addition(child->u, child->v);		
	} else {
		this->_graph.contract(child->u, child->v);
		this->_currentNode->depth++;
	}

	//>> Update the nullified list
	if (t == T_LINK) {
		for (int i=this->_nbVertices ; i < int(this->_columns.size()) ; i++) {
			if (   (this->_columns[i].find(child->u) != this->_columns[i].end() )
			    and(this->_columns[i].find(child->v) != this->_columns[i].end() ) ) {
				child->nullified.insert(i);
			}
		}
	} else {
		for (int i=this->_nbVertices ; i < int(this->_columns.size()) ; i++) {
			if (   (this->_columns[i].find(child->u) != this->_columns[i].end())
			    and(this->_columns[i].find(child->v) == this->_columns[i].end()) ) {
				child->nullified.insert(i);
			} else if (   (this->_columns[i].find(child->u) == this->_columns[i].end())
			           and(this->_columns[i].find(child->v) != this->_columns[i].end()) ) {
				child->nullified.insert(i);
			}
		}
	}
}

void BnP::backward() { /// ///////////////////////////////////////////////// ///
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
		this->_currentNode = this->_incumbent;
	}
}

void BnP::run(int timelimit, bool log) { /// /////////////////////////////// ///
	//>> Declaration
 	float UB    = float(this->_graph.size());
	float LB    = 0;
	float ub;
	float lb;
	int   i     = 1;
	int   dec   = 0;
	int   time  = -1;
	bool  timer = true;
	bool  timeout = false;
	if (timelimit == -1) {
		timer = false;
	}

	cout << uU  "Loaded graph:"  Uu << endl << this->_name << endl << endl;

	clock_t c0 = clock();

	//>> Branch & Price
	while (LB < UB) {

		//>> Solve current node
		this->_masterSolver.setParam(IloCplex::TiLim, timelimit);	
		this->solve();

	
		//>> Retrieve values
		lb = this->getCurrentLB();
		ub = this->getCurrentUB();

		//>> Update values
		if ((this->getCurrentDepth() == 0)and(lb > LB)) {
			LB = lb;
		}	
		if (ub <= UB) {
			UB = ub;
		}

		//>> Save if interesting integer solution
		if ((lb == ub)and(ub < this->_incumbent->ub)) {
			this->_incumbent = this->_currentNode;
		}

		//>> Print
		/*cout << "\33[H\33[2J";
		cout << uU  "Loaded graph:"  Uu << endl << this->_name << endl << endl;
		cout << uU "Global:" Uu << endl;
		cout << "UB = " << UB << endl;
		cout << "LB = " << LB << endl;
		this->_printCurrent();
		cout << endl;
		cout << uU "Data:" Uu << endl;
		cout << "Generated columns: " << int(this->_columns.size()) << endl;
		cout << "Currently useable columns: " << int(this->_columns.size()) - int(this->_currentNode->nullified.size()) << endl;
		cout << "Explored nodes: " << i << endl;
		cout << "CPU time: " << (clock()-c0)/(CLOCKS_PER_SEC) << endl;
		cout << endl;	
		this->_printSTD(false);*/
			

		//>> Find next edge
		pair<int,int> c = this->_choice(this->_graph, this->_columns, *(this->_currentNode));

		
		//>> Verify time limit
		time = (clock()-c0)/(CLOCKS_PER_SEC);
		if ((timer)and(time>=timelimit)){
			timeout = true;
			break;
		}

		//>> Select next step
		if ((UB > lb)and
		    (     (get<0>(c)!=0)
		        or(get<1>(c)!=0)
		   )) {
			this->forward(get<0>(c),get<1>(c));
		} else if (this->getCurrentDepth() != 0) {
			dec++;
			this->backward();
			while ((this->getCurrentStatus() == NS_2C)and(this->getCurrentDepth() != 0)) {
				this->backward();
			}
			this->forward();
		} else {
			//>> Clean graph
			while(this->_currentNode != this->_rootNode) {
				this->backward();
			}
			break;
		}
		i++;
	}
	if (log) {
		time = (clock()-c0)/(CLOCKS_PER_SEC);
		this->_printLOG(UB, LB, i, dec, time, timeout);
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

gc::ms_graph& BnP::getRefGraph() { /// ///////////////////////////////////// ///
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
	set<int> trivcol;
	trivcol.insert(i);
	this->_columns.push_back(trivcol);
}

void BnP::addCol(IloNumArray ilocol) { /// ///////////////////////////////// ///
	set<int> col;
	for (int i=0 ; i<ilocol.getSize() ; i++) {
		if (ilocol[i] == 1) {
			col.insert(i);
		}
	}
	this->_columns.push_back(col);
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

	//>> Create the ms_graph
	this->_graph = gc::ms_graph(V);

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
		this->_name = tks.back();
	} else {
		this->_name = "unknown";
	}
}

void BnP::_loadDOT(string filename) { /// ////////////////////////////////// ///
}

void BnP::_loadCLQ(string filename) { /// ////////////////////////////////// ///
}

void BnP::_loadDIMACS(string filename) { /// //////////////////////////////// ///
	//>> Declaration
	string         line;
	vector<string> tks;
	
	//>> Read the graph
	dimacs::read_graph(filename.c_str(),
		[&](int nv,  int) {this->_graph = gc::ms_graph(nv);},
		[&](int u, int v) {
			if (   (u != v)
			    and(find(this->_graph.matrix[u-1].begin(), 
			             this->_graph.matrix[u-1].end(),
			             v-1) ==  this->_graph.matrix[u-1].end()) ) {
				this->_graph.add_edge(u-1, v-1);
				
			}
		},
		[&](int, gc::weight) {});
	// Retrieve the file name.
	stringstream   sepline;
	sepline << filename;
	tks = vector<string>();
	while (getline(sepline, line, '/')) {
		tks.push_back(line);
	} 
	if (!tks.empty()){
		this->_name = tks.back();
	} else {
		this->_name = "unknown";
	}

}

void BnP::_printSTD(bool full) { /// //////////////////////////////// ///
	//>> Shortening names
	float         lb    = this->_incumbent->lb;
	float         ub    = this->_incumbent->ub;
	vector<float> wgt   = this->_incumbent->weights; 
	vector<int>   colI  = this->_incumbent->columns;	
	bool          cplt  = false;
	int           limit = 0;

	//>> Full or part ?
	if ((!full)and(int(colI.size()) > 5)) {
		limit = 5;
		cplt  = true;
	} else {
		limit = colI.size();
	}

	//>> Printing!
	cout << endl;
	cout << uU "Current best integer solution:" Uu << endl;
	cout << "Upper Bound: " << ub << endl;
	cout << "Lower Bound: " << lb << endl;
	cout << uU "Selected columns:" Uu << endl;
	for(int i=0 ; i<limit ; i++) {
		cout << uU "Column " << colI[i] Uu << "(" << wgt[i] << "x):" << endl;
		cout << "[";
		for(int j : this->_columns[colI[i]]) {
			cout << " " << j;
		}
		cout << " ]" << endl;
	}
	if (cplt) {
		cout << "..." << endl;
	}
}

void BnP::_printDOT() { /// //////////////////////////////////////////////// ///
	//>> Check if printable
	if ( int(this->_incumbent->columns.size()) > colorcount) {
		cout << wW "Warning: unable to print more than " << colorcount << " colors." Ww << endl;	
		return;
	}	

	//>> Open output file
	ofstream f(this->_name+"_out.dot");
	if (!f.is_open()) {
		cout << wW "Warning: Unable to open " << this->_name << "_out.log!" Ww << endl;	
		return;
	}

	//>> Header 
	f << "graph " << this->_name << " {" << endl;
	f << "    node[style=filled];" << endl;
	
	//>> Vertices
	for (int i=0 ; i < int(this->_graph.matrix.size()) ; i++) {
		for (int c=0 ; c < int(this->_incumbent->columns.size()) ; c++) {
			if (this->_columns[this->_incumbent->columns[c]].find(i)
			       != this->_columns[this->_incumbent->columns[c]].end()) {
				f << "    " << i << "[color=" << color[c] << "];" << endl;
				break;			
			}
		}
	}

	//>> Edges
	for (int i=0 ; i < int(this->_graph.matrix.size()) ; i++) {
		for (int j : this->_graph.matrix[i]) {
			if (i < j) {
				f << "    " << i << " -- " << j  << ";" << endl;
			}
		}
	}


	//>> Footer
	f << "}" << endl;


}

void BnP::_printSOL() { /// //////////////////////////////////////////////// ///
}

void BnP::_printCurrent() { /// ///////////////////////////////////////////// ///
	//>> Shortening names
	float lb           = this->_currentNode->lb;
	float ub           = this->_currentNode->ub;
	vector<float> wgt  = this->_currentNode->weights; 
	vector<int>   colI = this->_currentNode->columns;	

	//>> Printing!
	cout << endl;
	cout << uU "Current node:" Uu << endl;
	cout << "Upper Bound = " << ub << endl;
	cout << "Lower Bound = " << lb << endl;
}

void BnP::_printLOG(float UB, float LB, int nc, int dec, int time, bool timeout) { ///
	//>> Open as append
	ofstream log;
	log.open("log.csv", ios_base::app);

	//>> Append!
	//>>//>> Name
	log << this->_name << ";";
	
	//>>//>> K
	if (UB <= LB) {
		log << UB << ";";
	} else {
		log << "?" << ";";
	}

	//>>//>> LB
	log << LB << ";";
	
	//>>//>> UB
	log << UB << ";";

	//>>//>> Best Integer coloring	
	log << this->_incumbent->ub << ";";

	//>>//>> CPU time
	log << time << ";";

	//>>//>> timeout
	log << timeout << ";";

	//>>//>> Visited nodes
	log << nc << ";";

	//>>//>> Column generated
	log << int(this->_columns.size()) << ";";

	//>>//>> Dead-end
	log << dec << endl;
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

