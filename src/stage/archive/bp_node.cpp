//============================================================================//
//====// Includes //==========================================================//

#include "bp_node.hpp"

ILOSTLBEGIN

//============================================================================//
//====// Constructors //======================================================//

//>> Initialization of static attribute
int BP_node::_cpt = 0;

BP_node::BP_node() {}

BP_node::BP_node(BP_graph* g) { //============================================//
	//>> Reference to the graph
	this->_graph = g;

	//>> About the node
	this->_id                  = this->_cpt;
	this->_cpt                 = this->_cpt+1;
	this->_status              = bp::ND_UNINITIALIZED;
	this->_parent              = NULL;
	this->_children            = vector<BP_node*>();
	//this->_transformationTable = NULL;
	this->_CM                  = NULL;

	//>> About the solver
	this->_currentGraph = g->getAdjacencyGraph();	

	//>> About the generator
	this->_generator = NULL;
	this->_genMode   = bp::GM_SELF;
}

BP_node::BP_node(BP_node* p, bp::Transformation t) { //=======================//

}

//============================================================================//
//====// Methods //===========================================================//

void BP_node::init() { //=====================================================//
	//#TBD
}

void BP_node::begin() { //====================================================//
	this->_env = IloEnv();
}

void BP_node::createMasterProblem() { //======================================//
	//>> Declarations
	IloInt i;
	IloInt nbVertices(this->_graph->getVerticesNb());

	//>> Create a model
	this->_masterModel = IloModel(this->_env);
	
	//>> Create the objective
	this->_masterObj = IloMinimize(this->_env);
	IloAdd(this->_masterModel, this->_masterObj);
	
	//>> Create constraints : "At least one color by vertex."
	this->_masterRange = IloRangeArray(this->_env, nbVertices, 1, IloInfinity);
	IloAdd(this->_masterModel,this->_masterRange);

	//>> Create the decision vector
	this->_masterVector = IloNumVarArray(this->_env);
	
	//>> Retrieve columns
	// #TBD
	
	//>> Create trivial columns
	for (i=0 ; i<nbVertices ; i++) {
		this->_masterVector.add(IloNumVar( this->_masterObj(1) +
		                                  this->_masterRange[i](1) ));
	}

	//>> Create solver
	this->_masterSolver = IloCplex(this->_masterModel);	
}

void BP_node::createGeneratorProblem() { //===================================//
	//>> Declarations
	IloInt i, j;
	IloInt nbVertices(this->_graph->getVerticesNb());
	vector<vector<int>> mat = this->_graph->getAdjacencyMatrix();

	//>> Create a model
	this->_generatorModel = IloModel(this->_env);

	//>> Create the objective
	this->_generatorObj = IloMinimize(this->_env);
	IloAdd(this->_generatorModel, this->_generatorObj);

	//>> Create the decision vector
	this->_generatorVector = IloNumVarArray(this->_env, nbVertices, 
	                                       0, 1, ILOINT);
	
	//>> Create constraints : "Two adjacent vertices can not belong to
	//>> the same stable set."
	for (i=0 ; i<nbVertices ; i++) {
		if (mat[i].size() > 0) {
			for(auto j : mat[i]) {
				if (i <= j) {
					this->_generatorModel.add(
					            this->_generatorVector[i]
						 +  this->_generatorVector[j]
					         <= 1.0                      );
				}
			}
		}
	}

	//>> Create constraints : "A vertex either belongs to the maximal 
	//>> stable set or at least one of its neighbors belongs to it."
	for (i=0 ; i<nbVertices ; i++) {
		IloExpr expCstr(this->_env);
		expCstr += this->_generatorVector[i];
		if (mat[i].size() > 0) {
			for(auto j : mat[i]) {
				expCstr += this->_generatorVector[j];
			}
		}
		IloConstraint cstr = expCstr >= 1.0;
		this->_generatorModel.add(cstr);
	}

	//>> Create solver
	this->_generatorSolver = IloCplex(this->_generatorModel);
	
	//>> Update genMode
	this->_genMode = bp::GM_SELF;	
}

/*
void BP_node::setGenerator(Generator g) { //==================================//
	this->_generator = g;
	this->_gemMode   = GM_EXTERN;
}
*/

void BP_node::solve() { //====================================================//
	//>> Select the solving method
	(this->_genMode == bp::GM_SELF) ? this->_selfSolve() : this->_externSolve();
}

void BP_node::_selfSolve() { //===============================================//
	//>> Declaration
	IloInt      i;
	IloInt      nbVertices(this->_graph->getVerticesNb());
	IloNumArray price(this->_env, nbVertices);
	IloNumArray column(this->_env, nbVertices);

	//>> Solve to optimality
	for(;;) {
		//>> Solve master problem with current patterns
		this->_masterSolver.solve();

		//>> Retrieve costs / price
		for (i=0 ; i<nbVertices ; i++) {
			price[i] = -this->_masterSolver.getDual(this->_masterRange[i]);
		}

		//>> Apply these to the generator problem
		this->_generatorObj.setLinearCoefs(this->_generatorVector, price);
		
		//>> Solve generator problem
		this->_generatorSolver.solve();
		
		//>> Test the value if the new pattern may improve the master
		//>> pattern.
		if (   this->_generatorSolver.getValue(this->_generatorObj)+1 
		     > -RC_EPS) {
			break; // The best pattern is not usefull => break
		}

		//>> Add the new pattern to the master problem
		this->_generatorSolver.getValues(column, this->_generatorVector);
		this->_masterVector.add(IloNumVar( this->_masterObj(1)
		                                 + this->_masterRange(column) ));

		//>> Add the new pattern to the column manager
		//#TBD
	}
}

void BP_node::_externSolve() { //=============================================//
	//#TBD
}


Solution BP_node::retrieveSolution() const { //===============================//
	//>> The work is done by a constructor of Solution	
	this->_sol = &(Solution(this)); // <= Segmentation fault here
	                                // <= Must define this constructor
	//>> Return the solution _sol
	return *(this->_sol);
}

void BP_node::end() { //======================================================//
	this->_masterModel.end();
	this->_generatorModel.end();
	this->_masterObj.end();        
	this->_generatorObj.end();  
	this->_masterRange.endElements();   
	this->_masterVector.endElements();
	this->_generatorVector.endElements();
	this->_env.end();
}

void BP_node::branch(int u /*= -1*/, int v /*= -1*/) { //=====================//
	//>> Declaration
	int min;
	
	//>> Find any missing branch node
	if ((u==-1)||(v==-1)) {
		this->_selectBranch(&u,&v);
	}

	//>> find the min index :
	(u <= v) ? (min = u) : (min = v);

	//>> Create transformations
	bp::Transformation t1 = {};		
	bp::Transformation t2 = {};
	t1.type  = bp::TT_ADD;
	t1.node1 = u;   
	t1.node2 = v;
	t1.nodeF = min;
	t2.type  = bp::TT_MERGE;
	t2.node1 = u;   
	t2.node2 = v;
	t2.nodeF = min;	

	//>> Create children
	
	this->_children = vector<BP_node*>();
	this->_children.push_back(this->_graph->addNode(this, t1));
	this->_children.push_back(this->_graph->addNode(this, t2));
}

void BP_node::_selectBranch(int* u, int* v) {
	//#TBD
}


//============================================================================//
//====// Other Functions //===================================================//

//============================================================================//
//====// Getters & Setters //=================================================//

int BP_node::getID() const {
	return this->_id;
}

float BP_node::getUB() const {
	return this->_sol->getUB();
}

float BP_node::getLB() const {
	return this->_sol->getLB();
}

BP_node* BP_node::getParent() const {
	return this->_parent;
}

vector<BP_node*> BP_node::getChildren() const {
	return this->_children;
}

bp::NodeStatus BP_node::getStatus() const {
	return this->_status;
}

BP_graph* BP_node::getGraph() const {
	return this->_graph;
}

std::vector<float> BP_node::getVec() const {
	return std::vector<float>();
}

std::vector<std::vector<int>> BP_node::getCol() const {
	return std::vector<std::vector<int>>();
}

float BP_node::getLB2() const {
	return this->_masterSolver.getObjValue();
}

float BP_node::getUB2() const {
	int cptColor(0);
	for (IloInt j = 0; j < this->_masterVector.getSize(); j++) {
		if (this->_masterSolver.getValue(this->_masterVector[j]) != 0) {
			cptColor++;
   		}
   	}
	return float(cptColor);
}

void BP_node::setStatus(bp::NodeStatus ns) {
	this->_status = ns;
}

//============================================================================//

