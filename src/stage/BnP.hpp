#ifndef __BNP_HPP
#define __BNP_HPP

//============================================================================//
//====// Includes //==========================================================//

#include <vector>
#include <string>
#include <set>
#include <utility>

#include "../ca_graph.hpp"
#include "../algorithm.hpp"
	
#include <ilcplex/ilocplex.h>

//============================================================================//
//====// Defines //===========================================================//

//>> The threshold beyond which we consider that there's no improvement.
//>> (ie: [-RC_EPS < lambda < RC_EPS] <=> [lambda = 0])
#define RC_EPS 1e-6

//============================================================================//
//====// Namespace : begin //=================================================//

using namespace std;
namespace bnp {

//============================================================================//
//====// Enumerations //======================================================//

//>> Possible format for input files
enum InFormat {
	I_TGF = 0,
	I_DOT = 1,
	I_CLQ = 2
};

//>> Possible format for output files
enum OutFormat {
	O_STD = 0,
	O_SOL = 1,
	O_DOT = 2
};

//>> Possible status of a node
enum NodeStatus {
	NS_CREATED   = 0,
	NS_SOLVED    = 1,
	NS_1C        = 2,
	NS_2C        = 3
};

//>> Possible tranformation of a node
enum Transform {
	T_NONE  = 0,
	T_LINK  = 1,
	T_MERGE = 2
};

//============================================================================//
//====// Structures //========================================================//

//>> This data-structure describes a node.
//>> Used through a shared_ptr.
struct Node {
	NodeStatus    status;    // Status of the node, see NodeStatus
	set<int>      nullified; // Index of columns forced to 0
	float         ub;        // Upper bound
	float         lb;        // Lower bound
	int           depth;     // T_LINK depth
	vector<float> weights;   // Weights of used columns
	vector<int>   columns;   // Index   of used columns
	Transform     t;         // Transformation used to obtain this node
	int           u;         // Vertex used for |^|
	int           v;         // Vertex used for |^|
	int           su;        // Selected vertex to generate children
	int           sv;        // Selected vertex to generate children
};

//============================================================================//
//====// Typedef //===========================================================//

//>> Type of a function pointer toward a method/function which find
//>> new columns
typedef IloNumArray (*Generator)(IloNumArray, gc::ca_graph);

//>> Type of a function pointer toward a method/function which find
//>> the branching vertice pair
typedef pair<int,int> (*Choice)(gc::ca_graph);

//>> Type of a shared pointer toward a node
typedef shared_ptr<Node> pNode; 

//============================================================================//
//====// Class //=============================================================//

class BnP {
	
///////-/// ATTRIBUTES ///-/////////////////////////////////////////////////////

	public:
	//>> About the graph itself
	string              _name;        // Name of the loaded file
	gc::ca_graph        _graph;       // Adjacency Graph of the loaded file
	vector<pNode>       _nodes;       // List of node in the current trails
	pNode               _rootNode;    // The initial node
	pNode               _currentNode; // The currently loaded node
	pNode               _incumbent;   // The best node found so far
	//>> About the solver and problem
	IloEnv          _env;
	IloModel        _masterModel;
	IloObjective    _masterObj;
	IloRangeArray   _masterXRange;  // Master range for vertices
	IloRangeArray   _masterLRange;  // Master range for columns
	IloCplex        _masterSolver;  
	IloNumVarArray  _masterLVector; // Decision vectors 
	IloNumArray     price;
	IloNumArray     column;
	//>> About column generation
	Generator           _gen;     // The pointer towards the generation function
	Choice              _choice;  // The pointer towards the choice function
	vector<vector<int>> _columns; // All the discovered columns
	


///////-/// CONSTRUCTORS ///-///////////////////////////////////////////////////

	public:
	BnP(); // Only constructor	

///////-/// METHODS ///-////////////////////////////////////////////////////////

	public:
	void load(string filename, InFormat in = I_TGF); /*
	* Load a file given its name and format as a gc::ca_graph.
	* Initialize the first node and the master problem.
	*/

	void print(bnp::OutFormat out = bnp::O_SOL); /*
	* Print the current not in _name.ext.
	* ext is determined by out.
	*/

	void solve(); /*
	* Solve the current problem and update most of the current 
	* node's values. Update the incumbent.
	*/

	void forward(int u=-1, int v=-1); /*
	* Create a child node to *_currentNode through (t,u,v). The new pNode
	* replace _currentNode and is added to _nodes. The choice of t in 
	* Tranform depends of the _currentNode status.
	* NS_CREATED or NS_SOLVED => T_LINK
	* NS_1C                   => T_MERGE
	* NS_2C                   => do nothing and put a warning in cout.
	*/

	void backward(); /*
	* Pop the _currentNode from _nodes and free the allocated memory if 
	* it's not the _incumbent or _rootNode. The _currentNode is the parent 
	* node of the previous _currentNode.
	*/

	void selectIncumbent(); /*
	* Load the incumbent as the currentNode.
	*/

	void run(); /*
	* 
	*/

	private:
	void _loadTGF(string filename); // Called by load() to perform
	void _loadDOT(string filename); // the part of file reading 
	void _loadCLQ(string filename); // given the right format.
	void _printSTD(); // Called by print() to print the
	void _printDOT(); // result in the intented format 
	void _printSOL(); // and media.


///////-/// GETTERS, SETTERS, etc ///-//////////////////////////////////////////

	//>> Setters
	public:
	void setGenerator(Generator gen);
	void setChoice(Choice choice);
	void setDiscreetMode();
	void setNoisyMode();

	//>> Getters
	NodeStatus    getCurrentStatus() const;
	float         getCurrentLB()     const;
	float         getCurrentUB()     const;
	int           getCurrentDepth()  const;
	gc::ca_graph& getRefGraph();

	//>> Adders
	void addCol(IloInt i);
	void addCol(IloNumArray ilocol);
};

//============================================================================//
//====// Operators //=========================================================//

bool operator<(Node const& a, Node const& b);
bool operator<=(Node const& a, Node const& b);
bool operator==(Node const& a, Node const& b);
bool operator>=(Node const& a, Node const& b);
bool operator>(Node const& a, Node const& b);

//============================================================================//
//====// Namespace : end //===================================================//
}

#endif
