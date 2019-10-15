#ifndef __GRAPH_HPP
#define __GRAPH_HPP

//============================================================================//
//====// Includes //==========================================================//

#include <vector>
#include <string>

#include "../ca_graph.hpp"

#include "Generique.hpp"
#include "Node.hpp"

//============================================================================//
//====// Main Class //========================================================//

using namespace std;
class Graph {

	/// ATTRIBUTES ///
	private:
	//>> About the graph itself
	string              _name;
	gc::ca_graph        _graph;         // Loaded graph
	Node*               _rootNode;      // The root node of the tree
	Node*               _currentNode;   // The current node
	Node*               _incumbent;     // Best node so far
	vector<Node>        _nodes;         // All the nodes of the tree
	vector<vector<int>> _columns;	
	//>> About the solver and problem
	IloEnv          _env;
	IloModel        _masterModel;
	IloObjective    _masterObj;
	IloRangeArray   _masterXRange;  // Master range for vertices
	IloRangeArray   _masterLRange;  // Master range for columns
	IloCplex        _masterSolver;  
	IloNumVarArray  _masterLVector; // Decision vectors 
	//>> About column generation
	bnp::Generator  _gen;

	/// CONSTRUCTORS ///
	public:
	Graph(); /*
	* Create an empty graph. Call load to fill it. 
	*/

	/// METHODS ///
	public:
	void load(string filename, bnp::InFormat in = bnp::I_TGF); /*
	* Load a file given its name and format as a gc::ca_graph.
	* Initialize the first node and the master problem.
	*/

	void print(bnp::OutFormat out = bnp::O_SOL); /*
	* Print the current incumbent in _name.ext.
	* ext is determined by out.
	*/

	void solve(); /*
	* Solve the current problem and update most of the current 
	* node's values. Update the incumbent.
	*/

	void branch(int u=-1, int v=-1); /*
	* Create the children nodes of the current node over the vertices u and v.
	* Search for appropriate vertices if vertices u or v is -1 or not given.
	*/

	void forward(); /*
	* Select the appropriate child node to solve. This node becomes the
	* current node. Do nothing if there's no child node in _children or 
	* if these are all NS_COMPLETED.
	*/
	
	void backtrack(); /*
	* Select the parent node to solve. This node becomes the
	* current node. Do nothing if the current node is the root.
	*/

	void cut(); /*
	* The current node and its direct children (if any) become
	* NS_COMPLETED.
	*/

	private:
	void _loadTGF(string filename);
	void _loadDOT(string filename);
	void _loadCLQ(string filename);
	void _printSTD();
	void _printDOT();
	void _printSOL();

	/// GETTERS & SETTERS ///
	public:
	//>> Getters
	Node*       getRoot()        const;
	Node*       getIncumbent()   const;
	Node*       getCurrentNode() const;
	string      getName()        const;
	vector<int> getCol(int i)    const;
	int         getNCol()        const;

	//>> Adders
	void addCol(IloInt i);         // <= for trivial col
	void addCol(IloNumArray col);  // <= for other

	//>> Setters
	void setMethod(bnp::Generator gen); /*
	* Set the function pointer towards the methods which find new column.
	*/

	//>> Updaters
	void updateColRange();

};

#endif
