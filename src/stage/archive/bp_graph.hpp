#ifndef __BP_GRAPH_HPP
#define __BP_GRAPH_HPP

//============================================================================//
//====// Includes //==========================================================//

#include <string>
#include <vector>

#include "../ca_graph.hpp"

#include "generic.hpp"
#include "bp_node.hpp"
#include "solution.hpp"

//============================================================================//
//====// Enumerations and Structures //=======================================//

namespace bp {
	enum Format {
		F_TGF = 0,
		F_CLQ = 1,
		F_DOT = 2
	};

	enum BacktrackStatus {
		BTS_SUCCESS = 0,
		BTS_ROOT    = 1,
		BTS_FAILURE = 2
	};

	enum ForwardStatus {
		FWS_SUCCESS   = 0,
		FWS_VALIDATED = 1,
		FWS_FAILURE   = 2
	};
}

using namespace std;

//============================================================================//
//====// Main Class //========================================================//

class BP_graph {

	/// ATTRIBUTES ///
	private:
	string          _name;
	gc::ca_graph    _adjGraph;    // Loaded graph
	BP_node*        _root;        // The root node of the tree
	vector<BP_node> _nodes;       // All the nodes of the tree
	BP_node*        _currentNode; // The current node
	Solution*       _incumbent;   // Best solution so far
	
	/// CONSTRUCTORS ///
	public:
	BP_graph(); /*
	* Only constructor. 
	*/

	/// METHODS ///
	public:
	void load(string filename, bp::Format ext); /*
	* >> Wrapper function <<
	*
	* Call either _loadCLQ() or _loadTGF() depending on the value of
	* ext. It only changes the read format.
	* 
	* Both wrapped functions do :
	* Create _adjGraph from the data, create the _root, empty _nodes and 
	* then add _root to it, retrieve _name and initialize the _incumbent.
	*/



	void run(bp::GenMode gm = bp::GM_SELF, bp::Generator g = NULL); /*
	* Run the branch and price algorithm over the loaded problem.
	* Use the appropriate generation method according to gm.	
	*/

	bp::BacktrackStatus backtrack(); /*
	* _ currentNode becomes the node _parent of _currentNode. If _parent is 
	* NULL then this method return BTS_ROOT else it returns BTS_SUCCESS.
	* If _currentNode is NULL when you call this function, BTS_FAILURE
	* will be returned and an error message will be displayed on the shell.
	*/

	bp::ForwardStatus forward(); /*
	* _currentNode becomes the first node of _child of the _currentNode
	* that is neither ND_VALIDATED nor ND_TERMINATED and this function 
	* return FWS_SUCCESS. 
	* If both _child are ND_VALIDATED (or the vector is empty) and 
	* _currentNode has a _status of ND_TERMINATED then _currentNode is not
	*  updated, the current node _status becomes ND_VALIDATED and this 
	* function return FWS_VALIDATED.
	* If _currentNode is already ND_VALIDATED then this methods do nothing
	* and return FWS_VALIDATED.
	* In any other case, this function is aborted and return FWS_FAILURE 
	* and an error message will be displayed on the shell.
	*/	

	void evaluateCurrentNode(bp::GenMode gm = bp::GM_SELF, bp::Generator g = NULL); /*
	* Solve the problem for the _currentNode.
	* Also do any necessary step from init() to end().
	* Use the appropriate generation method according to gm.
	*/

	void cut(); /*
	* _currentNode._status becomes ND_VALIDATED.	
	*/

	BP_node* addNode(BP_node* n, const bp::Transformation t); /*
	* Add to the list of node a child node of n by t.
	* Return a pointer to the new node.
	*/

	private:
	void _loadTGF(string filename); // see load() for more infos
	void _loadCLQ(string filename); // on these methods
	void _loadDOT(string filename); // |^/

	/// GETTERS & SETTERS ///
	public:
	gc::ca_graph        getAdjacencyGraph()  const;
	vector<vector<int>> getAdjacencyMatrix() const;
	string              getName()            const;
	Solution            getIncumbent()       const;
	int                 getVerticesNb()      const;

};

//============================================================================//

#endif
