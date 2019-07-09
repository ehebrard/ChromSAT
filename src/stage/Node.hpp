#ifndef __NODE_HPP
#define __NODE_HPP

//============================================================================//
//====// Includes //==========================================================//

#include <vector>

#include "Generique.hpp"

//============================================================================//
//====// Main Class //========================================================//

using namespace std;
class Node {

	/// ATTRIBUTES ///
	public:
	bnp::NodeStatus     _status;
	vector<int>         _nullified; // Nullified columns index
	bnp::Transform      _t;         // Transform
	int                 _u;         // First transformed node
	int                 _v;         // Second transformed node
	Node*               _parent;
	vector<Node*>       _children;
	float               _ub;        // Upper bound
	float               _lb;        // Lower bound
	vector<float>       _w;         // Weights
	vector<int>         _c;         // Columns
	
	/// CONSTRUCTORS ///
	public:
	Node(); /*
	* Create an root-node.
	*/

	Node(Node n, bnp::Transform t, int u, int v); /*
	* Create a node which is the child of n by t.
	*/
	
	/// METHODS ///
	public:
	bool operator<(const Node& other);
	bool operator<=(const Node& other);
	bool operator==(const Node& other);
	bool operator>=(const Node& other);
	bool operator>(const Node& other);

	
	/// GETTERS & SETTERS ///
	public:	
	//>> Getters
	float               getLB()         const;
	float               getUB()         const;
	vector<float>       getWeights()    const;
	vector<int>         getColumns()    const;
	vector<int>         getNullified()  const;
	Node*               getParent()     const;
	Node*               getChild(int i) const;
	int                 getNChildren()  const;
	bnp::NodeStatus     getStatus()     const;
	//>> Setters
	void                setUB(float ub);
	void                setLB(float lb);
	void                setStatus(bnp::NodeStatus ns);
	//>> Adders
	void                addNullified(int i);
	void                addChild(Node* n);
	void                addWeight(float w);
	void                addColumn(int c);
	//>> Clearers
	void                clearChild();
	//>> Updaters
	void                updateStatus();
	void                updateNullified(vector<vector<int>> columns);

	
	
};

#endif
