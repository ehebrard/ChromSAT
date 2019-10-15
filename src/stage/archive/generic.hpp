#ifndef __GENERIC_HPP
#define __GENERIC_HPP

#include <ilcplex/ilocplex.h>
#include "../ca_graph.hpp"

class BP_graph;
class BP_node;
class TransformationTable;
class ColumnManager;
class Solution;

namespace bp {
	//using   Generator = IloNumArray(gc::ca_graph, IloNumArray);
	typedef IloNumArray (*Generator)(gc::ca_graph, IloNumArray);
	//    |^|New Column|^|         |^|_currentGraph, reduced|^|
	//                             |^|  cost of the master  |^|
	//                             |^|         problem      |^|	

	enum GenMode {
		GM_SELF   = 0, // Use this class methods to find a new column. (default)
		GM_EXTERN = 1  // Use a given extern function to find a new column.
	};

	enum NodeStatus {
		ND_UNINITIALIZED = 0, // The node have just been created and you should
			              // call init() next.
		ND_INITIALIZED   = 1, // The node initialiszed but not the models and
			              // solvers have to be initialize.
		ND_READY         = 2, // The node is ready. You can call solve().
		ND_SOLVED        = 3, // The problem is solved!
		ND_TERMINATED    = 4, // The environment has been destroyed. #TooReal
		ND_VALIDATED     = 5, // TERMINATED and all _child VALIDATED.
		ND_ERROR         = 6  // An error has occured at some point.
	};

	enum TransformationType {
		TT_ADD   = 0, // Add an edge
		TT_MERGE = 1  // Merge two vertices
	};

	typedef struct {
		TransformationType type;
		int                node1;
		int                node2;
		int                nodeF;
	} Transformation ;
}

#endif
