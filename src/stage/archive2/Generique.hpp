#ifndef __GENERIQUE_HPP
#define __GENERIQUE_HPP

//============================================================================//
//====// Includes & Define //=================================================//

#include <ilcplex/ilocplex.h>
#include "../ca_graph.hpp"

#define RC_EPS 1e-6

//============================================================================//
//====// Enumerations and Structures //=======================================//

namespace bnp {
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

	//>> Possible transformation for nodes
	enum Transform {
		T_NONE  = 0,
		T_LINK  = 1,
		T_MERGE = 2
	};

	//>> Status of a node
	enum NodeStatus {
		NS_CREATED   = 0, 
		NS_SOLVED    = 1, 
		NS_COMPLETED = 2, // NS_SOLVED && (_children is NS_SOLVED)
		NS_PERFECT   = 3  // NS_SOLVED && (_lb == _ub)
	};
	
	//>> Type of a function pointer toward a method/function which find
	//>> new columns
	typedef IloNumArray (*Generator)(IloNumArray, gc::ca_graph);
}

#endif
