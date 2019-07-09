#ifndef __TRANSFORM_HPP
#define __TRANSFORM_HPP

//============================================================================//
//====// Includes //==========================================================//

#include <list>
#include <ilcplex/ilocplex.h>

#include "generic.hpp"
#include "bp_graph.hpp"

//============================================================================//
//====// Enumerations and Structures //=======================================//

/* Here for documentation purpose.
 * Properly declared in generic.hpp.
 *//*
namespace bp {
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
*/
//============================================================================//
//====// Main Class //========================================================//


class TransformationTable {

	/// ATTRIBUTES ///
	private:
	std::list<bp::Transformation> _transformations;
	BP_graph*       _ref;
	std::list<int>            _initToFinal;
	std::list<std::list<int>> _finalToInit;
	int                       _currentSize;
	
	/// CONSTRUCTORS ///
	public:
	TransformationTable(BP_graph* g); /*
	* Used for root.
	*/

	TransformationTable(TransformationTable tt, bp::Transformation t); /*
	* Copy tt and add t
	*/

	/// METHODS ///
	public:
	void add(bp::Transformation t); /*
	* Add t to the table.
	*/

	std::list<int> degenerate(std::list<int> col); /*
	* Transform a valid column in the initial graph into a valid column
	* for the current transformed graph.
	*/

	std::list<int> regenerate(std::list<int> col); /*
	* Transform a valid column in the current graph into a valid column
	* for the initial graph. This new column has no guaranty to be a 
	* maximal stable set.
	*/

	bool validate(std::list<int> col); /*
	* Check if whether or not col can be used as a valid column after
	* degenerate() assuming that it is valid for the initial graph.
	*/	

	void revert(std::list<int> col);

	/// GETTERS & SETTERS ///
	public:
	BP_graph* getRef();

};

//============================================================================//

#endif
