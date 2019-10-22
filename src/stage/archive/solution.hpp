#ifndef __SOLUTION_HPP
#define __SOLUTION_HPP

//============================================================================//
//====// Includes //==========================================================//

#include <string>
#include <vector>

#include "bp_node.hpp"

using namespace std;

//============================================================================//
//====// Enumerations and Structures //=======================================//

enum SolutionOutputMode {
	SOM_DOT = 0, // Create a _name.dot file.
	SOM_SOL = 1, // Create a _name.sol file.
	SOM_STD = 2  // Print in shell.
};

//============================================================================//
//====// Main Class //========================================================//

class Solution {

	/// ATTRIBUTES ///
	private:
	string              _name;
	vector<float>       _vWeights; // List of choosen weights for each set
	vector<vector<int>> _vStables; // List of chosen stable set
	float               _lb;       // Lower bound
	float               _ub;       // Upper bound
		
	/// CONSTRUCTORS ///
	public:
	Solution(); /* 
	* Default constructor. No use so far.
	*/

	Solution(BP_node* n); /*
	* Extract a solution from a bp_node.
	*/

	/// METHODS ///
	public:
	void print(SolutionOutputMode som = SOM_DOT); /*
	* >> Wrapper function <<
	*
	* Call either _printDOT() or _printSOL() or _printSTD().
	* All of the wrapped functions print their content :
	* - SOM_DOT in a dot file with color in _tatus is OPT_INT.
	* - SOM_SOL in a sol file.
	* - SOM_STD in the shell.
	*/

	private:
	void _printDOT(); // See print() for more 
	void _printSOL(); // infos on these 3 functions.
	void _printSTD(); // |^/


	/// GETTERS & SETTERS ///
	public:
	float getLB() const;
	float getUB() const;
	
};

#endif 
