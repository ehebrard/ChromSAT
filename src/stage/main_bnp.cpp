//============================================================================//
//====// Include //===========================================================//

#include "BnP.hpp"

//============================================================================//
//====// Namespace //=========================================================//

ILOSTLBEGIN;
using namespace bnp;

//============================================================================//
//====// Define //============================================================//

#define wW "\033[33;1m"<<
#define Ww <<"\033[0m"
#define eE "\033[31;1m"<<
#define Ee <<"\033[0m"
#define uU "\033[1;4m"<<
#define Uu <<"\033[0m"
#define gG "\033[32;1m"<<
#define Gg <<"\033[0m"

//============================================================================//
//====// Generator //=========================================================//

IloNumArray gensolve(IloNumArray price, gc::ca_graph graph) {
	//>> Declarations
	IloInt i, j;
	vector<vector<int>> mat = graph.matrix;
	IloInt nbVertices(mat.size());
	
	//>> Create an environment
	IloEnv env;
	
	try {		
			
		//>> Output storage
		IloNumArray col(env);

		//>> Create a model
		IloModel generatorModel = IloModel(env);

		//>> Create the objective
		IloObjective generatorObj = IloMinimize(env);
		IloAdd(generatorModel, generatorObj);

		//>> Create the decision vector
		IloNumVarArray generatorVector = IloNumVarArray(env, nbVertices, 0, 1, ILOINT);
		
		//>> Create constraints : "Two adjacent vertices can not belong to
		//>> the same stable set."
		for (i=0 ; i<nbVertices ; i++) {
			if (mat[i].size() > 0) {
				for(auto j : mat[i]) {
					if (i <= j) {
						generatorModel.add(   generatorVector[i]
								   +  generatorVector[j]
							           <= 1.0 );
					}
				}
			}
		}

		//>> Create constraints : "A vertex either belongs to the maximal 
		//>> stable set or at least one of its neighbors belongs to it."
		for (i=0 ; i<nbVertices ; i++) {
			IloExpr expCstr(env);
			expCstr += generatorVector[i];
			if (mat[i].size() > 0) {
				for(auto j : mat[i]) {
					expCstr += generatorVector[j];
				}
			}
			IloConstraint cstr = expCstr >= 1.0;
			generatorModel.add(cstr);
		}

		//>> Create solver
		IloCplex generatorSolver = IloCplex(generatorModel);
		
		generatorSolver.setOut(env.getNullStream());

		//>> Apply price to objective
		generatorObj.setLinearCoefs(generatorVector, price);

		
		//>> Solve generator problem
		generatorSolver.solve();


		//>> Retrieve pattern
		generatorSolver.getValues(col, generatorVector);

		//>> Output
		return col;
		
	} catch (IloWrongUsage e) {
		cout << e << endl;
		exit(EXIT_FAILURE);		
	}

}

//============================================================================//
//====// Choice //============================================================//

/*GLOBAL :*/gc::coloring_algorithm<gc::ca_graph>* choicer;

pair<int, int> makechoice(gc::ca_graph g) {
	int u = rand()%int(g.size());
	int v = rand()%int(g.size());
	while (find(g.matrix[u].begin(), g.matrix[u].end(), v) != g.matrix[u].end()) {
		u = rand()%int(g.size());
		v = rand()%int(g.size());
	}
	return make_pair(u,v);
}

//============================================================================//
//====// Main //==============================================================//

int main(int argc, char * argv[]) {

	//>> Create the BnP solver
	BnP bnp;
	bnp.load(string(argv[1]), I_TGF);
	bnp.setDiscreetMode();

	//>> Set the modular functions
	bnp.setChoice(makechoice);
	bnp.setGenerator(gensolve);

	//>> Run the solver
	bnp.run();

	//>> Print the result
	bnp.print(O_STD);
}

int mainee(int argc, char * argv[]) {
	//>>
	int count = 1000;

	//>> Create the BnP solver
	BnP bnp;
	bnp.load("tgf/map_75_100900ppm.tgf", I_TGF);
	bnp.setDiscreetMode();

	//>> Set the modular functions
	bnp.setChoice(makechoice);
	bnp.setGenerator(gensolve);

	//>>
	for(vector<int> v : bnp._graph.matrix) {
		cout << "[";		
		for(int i : v) {
			cout << i << " ";
		}
		cout << "]" << endl;
	}


	//>>
	for(int i=0 ; i<count ; i++) {
		if(i%100 == 0) {
			cout << i << endl;
		}
		bnp.solve();
		bnp.forward();
	}
}
