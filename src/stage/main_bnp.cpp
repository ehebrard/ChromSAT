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
	

	try {
		//>> Create an environment
		IloEnv env;
			
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
		
		//>> print the col
		cout << wW "";		
		for(j=0 ; j<col.getSize() ; j++) {
			cout << col[j] << " ";
		}
		cout << "" Ww << endl;

		
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
	/*cout << gG int(choicer->G.size()) Gg << endl;
	++choicer->depth;
	choicer->G.trail.push_back(choicer->N);
	gc::arc e{choicer->s.select()};
	
	return make_pair(e[0],e[1]);*/
	cout << gG "???" Gg << endl;
	gc::arc e = choicer->s.make_choice();
		cout << gG "???" Gg << endl;
	cout << gG e[0] << " " << e[1] Gg << endl;
}

//============================================================================//
//====// Main //==============================================================//

int main(int argc, char * argv[]) {

	cout << wW "1" Ww << endl;

	//>> Create the BnP solver
	BnP bnp;
	bnp.load(string(argv[1]), I_TGF);
	bnp.setNoisyMode();

	cout << wW "2" Ww << endl;

	//>> Create the choicer
	gc::statistics stat(int(bnp.getRefGraph().size()));

	cout << wW "2.5" Ww << endl;

	gc::options opt = gc::parse(argc, argv);

	cout << wW "3" Ww << endl;

	choicer = new gc::coloring_algorithm<gc::ca_graph>(bnp.getRefGraph(), stat, opt);

	cout << wW "4" Ww << endl;

	//>> Set the modular functions
	bnp.setChoice(makechoice);
	bnp.setGenerator(gensolve);

	cout << wW "5" Ww << endl;

	//>> Run the solver
	bnp.run();

	cout << wW "6" Ww << endl;

	//>> Print the result
	bnp.print(O_STD);

	cout << wW "7" Ww << endl;

	//>> Destroy the choicer
	delete choicer;

	cout << wW "8" Ww << endl;
}
