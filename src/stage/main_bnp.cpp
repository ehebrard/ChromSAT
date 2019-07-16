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

		/*cout << "Nouvelle colonne :\n[" ;
		for(IloInt i=0 ; i<col.getSize() ; i++) {
			if(col[i] !=0) {
				cout << "(" << i << ":" << col[i] << ") ";
			}
		}
		cout << "]" << endl;*/
		//>> Output
		return col;
		
	} catch (IloWrongUsage e) {
		cout << e << endl;
		exit(EXIT_FAILURE);		
	}

}

//============================================================================//
//====// Choice //============================================================//


pair<int, int> makechoice(gc::ca_graph g) {
	gc::arc a = g.any_non_edge();
	return make_pair(a[0],a[1]);
}

//============================================================================//
//====// Print //=============================================================//

void printNode(pNode n) {
	cout << gG endl;
	cout << "Size of .nullified: " << int(n->nullified.size()) << endl;
	cout << "Size of .weights  : " << int(n->weights.size())   << endl;
	cout << "Size of .columns  : " << int(n->columns.size())   << endl,
	cout Gg;

}


//============================================================================//
//====// Main //==============================================================//


void printG(gc::ca_graph g) {
	cout << endl;
	for (vector i : g.matrix) {
		cout << "[ ";
		for(int j : i) {
			cout << j << " ";
		}
		cout << "]" << endl;;
	}
	cout << endl;
}


int maieyn() {
	gc::ca_graph g(4);
	printG(g);
	g.contract(0,1);
	printG(g);
	g.addition(3,2);
	printG(g);
	g.addition(0,2);
	printG(g);
	g.addition(1,3);
	printG(g);
}

int main(int argc, char * argv[]) {

	cout << eE "1" Ee << endl;

	//>> Create the BnP solver
	BnP bnp;
	bnp.load(string(argv[1]), I_TGF);
	bnp.setNoisyMode();

	cout << eE "2" Ee << endl;

	//>> Set the modular functions
	bnp.setChoice(makechoice);
	bnp.setGenerator(gensolve);

	cout << eE "3" Ee << endl;

	//>> Run the solver
	bnp.run();

	cout << eE "4" Ee << endl;

	//>> Print the result
	bnp.selectIncumbent();
	bnp.print(O_STD);

	cout << eE "END" Ee << endl;
}

int maini(int argc, char * argv[]) {
	std::srand(std::time(nullptr));

	//>>
	int count = 10000;

	//>> Create the BnP solver
	BnP bnp;
	bnp.load("tgf/map_75_100900ppm.tgf", I_TGF);
	bnp.setDiscreetMode();

	//>> Set the modular functions
	bnp.setChoice(makechoice);
	bnp.setGenerator(gensolve);

	//>>
	for(int i=0 ; i<count ; i++) {
		if(i%100 == 0) {
			cout << i << endl;
		}
		//bnp.solve();
		bnp.forward();
		printNode(bnp._currentNode);
		cout << gG "Nb of node in memory: " << int(bnp._nodes.size()) << endl;
	}

	cout << eE "Fin pour " << count << " forward" Ee << endl;
}

