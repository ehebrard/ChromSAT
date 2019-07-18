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

vector<int> gensolve(IloNumArray price, gc::ca_graph graph) {
	//>> Declarations
	IloInt i, j;
	vector<vector<int>> mat = graph.matrix;
	IloInt nbVertices(mat.size());
	

	
	try {	
		//>> Create an environment
		IloEnv env;
	
		//>> Output storage
		IloNumArray col(env);
		vector<int> out;

		//>> Create a model
		IloModel generatorModel = IloModel(env);

		//>> Create the objective
		IloObjective generatorObj = IloMinimize(env);
		IloAdd(generatorModel, generatorObj);

		//>> Create the decision vector
		IloNumVarArray generatorVector = IloNumVarArray(env, int(graph.size()), 0, 1, ILOINT);
		
		//>> Create constraints : "Two adjacent vertices can not belong to
		//>> the same stable set."
		int l =0;
		for (int k =0 ; k<graph.size() ; k++) {
			i = graph.nodes[k];
			if (mat[i].size() > 0) {
				for(auto j : mat[i]) {
					if ((i <= j)and(graph.nodes.contain(j))) {
						l = graph.nodes.index(j);
						generatorModel.add(   generatorVector[k]
								   +  generatorVector[l]
							           <= 1.0 );
					}
				}
			}
		}

		//>> Create constraints : "A vertex either belongs to the maximal 
		//>> stable set or at least one of its neighbors belongs to it."
		
		for (int k =0 ; k<graph.size() ; k++) {
			i = graph.nodes[k];
			IloExpr expCstr(env);
			expCstr += generatorVector[k];
			if (mat[i].size() > 0) {
				for(auto j : mat[i]) {
					if (graph.nodes.contain(j)) {
						l = graph.nodes.index(j);
						expCstr += generatorVector[l];
					}
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

		

		/*cout << "New column! :\n[" ;
		for(IloInt i=0 ; i<col.getSize() ; i++) {
			if(col[i] !=0) {
				cout << "(" << i << ":" << col[i] << ") ";
			}
		}
		cout << "]" << endl;*/

		for(int w=0 ; w<nbVertices ; w++) {
			int p = graph.parent[w];
			while ( p != graph.parent[p] ) {
				p = graph.parent[p];
			}
			out.push_back(col[graph.nodes.index(p)]);
		}

		env.end();

		//>> Output
		return out;
		
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
//====// Usage //=============================================================//

void exitWithUsage() {
	cout << wW uU "USAGE:" Uu Ww;
	cout <<          wW " bnp <filename> [-i <input format>  ||" << endl;
	cout <<       "                       -o <output format> ||" << endl;
	cout <<       "                       -d ]" Ww               << endl;
	cout << wW uU "with :" Uu Ww                                  << endl;
	cout <<    wW "<input format> among tgf, dot, clq and dimacs" << endl;
	cout <<       "<output format> among dot, sol and std"        << endl;
	cout <<       "-d limits the output printed on the shell"     << endl;
	cout Ww;
}

//============================================================================//
//====// Main //==============================================================//


int main(int argc, char * argv[]) {
	//>> Declaration
	string filename;
	InFormat in = I_TGF;
	OutFormat out = O_STD;
	bool noise = true;

	//>> Sort argv
	//>>//>> filename
	if (argc < 2) {
		exitWithUsage();
	} else {
		filename = string(argv[1]);	
	}
	//>>//>> other parameters
	for(int i=2 ; i<argc ; i++) {
		if (string(argv[i]) == "-i") {
			if (i+1 < argc) {
				i++;
				if (string(argv[i]) == "tgf") {
					in = I_TGF;
				} else if (string(argv[i]) == "dot"){
					in = I_DOT;					
				} else if (string(argv[i]) == "clq"){
					in = I_DOT;
				} else if (string(argv[i]) == "dimacs"){
					in = I_DIMACS;
				} else {
					exitWithUsage();
				}
			} else {
				exitWithUsage();
			}
		} else if (string(argv[i]) == "-o") {
			if (i+1 < argc) {
				i++;
				if (string(argv[i]) == "std") {
					out = O_STD;
				} else if (string(argv[i]) == "dot"){
					out = O_DOT;					
				} else if (string(argv[i]) == "sol"){
					out = O_SOL;
				} else {
					exitWithUsage();
				}
			} else {
				exitWithUsage();
			}

		} else if (string(argv[i]) == "-d") {
			noise = false;
		} else {
			exitWithUsage();
		}
	}
	//>> COUT
	cout << gG "Loading data, please wait." Gg << endl;


	//>> Create the BnP solver
	BnP bnp;
	bnp.load(filename, in);
	if (noise) {
		bnp.setNoisyMode();
	} else {
		bnp.setDiscreetMode();
	}


	//>> Set the modular functions
	bnp.setChoice(makechoice);
	bnp.setGenerator(gensolve);

	//>> Run the solver
	bnp.run();

	//>> Print the result
	/*bnp.selectIncumbent();
	bnp.print(out);*/
}

