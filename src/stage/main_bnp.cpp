//============================================================================//
//====// Include //===========================================================//

#include "BnP.hpp"

#include <cmath>

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
//====// Global //============================================================//

gc::SN_MODE snmode = gc::SN_MODE_MIN_SWAP;

//============================================================================//
//====// Generator //=========================================================//

vector<int> othergen(IloNumArray price, gc::ms_graph graph) {
	
	vector<float> input_price;
	for (IloInt i=0 ; i<price.getSize() ; i++) {
		input_price.push_back(price[i]);
	}
	vector<int> out = graph.ms_find_set(input_price, -1, snmode);
	return out;
}

vector<int> gensolve(IloNumArray price, gc::ms_graph graph) {
	//>> Declarations
	IloInt i;
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
		for (int k =0 ; k<int(graph.size()) ; k++) {
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
		//generatorSolver.setParam(IloCplex::TiLim, 300);	

		//>> Apply price to objective
		generatorObj.setLinearCoefs(generatorVector, price);

		
		//>> Solve generator problem
		generatorSolver.solve();


		//>> Retrieve pattern
		generatorSolver.getValues(col, generatorVector);

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

pair<int, int> makechoice(gc::ms_graph g, vector<set<int>> c, Node n) {
	gc::arc a = g.any_non_edge();
	return make_pair(a[0],a[1]);
}

pair<int, int> otherchoice(gc::ms_graph g, vector<set<int>> c, Node n) {

	//>> Alias
	vector<int>   valid   = n.columns;
	vector<float> weights = n.weights;
	int u = -1;
	int v = -1;

	//>> Select the most fractionnal column
	set<int> col1;
	int select = -1;
	float gap = 0.51;	
	for(int i=0 ; i<int(weights.size()) ; i++) {
		if (abs(0.5-weights[i]) < gap) {
			gap = abs(0.5-weights[i]);
			select = valid[i];
		}
	}
	col1 = c[select];

	for(int k : col1) {
		//>> Find the first row
		u = k;

		//>> Find the second row
		for(int i : valid) {
			if (c[i].find(u) != c[i].end()) {
				for(int j : c[i]) {
					if (col1.find(j) == col1.end()) {
						v = j;
						break;
					} 
				}
				if (v != -1) break;
			}	
		}
		if (v != -1) break;
	}
	
	//>> Convert indices
	while (u != g.parent[u]) {
		u = g.parent[u];	
	}
	while (v != g.parent[v]) {
		v = g.parent[v];	
	}

	//>> Return
	return make_pair(u, v);
	

}

//============================================================================//
//====// Sort //==============================================================//

bool defaultsort (gc::ms_Node n1, gc::ms_Node n2) {
	if (n1.weight == n2.weight) {
		return n1.degree < n2.degree;
	} else {
		return n1.weight < n2.weight;
	}
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
	cout <<          wW " bnp <filename> [-i  <input format>   ||" << endl;
	cout <<       "                       -o  <output format>  ||" << endl;
	cout <<       "                       -sm <selector mode>  ||" << endl;
	cout <<       "                       -gm <generator mode> ||" << endl;
	cout <<       "                       -d ]" Ww                << endl;
	cout << wW uU "with :" Uu Ww                                  << endl;
	cout <<    wW "<input format> among tgf, dot, clq and dimacs" << endl;
	cout <<       "<output format> among dot, sol and std"        << endl;
	cout <<       "<selector mode> among swap, sort and near"     << endl;
	cout <<       "<generator mode> among cplex and greedy"       << endl;
	cout <<       "-d limits the output printed on the shell"     << endl;
	cout Ww;
	exit(EXIT_FAILURE);
}

//============================================================================//
//====// Main //==============================================================//


int main(int argc, char * argv[]) {
	//>> Declaration
	string filename;
	InFormat in = I_TGF;
	OutFormat out = O_STD;
	bool noise = true;
	Generator gen = gensolve;

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

		} else if (string(argv[i]) == "-gm") {
			if (i+1 < argc) {
				i++;
				if (string(argv[i]) == "cplex") {
					gen = gensolve;
				} else if (string(argv[i]) == "greedy"){
					gen = othergen;
				} else {
					exitWithUsage();
				}
			} else {
				exitWithUsage();
			}

		} else if (string(argv[i]) == "-sm") {
			if (i+1 < argc) {
				i++;
				if (string(argv[i]) == "swap") {
					snmode = gc::SN_MODE_MIN_SWAP;
				} else if (string(argv[i]) == "sort"){
					snmode = gc::SN_MODE_SORT;
				} else if (string(argv[i]) == "near"){
					snmode = gc::SN_MODE_NEAR;
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
	cout << gG "Loading data, please wait..." Gg << endl;

	//>> Create the BnP solver
	BnP bnp;
	bnp.load(filename, in);
	if (noise) {
		bnp.setNoisyMode();
	} else {
		bnp.setDiscreetMode();
	}

	//>> Set the modular functions
	bnp.setChoice(otherchoice);
	bnp.setGenerator(gen);
	bnp.setGraphNodeSorter(defaultsort);

	//>> Run the solver
	cout << gG "Graph loaded. Searching for starting point..." Gg << endl;
	bnp.run(-1, true);

	//>> Print the result
	bnp.print(out);

	//>> Return
	return 0;
}

