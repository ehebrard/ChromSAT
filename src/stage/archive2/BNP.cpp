//============================================================================//
//====// Includes //==========================================================//

#include "Graph.hpp"
#include "Node.hpp"
#include <string>

ILOSTLBEGIN;

//============================================================================//
//====// Minor functions //===================================================//

void usage_n_quit() {
	string usage = "Usage: bnp <inputfile> [-i <EXT> || -o <EXT>]\n    -i <EXT>: input  format with <EXT> in {TGF, DOT, CLQ}\n    -o <EXT>: output format with <EXT> in {STD, SOL, DOT}";
	cout << usage << endl;
	exit(EXIT_FAILURE);
}

//============================================================================//
//====// Main //==============================================================//

int main(int argc, char * argv[]) {
	//>> Declarations and default values
	bnp::InFormat  iF = bnp::I_TGF;
	bnp::OutFormat oF = bnp::O_STD;	
	string  filename;
	
	//>> Test the argc :
	if ((argc < 2)or(argc%2 != 0)) {
		usage_n_quit();
	}

	//>> Retrieve filename
	filename = argv[1];

	//>> Check the flags
	for(int i=2 ; i<argc ; i++) {
		if(string(argv[i]) == "-i") {
			if(string(argv[i+1]) == "TGF") {
				iF = bnp::I_TGF;
			} else if (string(argv[i+1]) == "DOT") {
				iF = bnp::I_DOT;
			} else if (string(argv[i+1]) == "CLQ") {
				iF = bnp::I_CLQ;
			} else {
				usage_n_quit();
			}
			i++;
		} else if (string(argv[i]) == "-o"){
			if(string(argv[i+1]) == "STD") {
				oF = bnp::O_STD;
			} else if (string(argv[i+1]) == "DOT") {
				oF = bnp::O_DOT;
			} else if (string(argv[i+1]) == "SOL") {
				oF = bnp::O_SOL;
			} else {
				usage_n_quit();
			}
			i++;
		} else {
			usage_n_quit();	
		}		
	}


	//>> Create and load a graph
	Graph g;
	g.load(filename, iF);

	
}
