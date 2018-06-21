#include <iostream>

#include "brancher.hpp"
#include "dimacs.hpp"

#include "sparse_graph.hpp"
#include "vertices_vec.hpp"

#include "mycielski.hpp"
#include "options.hpp"
#include "prop.hpp"
#include "rewriter.hpp"
#include "statistics.hpp"
#include "utils.hpp"
#include <minicsp/core/cons.hpp>
#include <minicsp/core/solver.hpp>
#include <minicsp/core/utils.hpp>
#include <time.h>

int main(int argc, char* argv[]){

/* vertices_vec tests */
/*
	gc::vertices_vec v1;
	v1.push_back(1);
	v1.push_back(2);
	v1.push_back(3);

	gc::vertices_vec v2;
	v2.push_back(4);
	v2.push_back(5);
	v2.push_back(2);
	
	gc::vertices_vec v3;
	v3.push_back(1);
	v3.push_back(2);
	v3.push_back(3);

	v1.safe_intersect_with(v2);
	v1.display();
	v2.display();
	v3.safe_union_with(v2);
	v3.display();
*/

	/**********************/
	/* Test reading graph */
	/**********************/
	
	// Command line parsing
	auto options = gc::parse(argc, argv);
	// and printing
    options.describe(std::cout);
	gc::sparse_graph g;

	// Check this
    dimacs::read_graph(options.instance_file.c_str(),
        [&](int nv, int) { g = gc::sparse_graph{nv}; }, // [&] clause de capture (variables accessibles par references)
        [&](int u, int v) {
			if (u<v) 
				g.add_edge(u - 1, v - 1);
            else if (u > v)
                g.safe_add_edge(u - 1, v - 1);
        },
        [&](int, gc::weight) {}); // Last argument ? (add node an)
    g.describe(std::cout);
	std::cout << "#edges : " << g.num_edges() << std::endl;
	//g.nodes.display(std::cout);
	//std::cout << "\n";	
	g.display_adjacency(std::cout);
	std::cout << "Density : " << g.sparsity()*10000 << "%°°" << std::endl;
	
	/**********************/
	/* Test clique search */
	/**********************/
/*
	gc::sparse_clique_finder cf{g};
	int lb{cf.find_cliques(g.nodes)};
	//cf.display();
	std::cout << "LB = " << lb << std::endl;
*/

	/**********************/
	/*   Test BK search   */
	/**********************/
/*
	gc::BronKerbosch bk{g};
	bk.find_cliques_withoutPivot(bk.actual_clique, bk.candidates, bk.banned);
	//bk.find_cliques_withoutPivot(bk.actual_clique, bk.candidates, bk.banned);
	
	// Display cliques
	bk.display();

	//bk.order_by_degree();
	
*/
	/**********************/
	/*   Test MDED	      */
	/**********************/
	
	gc::Min_deg_elim_game mdeg{g};
	mdeg.elimination_game(mdeg.reverse_ordering);

	mdeg.g_filled.describe(std::cout);
	std::cout << "Filled #edges : " << mdeg.g_filled.num_edges() << std::endl;	
	mdeg.g_filled.display_adjacency(std::cout);
	std::cout << "Filled Density : " << mdeg.g_filled.sparsity()*10000 << "%°°" << std::endl;
	
	/* Printer of vertices_vec, works
	gc::print_container<gc::vertices_vec> my_p(g.adjacency[0]);
	std::cout << my_p; 
	*/

    //printf("Execution time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	return 0;
}
