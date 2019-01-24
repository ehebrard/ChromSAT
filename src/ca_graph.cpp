#include "ca_graph.hpp"



const gc::ca_graph::find(const int v) {
	auto x{v};
	while(parent[x] != x) x = parent[x];
	return x;
}