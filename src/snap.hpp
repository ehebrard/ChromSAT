#ifndef __VC_SNAP_HH
#define __VC_SNAP_HH

#include <assert.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

namespace snap
{

template <typename setsize, typename add_edge, typename add_node>
void read_graph(const char* fn, setsize ss, add_edge ae, add_node an)
{
    using std::cerr;
    try {
        std::ifstream ifs(fn);
        if (!ifs)
            throw std::runtime_error("Could not open file for reading");
				
				std::map<int,int> node_map;
				std::vector<std::pair<int,int>> edges;
				int num_nodes = 0;
				
        for (std::string line; getline(ifs, line);) {
            std::istringstream iss(line);

            // edge
            int x, y;
            iss >> x >> y;

            if (x == y || line[0] == '%')
                continue;
						
						
						if(node_map.find( x ) == node_map.end())
								node_map[x] = num_nodes++;
						
						if(node_map.find( y ) == node_map.end())
								node_map[y] = num_nodes++;

						std::pair<int,int> e{node_map[x],node_map[y]};
						
						edges.push_back(e);
            // ae(node_map[x], node_map[y]);
        }
				
				ss(num_nodes, edges.size());
				for(auto e : edges) 
					ae(e.first, e.second);
						
    } catch (std::exception& e) {
        std::cout.flush();
        cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
}
}

#endif
