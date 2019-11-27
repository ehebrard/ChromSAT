

#include "coloring_heuristic.hpp"


namespace gc
{
	
	const std::vector<int>& coloring_heuristic::get_coloring() const
	{
	    return color;
	}

	// satur[y] was d, and is now d+1
	void coloring_heuristic::move_up(const int y, const int d)
	{
	    // swap y with *last_vertex[d]
	    auto l{*last_vertex[d]};

	    rank[l] = rank[y];
	    rank[y] = last_vertex[d];

	    *rank[y] = y;
	    *rank[l] = l;

	    ++last_vertex[d];
	}

	// satur[y] was d+1, and is now d
	void coloring_heuristic::move_down(const int y, const int d)
	{

	    assert(d >= 0);
	    assert(d < last_vertex.size());
	    assert(last_vertex[d] != begin(order));

	    // swap y with *last_vertex[d]-1
	    auto l{*(--last_vertex[d])};

	    rank[l] = rank[y];
	    rank[y] = last_vertex[d];

	    *rank[y] = y;
	    *rank[l] = l;
	}
	
	void coloring_heuristic::clear()
	{
	    last_vertex.clear();
	    color.clear();
	    if (neighbor_colors.size() > 0)
	        for (auto v : order)
	            neighbor_colors[v].clear();
	    order.clear();
	}
	

} // namespace gc

