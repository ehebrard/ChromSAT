

#include <algorithm>
#include <cassert>
#include <iomanip>

#include "ca_graph.hpp"
// #include "graph.hpp"
// #include "dsatur.hpp"

namespace gc
{

// let's just do it greedily because I'm not sure about the objective function
// anyway
void ca_graph::get_subproblem(std::vector<int>& vertices, const int size_limit)
{
    bi_graph B;
    degeneracy_finder df(*this);
    df.degeneracy_ordering();
    cliquer cq(*this);

    gc::bitset visited(0, nodes.capacity() - 1, gc::bitset::empt);

    auto clique_sz = cq.find_cliques(rbegin(df.order), rend(df.order));

    std::vector<double> unmatched_score(cq.num_cliques, 0);
		

    for (auto i{0}; i < cq.num_cliques; ++i) {
        std::cout << std::setw(2) << cq.cliques[i].size() << " " << std::setw(3)
                  << i << " ";

        for (auto j{0}; j < cq.num_cliques; ++j) 
					if(i!=j)
				{

            B.get_from_cliques(*this, cq.cliques[i], cq.cliques[j]);
            auto mm{B.hopcroftKarp()};

            auto mm_bound
                = (cq.cliques[i].size() + cq.cliques[j].size() - B.I - mm);

            if (mm_bound - std::max(cq.cliques[i].size(), cq.cliques[j].size())
                > 0)
                std::cout //<< " " << std::setw(2)
                    << (mm_bound - std::max(cq.cliques[i].size(),
                                       cq.cliques[j].size()));
            else
                std::cout << " ";

        } else
            std::cout << " ";

        std::cout << std::endl;
    }

    // contains the cliques for which we computed the matching to every other
    // cliques
    intstack explored_cliques;
    explored_cliques.reserve(cq.num_cliques);

    //
    intstack chosen_cliques;
    chosen_cliques.reserve(cq.num_cliques);

    std::vector<std::vector<int>> matching_debt;

    std::vector<int> largest_cliques;
    for (auto i{0}; i < cq.num_cliques; ++i) {
        if (cq.cliques[i].size() == clique_sz) {
            largest_cliques.push_back(i);
        }
    }

    // while (vertices.size() < size_limit) {
    for (auto i : largest_cliques) {

        // std::cout << std::setw(2) << cq.cliques[i].size() << " " <<
        // std::setw(2)
        //           << i;

        matching_debt.resize(matching_debt.size() + 1);
        for (auto j{0}; j < cq.num_cliques; ++j) {
            if (i != j) {
                if (explored_cliques.contain(j)) {

                    matching_debt.back().push_back(
                        matching_debt[explored_cliques.index(j)][i]);
										
										// std::cout << " .";
										
                } else {
                    B.get_from_cliques(*this, cq.cliques[i], cq.cliques[j]);
                    int mm{B.hopcroftKarp()};

                    int mm_bound{static_cast<int>(cq.cliques[i].size()
                                     + cq.cliques[j].size())
                        - B.I - mm};

                    int cq_bound{static_cast<int>(
                        std::max(cq.cliques[i].size(), cq.cliques[j].size()))};

                    assert(mm_bound >= cq_bound);
                    matching_debt.back().push_back(mm_bound - cq_bound);

                    unmatched_score[i] += std::pow(
                        (double)(cq.num_cliques), (mm_bound - cq_bound));
										
										// std::cout << " *";
                }

            } else
                matching_debt.back().push_back(0);

            // std::cout << " " << std::setw(2) << matching_debt.back()[j];
        }

        explored_cliques.add(i);

        // std::cout << " -> " << unmatched_score[i] << std::endl;
    }

    // choose the best next clique
    auto choice{*std::max_element(begin(largest_cliques), end(largest_cliques),
        [&](int x, int y) { return unmatched_score[x] < unmatched_score[y]; })};

    chosen_cliques.add(choice);
    for (auto u : cq.cliques[choice]) {
        visited.add(u);
        vertices.push_back(u);
    }

    // std::cout << " -> CHOOSE " << choice << ": " << visited << std::endl;
    // for (int i{0}; i < cq.num_cliques; ++i)
    //     std::cout << std::setw(5) << i;
    // std::cout << std::endl;

    while (vertices.size() < size_limit) {

        for (auto cl{chosen_cliques.begin_not_in()};
             cl != chosen_cliques.end_not_in(); ++cl) {
            auto i{*cl};
            unmatched_score[i] = 0;

            for (auto j : chosen_cliques)
                unmatched_score[i]
                    += (double)(matching_debt[explored_cliques.index(j)][i]);
        }

        // for (int i{0}; i < cq.num_cliques; ++i)
        //     if (!chosen_cliques.contain(i))
        //         std::cout << std::setw(5) << " ";
        //     else
        //         std::cout << std::setw(5) << (int)(unmatched_score[i]);

        choice = *std::max_element(chosen_cliques.begin_not_in(),
            chosen_cliques.end_not_in(), [&](int x, int y) {
                return unmatched_score[x] < unmatched_score[y];
            });

        // std::cout << std::endl;

        // compute the matching with the chosen clique

        matching_debt.resize(matching_debt.size() + 1);
        for (int i{0}; i < cq.num_cliques; ++i)
            if (i != choice and !chosen_cliques.contain(i)) {
                if (explored_cliques.contain(i)) {

                    matching_debt.back().push_back(
                        matching_debt[explored_cliques.index(i)][choice]);
                } else {

                    B.get_from_cliques(*this, cq.cliques[choice], cq.cliques[i]);
                    int mm{B.hopcroftKarp()};

                    int mm_bound{static_cast<int>(cq.cliques[choice].size()
                                     + cq.cliques[i].size())
                        - B.I - mm};

                    int cq_bound{static_cast<int>(
                        std::max(cq.cliques[choice].size(), cq.cliques[i].size()))};

                    assert(mm_bound >= cq_bound);
                    matching_debt.back().push_back(mm_bound - cq_bound);
                }

            } else
                matching_debt.back().push_back(0);

        chosen_cliques.add(choice);
        if (!explored_cliques.contain(choice))
            explored_cliques.add(choice);
        for (auto u : cq.cliques[choice])
            if (!visited.contain(u)) {
                visited.add(u);
                vertices.push_back(u);
            }

        // std::cout << " -> CHOOSE " << choice << ": " << visited << std::endl;
    }

    // exit(1);
}

int ca_graph::find(const int v) const
{
    auto x{v};
    while (parent[x] != x)
        x = parent[x];
    return x;
}

void ca_graph::add_neighbor(const int u, const int v)
{
    matrix[u].push_back(v);
    auto vptr = rbegin(matrix[u]);
    auto stop = rend(matrix[u]);
    while (++vptr != stop) {
        if (*(vptr - 1) > *vptr)
            break;
        std::swap(*vptr, *(vptr - 1));
    }
}

void ca_graph::addition(const int u, const int v)
{
    assert(parent[u] == u);
    assert(parent[v] == v);
    add_neighbor(u, v);
    add_neighbor(v, u);

    trail.push_back(-u - 1);
    trail.push_back(-v - 1);

    ++num_edges;
}

void ca_graph::add_edge(const int u, const int v)
{
    matrix[u].push_back(v);
    matrix[v].push_back(u);
    ++num_edges;
}

void ca_graph::canonize()
{
    num_edges = 0;
    for (auto v : nodes) {
        std::sort(begin(matrix[v]), end(matrix[v]));
        matrix[v].erase(
            std::unique(begin(matrix[v]), end(matrix[v])), end(matrix[v]));
        num_edges += matrix[v].size();
    }
    num_edges /= 2;
}

void ca_graph::remove_neighbor(const int u, const int v)
{
    auto it{end(matrix[u])};
    do {
        --it;
    } while (*it != v);
    matrix[u].erase(it);
}

void ca_graph::swap_neighbor(const int u, const int v, const int w)
{
    auto r_v{gc::find_in(begin(matrix[u]), end(matrix[u]), v)};

    *r_v = w;
    if (w > v) {
        while (++r_v != end(matrix[u]) and *r_v < *(r_v - 1)) {
            std::swap(*r_v, *(r_v - 1));
        }
    } else {
        while (r_v != begin(matrix[u]) and *r_v < *(r_v - 1)) {
            std::swap(*r_v, *(r_v - 1));
            --r_v;
        }
    }
}

void ca_graph::contract(const int u, const int v)
{
    // check_consistency("before contract");

    parent[v] = u;
    rank[u] += rank[v];
    nodes.remove(v);
    nodeset.remove(v);

    union_buffer.clear();
    inter_buffer.clear();

    // keep track of the current trail's tail
    auto sz{trail.size()};

    // union the neighborhood and put N(v) \ N(u) in the trail
    gc::set_union_delta(begin(matrix[u]), end(matrix[u]), begin(matrix[v]),
        end(matrix[v]), std::back_inserter(union_buffer),
        std::back_inserter(trail), std::back_inserter(inter_buffer));

    num_edges -= inter_buffer.size();

    std::swap(matrix[u], union_buffer);

    // vertices in N(v) \ N(u)
    for (auto w{begin(trail) + sz}; w != end(trail); ++w)
        swap_neighbor(*w, v, u);

    trail.push_back(sz);
    sz = trail.size();

    // vertices in N(u) \cap N(v)
    for (auto w : inter_buffer) {
        remove_neighbor(w, v);
        trail.push_back(w);
    }

    trail.push_back(sz);
    trail.push_back(v);

    // check_consistency("after contract");
}

arc ca_graph::undo()
{
    arc pdecision;
    if (trail.back() < 0) {
        undo_addition();
    } else {
        pdecision = undo_contraction();
    }

    return pdecision;
}

arc ca_graph::backtrack(const int plevel)
{
    arc pdecision;
    while (trail.back() != plevel) {
        pdecision = undo();
    }
    trail.pop_back();

    return pdecision;
}

arc ca_graph::undo_addition()
{
    int v{-trail.back() - 1};
    trail.pop_back();
    int u{-trail.back() - 1};
    trail.pop_back();

    remove_neighbor(u, v);
    remove_neighbor(v, u);

    --num_edges;

    return arc{u, v};
}

arc ca_graph::undo_contraction()
{
    int v{trail.back()};
    trail.pop_back();
    int u{parent[v]};
    parent[v] = v;
    rank[u] -= rank[v];
    nodes.add(v);
    nodeset.add(v);

    int sz = trail.back();
    trail.pop_back();
    for (auto w{begin(trail) + sz}; w != end(trail); ++w)
        add_neighbor(*w, v);

    num_edges += (trail.size() - sz);

    trail.resize(sz);
    sz = trail.back();
    trail.pop_back();

    for (auto w{begin(trail) + sz}; w != end(trail); ++w)
        swap_neighbor(*w, u, v);

    union_buffer.clear();
    std::set_difference(begin(matrix[u]), end(matrix[u]), begin(trail) + sz,
        end(trail), std::back_inserter(union_buffer));
    std::swap(matrix[u], union_buffer);
    trail.resize(sz);

    return arc{u, v};
}

arc ca_graph::any_non_edge()
{
    for (auto u : nodes) {
        if (matrix[u].size() < nodes.size() - 1) {
            int w{0};
            while (!nodes.contain(w) or u == w)
                ++w;
            for (auto v : matrix[u]) {
                if (v > w) {
                    return arc{u, w};
                } else {
                    w = v + 1;
                    while (!nodes.contain(w) or u == w)
                        ++w;
                }
            }
        }
    }
    return arc{0, 0};
}

void ca_graph::old_search(gc::statistics& stats, gc::options& options)
{

    // exit(1);

    int limit{static_cast<int>(nodes.capacity())};

    check_consistency("beg search");

    int depth{0};

    bi_graph B;

    gc::bitset V(0, nodes.capacity() - 1, gc::bitset::empt);
    // gc::bitset ub_core(0, nodes.capacity() - 1, gc::bitset::empt);
    // gc::bitset lb_core(0, nodes.capacity() - 1, gc::bitset::empt);
    gc::bitset max_clique(0, nodes.capacity() - 1, gc::bitset::empt);

    degeneracy_finder df(*this);
    dsatur brelaz;
    brelaz.use_recolor = false;
    df.degeneracy_ordering();

    std::vector<int> degeneracy(df.degrees);
    for (auto d{begin(df.order) + 1}; d < end(df.order); ++d) {
        degeneracy[*d] = std::max(degeneracy[*(d - 1)], degeneracy[*d]);
        // std::cout << std::setw(3) << *d << " " << std::setw(3) <<
        // df.degrees[*d] << " " << std::setw(3) << degeneracy[*d] << std::endl;
    }

    std::vector<int> dg_rank(nodes.capacity());

    int i{0}, lb_frontier{0},
        ub_frontier{static_cast<int>(df.order.size()) - 1};
    for (auto u : df.order) {
        dg_rank[u] = i++;
        // if(df.degree[u] == )
    }

    cliquer cq(*this);

    // assert(df.degeneracy == *std::max_element(begin(df.degrees),
    // end(df.degrees)));

    int lb{1 + (num_edges > 0)};
    stats.notify_lb(lb);

    int ub{df.degeneracy + 1};		
    stats.notify_ub(ub);

    // std::cout << "ub = " << ub << " / " << limit << std::endl;

    std::vector<int> coloring(limit, 0);
    // int pedge[2] = {0, 0};
    bool contraction{true};
    arc pedge{0, 0};
    // decision prev{0,0};
    int nub{ub};

    int cur_lb{lb};

    int period = std::pow(10, 6 - options.verbosity);

    // auto largest_degree_criterion = [&](int x, int y) {return
    // matrix[x].size() > matrix[y].size();};
    // auto degeneracy_rank_criterion
    //     = [&](int x, int y) { return dg_rank[x] > dg_rank[y]; };
    auto global_criterion = [&](int better, int worse) {

        // 1/ start with the max clique
        auto better_in_maxclique{max_clique.contain(better)};
        auto worse_in_maxclique{max_clique.contain(worse)};
        if (better_in_maxclique and !worse_in_maxclique)
            return true;
        if (worse_in_maxclique != better_in_maxclique)
            return false;

        // 2/ priority on vertices of the maximum ub-core
        // auto better_in_ub_core{dg_rank[better] >= ub_frontier};
        // auto worse_in_ub_core{dg_rank[worse] >= ub_frontier};
        auto better_in_ub_core{degeneracy[better] >= ub};
        auto worse_in_ub_core{degeneracy[worse] >= ub};
        if (better_in_ub_core and !worse_in_ub_core)
            return true;
        if (better_in_ub_core != worse_in_ub_core)
            return false;

        // 3/ priority on vertices of the maximum lb-core
        // auto better_in_lb_core{dg_rank[better] >= lb_frontier};
        // auto worse_in_lb_core{dg_rank[worse] >= lb_frontier};
        auto better_in_lb_core{degeneracy[better] >= lb};
        auto worse_in_lb_core{degeneracy[worse] >= lb};
        if (better_in_lb_core and !worse_in_lb_core)
            return true;
        if (better_in_lb_core != worse_in_lb_core)
            return false;

        // degeneracy [STATIC]
        if (degeneracy[better] > degeneracy[worse])
            return true;
        if (degeneracy[better] < degeneracy[worse])
            return false;

        // degree [DYNAMIC]
        if (matrix[better].size() > matrix[worse].size())
            return true;

        return false;

    };

    for (auto u : nodes)
        brelaz.order.push_back(u);

    std::sort(begin(brelaz.order), end(brelaz.order), global_criterion);

    // for (auto u : brelaz.order)
    //     std::cout << std::setw(3) << u << " " << std::setw(3)
    //               << matrix[u].size() << " " << std::setw(3) << degeneracy[u]
    //               << " " << (degeneracy[u] >= ub) << (degeneracy[u] >= lb)
    //               << (max_clique.contain(u)) << "\n";

    std::vector<int> largest_cliques;

    while (lb < ub) {

        // check_consistency("search loop");
        // if (stats.total_conflicts > 100)
        //     exit(1);

        if (options.verbosity > 1
            && stats.notify_iteration(depth) % period == 0)
            stats.custom_force_display(std::cout);

        double tbefore = minicsp::cpuTime();
        cq.clear();
        auto clique_sz = cq.find_cliques(
            begin(brelaz.order), end(brelaz.order), options.cliquelimit);

        stats.notify_nclique(cq.num_cliques);
        stats.notify_clique_time(minicsp::cpuTime() - tbefore);

        // std::cout << clique_sz << std::endl;

        if (clique_sz > cur_lb) {
            cur_lb = clique_sz;
        }

        tbefore = minicsp::cpuTime();
        largest_cliques.clear();
        // largest_cliques.push_back(0);

        for (auto i{0}; i < cq.cliques.size(); ++i) {
            if (cq.cliques[i].size() >= clique_sz - 1) {
                largest_cliques.push_back(i);
            }
        }

        max_clique.clear();
        max_clique.union_with(cq.cliques[largest_cliques[0]]);

        int matching_bound = clique_sz;

        // std::cout << "     ";
        // for (auto i{1}; i < largest_cliques.size(); ++i) {
        //     std::cout << " " << std::setw(2) << i;
        // }
        // std::cout << std::endl;

        if (ub - clique_sz < 4 and largest_cliques.size() > 1) {

            for (auto i{0}; i < largest_cliques.size(); ++i) {
                auto c1{largest_cliques[i]};

                // std::cout << std::setw(2) << cq.cliques[c1].size() << " "
                //           << std::setw(2) << i;
                //
                // for (auto j{0}; j < i; ++j)
                //     std::cout << "   ";

                for (auto j{i + 1}; j < largest_cliques.size(); ++j) {
                    auto c2{largest_cliques[j]};

                    B.get_from_cliques(*this, cq.cliques[c1], cq.cliques[c2]);
                    auto mm{B.hopcroftKarp()};

                    auto mm_bound = (cq.cliques[c1].size()
                        + cq.cliques[c2].size() - B.I - mm);

                    // std::cout << " " << std::setw(2)
                    //           << (mm_bound - std::max(cq.cliques[c1].size(),
                    //                              cq.cliques[c2].size()));
                    // // << cq.cliques[c1].size() << "," <<
                    // cq.cliques[c2].size()
                    // // << ":"
                    // // << mm_bound;

                    if (mm_bound > matching_bound) {
                        matching_bound = mm_bound;
                    }
                }

                // std::cout << std::endl;
            }
            // exit(1);

            // std::cout << std::endl;

            // if(clique_sz < matching_bound-1)
            // 	std::cout << (matching_bound - clique_sz) << std::endl;

            // // if(clique_sz == matching_bound)
            // std::cout << clique_sz << " " << (matching_bound - clique_sz)
            //           << std::endl;

            stats.notify_bound_delta(clique_sz, matching_bound);
            if (cur_lb < matching_bound) {
                cur_lb = matching_bound;
            }
        }
        stats.notify_matching_time(minicsp::cpuTime() - tbefore);

        tbefore = minicsp::cpuTime();
        brelaz.clear();
        nub = brelaz.brelaz_color_score(
            *this, ub - 1, global_criterion, size(), 12345);

        if (nub != size())
            std::cout << "YEEPEE!\n";

        // cur_lb = 0;
        clique_sz = 0;
        for (auto v : brelaz.order) {
            // std::cout << v << " " << brelaz.color[v] << std::endl;
            if (brelaz.color[v] < clique_sz) {
                //
                if (brelaz.color[brelaz.order[brelaz.color[v]]]
                    != brelaz.color[v]) {
                    std::cout << "col[" << v << "] = " << brelaz.color[v]
                              << " < " << clique_sz << std::endl;
                }

                assert(brelaz.color[brelaz.order[brelaz.color[v]]]
                    == brelaz.color[v]);

                // std::cout << "--> branch on " << v << "," <<
                // brelaz.order[brelaz.color[v]] << std::endl;
                pedge = arc{v, brelaz.order[brelaz.color[v]]};
                break;
            } // else {
            // 							std::cout << "
            // "
            // <<
            // v;
            // 							max_clique.add(v);
            // 						}
            if (++clique_sz == ub)
                break;
        }
        stats.notify_dsatur_time(minicsp::cpuTime() - tbefore);

        // std::cout << clique_sz << ":";
        // for (int i = 0; i < clique_sz; ++i) {
        //     std::cout << " " << brelaz.order[i];
        // }
        //
        // std::cout << std::endl
        //           << max_clique.size() << ": " << max_clique <<
        //           std::endl;

        // std::cout << " nub = " << nub << "/" << size() << std::endl;

        // std::cout << "brelaz clique = " << clique_sz << std::endl;

        if (clique_sz > cur_lb)
            cur_lb = clique_sz;

        // std::cout << "probe clique = " << clique_sz << " in " << V.size()
        // <<
        // "/"
        //           << size() << std::endl;

        // // std::cout << cur_lb << ": " << max_clique << std::endl;
        //
        // std::cout << cq.cliques[lgst].size() << ":" ;
        //         for (auto u : cq.cliques[lgst]) {
        // 		std::cout << " " << u;
        //         }
        // std::cout << std::endl ;
        //
        // // assert(max_clique.size() == cur_lb);

        if (depth == 0 and cur_lb > lb) {
            lb = cur_lb;
            stats.notify_lb(lb);
            while (degeneracy[df.order[lb_frontier]] < lb)
                ++lb_frontier;
        }

        if (ub > nub) {
            ub = nub;
            stats.notify_ub(ub);
            while (ub_frontier > 0 and degeneracy[df.order[ub_frontier]] >= ub)
                --ub_frontier;
        }

        // cq.clear();

        // if (size() * (size() - 1) == 2 * num_edges)
        //     cur_lb = size();
        // else
        //     cur_lb = 0;

        if (options.verbosity > 4) {
            if (options.verbosity > 5)
                std::cout << std::endl << *this;

            std::cout //<< std::endl << *this
                << nodes.size() << "/" << num_edges << ": [" << lb << ".."
                << cur_lb << ".." << ub << "]" << std::endl;
        }

        //

        if (ub > cur_lb) {
            // pedge = any_non_edge();

            ++depth;

            if (options.verbosity > 3) {
                std::cout << std::setw(3) << std::left << depth << std::right;
                for (int i = 0; i < depth; i++)
                    std::cout << " ";
                std::cout << "/(" << pedge[0] << "," << pedge[1] << ")\n";
            }

            trail.push_back(limit);
            contract(pedge[0], pedge[1]);
            contraction = true;
            // std::cout << *this << std::endl;
        } else if (depth > 0) {

            ++stats.total_conflicts;

            --depth;
            pedge = backtrack(limit);
            addition(pedge[0], pedge[1]);
            contraction = false;

            if (options.verbosity > 3) {
                std::cout << std::setw(3) << std::left << depth << std::right;
                for (int i = 0; i < depth; i++)
                    std::cout << " ";
                std::cout << "+(" << pedge[0] << "," << pedge[1] << ")\n";
            }

            cur_lb = ub - 1;

        } else {
            break;
        }

        // if (stats.total_iteration == 10000)
        //     break;
    }

    stats.custom_force_display(std::cout);
}

void ca_graph::color_preprocess(gc::statistics& stats, gc::options& options)
{

    // exit(1);

    int limit{static_cast<int>(nodes.capacity())};

    check_consistency("beg search");

    bi_graph B;

    gc::bitset V(0, nodes.capacity() - 1, gc::bitset::empt);
    // gc::bitset ub_core(0, nodes.capacity() - 1, gc::bitset::empt);
    // gc::bitset lb_core(0, nodes.capacity() - 1, gc::bitset::empt);
    gc::bitset max_clique(0, nodes.capacity() - 1, gc::bitset::empt);

    degeneracy_finder df(*this);
    dsatur brelaz;
    brelaz.use_recolor = false;
    df.degeneracy_ordering();

    int ub{df.degeneracy + 1};
    // gc::bitset colors(0, ub - 1, gc::bitset::empty);

    brelaz.greedy(*this, rbegin(df.order), rend(df.order), ub);

    std::vector<int> degeneracy(df.degrees);

    std::vector<int> coloring(limit, 0);

    // std::cout << std::setw(3) << *begin(df.order) << " " << std::setw(3) <<
    // df.degrees[*begin(df.order)] << " " << std::setw(3) <<
    // degeneracy[*begin(df.order)] << std::endl;
    //
    for (auto d{begin(df.order)}; d < end(df.order); ++d) {
        auto v{*d};

        coloring[*d] = brelaz.color[*d];

        if (d != begin(df.order))
            degeneracy[v] = std::max(degeneracy[*(d - 1)], degeneracy[v]);
        std::cout << std::setw(3) << v << " " << std::setw(3) << df.degrees[v]
                  << " " << std::setw(3) << degeneracy[v] << std::endl;

        // for()
    }

    // brelaz.greedy(rbegin(df.order), rend(df.order));

    std::vector<int> dg_rank(nodes.capacity());

    int i{0};
    for (auto u : df.order) {
        dg_rank[u] = i++;
    }

    cliquer cq(*this);

    int lb{1 + (num_edges > 0)};
    stats.notify_lb(lb);

    // int ub{df.degeneracy + 1};
    stats.notify_ub(ub);

    int nub{ub};


    // auto largest_degree_criterion = [&](int x, int y) {return
    // matrix[x].size() > matrix[y].size();};
    auto global_criterion = [&](int better, int worse) {

        // 1/ start with the max clique
        auto better_in_maxclique{max_clique.contain(better)};
        auto worse_in_maxclique{max_clique.contain(worse)};
        if (better_in_maxclique and !worse_in_maxclique)
            return true;
        if (worse_in_maxclique != better_in_maxclique)
            return false;

        // 2/ priority on vertices of the maximum ub-core
        // auto better_in_ub_core{dg_rank[better] >= ub_frontier};
        // auto worse_in_ub_core{dg_rank[worse] >= ub_frontier};
        auto better_in_ub_core{degeneracy[better] >= ub};
        auto worse_in_ub_core{degeneracy[worse] >= ub};
        if (better_in_ub_core and !worse_in_ub_core)
            return true;
        if (better_in_ub_core != worse_in_ub_core)
            return false;

        // 3/ priority on vertices of the maximum lb-core
        // auto better_in_lb_core{dg_rank[better] >= lb_frontier};
        // auto worse_in_lb_core{dg_rank[worse] >= lb_frontier};
        auto better_in_lb_core{degeneracy[better] >= lb};
        auto worse_in_lb_core{degeneracy[worse] >= lb};
        if (better_in_lb_core and !worse_in_lb_core)
            return true;
        if (better_in_lb_core != worse_in_lb_core)
            return false;

        // degeneracy [STATIC]
        if (degeneracy[better] > degeneracy[worse])
            return true;
        if (degeneracy[better] < degeneracy[worse])
            return false;

        // degree [DYNAMIC]
        if (matrix[better].size() > matrix[worse].size())
            return true;

        return false;

    };

    // brelaz.full = true;
    // brelaz.clear();

    nub = brelaz.brelaz_color_score(
        *this, ub - 1, global_criterion, size(), 12345);

    // std::vector<int> coloring(capacity(), -1);

    if (nub < ub)
        for (auto v{0}; v < capacity(); ++v) {
            coloring[v] = brelaz.color[v];
            std::cout << " " << coloring[v];
        }
    else
        exit(1);

    // std::cout << " -> " << nub << std::endl;

    brelaz.local_search(*this, coloring, stats, options, begin(brelaz.order),
        end(brelaz.order));

    // brelaz.check_full_consistency(*this, "after local search");

    // for (auto v : nodes) {
    //     coloring[v] = brelaz.color[v];
    //     std::cout << " " << coloring[v];
    // }
    //
    // std::cout << " -> " << nub << std::endl;

    brelaz.restart(*this, coloring);

    // brelaz.check_full_consistency(*this, "after restart");

    // brelaz.full = true;
    // brelaz.findpath_out

    // brelaz.clear();
    //
    //     brelaz.local_search(*this, coloring, stats, options,
    //     begin(brelaz.order),
    //         end(brelaz.order));
    //
    // brelaz.clear();
    //
    //     brelaz.init_local_search(*this, coloring, begin(brelaz.order),
    //         end(brelaz.order));

    brelaz.massacre_vertices(*this, options, stats);
}

void ca_graph::search(gc::statistics& stats, gc::options& options)
{

    // exit(1);

    int limit{static_cast<int>(nodes.capacity())};

    check_consistency("beg search");

    int depth{0};

    /* some structs */
    bi_graph B;
    gc::bitset V(0, nodes.capacity() - 1, gc::bitset::empt);
    gc::bitset max_clique(0, nodes.capacity() - 1, gc::bitset::empt);

    /*** compute degeneracy and k-cores, and initialise ub ***/
    degeneracy_finder df(*this);
    df.degeneracy_ordering();
    int ub{df.degeneracy + 1};

    /*** compute the coloring corresponding to the degeneracy ordering ***/
    dsatur brelaz;
    brelaz.use_recolor = false;
    brelaz.greedy(*this, rbegin(df.order), rend(df.order), ub);
    std::vector<int> degeneracy(df.degrees);
    std::vector<int> coloring(limit, 0);
    for (auto d{begin(df.order)}; d < end(df.order); ++d) {
        auto v{*d};
        coloring[*d] = brelaz.color[*d];
        if (d != begin(df.order))
            degeneracy[v] = std::max(degeneracy[*(d - 1)], degeneracy[v]);
        std::cout << std::setw(3) << v << " " << std::setw(3) << df.degrees[v]
                  << " " << std::setw(3) << degeneracy[v] << std::endl;
    }

    /*** compute the rank in the degeneracy ordering and the frontier between
     * (ir)relevant vertices ***/
    std::vector<int> dg_rank(nodes.capacity());
    int i{0}, lb_frontier{0},
        ub_frontier{static_cast<int>(df.order.size()) - 1};
    for (auto u : df.order) {
        dg_rank[u] = i++;
    }

    int lb{1 + (num_edges > 0)};
    stats.notify_lb(lb);
    stats.notify_ub(ub);

    bool contraction{true};
    arc pedge{0, 0};
    int nub{ub};

    int cur_lb{lb};

    int period = std::pow(10, 6 - options.verbosity);

    // auto largest_degree_criterion = [&](int x, int y) {return
    // matrix[x].size() > matrix[y].size();};
    // auto degeneracy_rank_criterion
    //     = [&](int x, int y) { return dg_rank[x] > dg_rank[y]; };
    auto global_criterion = [&](int better, int worse) {

        // 1/ start with the max clique
        auto better_in_maxclique{max_clique.contain(better)};
        auto worse_in_maxclique{max_clique.contain(worse)};
        if (better_in_maxclique and !worse_in_maxclique)
            return true;
        if (worse_in_maxclique != better_in_maxclique)
            return false;

        // 2/ priority on vertices of the maximum ub-core
        // auto better_in_ub_core{dg_rank[better] >= ub_frontier};
        // auto worse_in_ub_core{dg_rank[worse] >= ub_frontier};
        auto better_in_ub_core{degeneracy[better] >= ub};
        auto worse_in_ub_core{degeneracy[worse] >= ub};
        if (better_in_ub_core and !worse_in_ub_core)
            return true;
        if (better_in_ub_core != worse_in_ub_core)
            return false;

        // 3/ priority on vertices of the maximum lb-core
        // auto better_in_lb_core{dg_rank[better] >= lb_frontier};
        // auto worse_in_lb_core{dg_rank[worse] >= lb_frontier};
        auto better_in_lb_core{degeneracy[better] >= lb};
        auto worse_in_lb_core{degeneracy[worse] >= lb};
        if (better_in_lb_core and !worse_in_lb_core)
            return true;
        if (better_in_lb_core != worse_in_lb_core)
            return false;

        // degeneracy [STATIC]
        if (degeneracy[better] > degeneracy[worse])
            return true;
        if (degeneracy[better] < degeneracy[worse])
            return false;

        // degree [DYNAMIC]
        if (matrix[better].size() > matrix[worse].size())
            return true;

        return false;

    };

    /*** compute a better coloring using dsatur with complex tie breaking ***/
    nub = brelaz.brelaz_color_score(
        *this, ub - 1, global_criterion, size(), 12345);

    if (nub < ub)
        for (auto v{0}; v < capacity(); ++v) {
            coloring[v] = brelaz.color[v];
            std::cout << " " << coloring[v];
        }

    /*** improve the coloring with local search ***/
    brelaz.local_search(*this, coloring, stats, options, begin(brelaz.order),
        end(brelaz.order));

    std::vector<int> largest_cliques;

    cliquer cq(*this);
    while (lb < ub) {

        // check_consistency("search loop");

        if (options.verbosity > 1
            && stats.notify_iteration(depth) % period == 0)
            stats.custom_force_display(std::cout);

        double tbefore = minicsp::cpuTime();
        cq.clear();
        auto clique_sz = cq.find_cliques(
            begin(brelaz.order), end(brelaz.order), options.cliquelimit);

        stats.notify_nclique(cq.num_cliques);
        stats.notify_clique_time(minicsp::cpuTime() - tbefore);

        if (clique_sz > cur_lb) {
            cur_lb = clique_sz;
        }

        tbefore = minicsp::cpuTime();
        largest_cliques.clear();

        for (auto i{0}; i < cq.cliques.size(); ++i) {
            if (cq.cliques[i].size() >= clique_sz-1) {
                largest_cliques.push_back(i);
            }
        }

        max_clique.clear();
        max_clique.union_with(cq.cliques[largest_cliques[0]]);

        int matching_bound = clique_sz;

        if (ub - clique_sz < 4 and largest_cliques.size() > 1) {

            for (auto i{0}; i < largest_cliques.size(); ++i) {
                auto c1{largest_cliques[i]};

                for (auto j{i + 1}; j < largest_cliques.size(); ++j) {
                    auto c2{largest_cliques[j]};

                    B.get_from_cliques(*this, cq.cliques[c1], cq.cliques[c2]);
                    auto mm{B.hopcroftKarp()};

                    auto mm_bound = (cq.cliques[c1].size()
                        + cq.cliques[c2].size() - B.I - mm);

                    if (mm_bound > matching_bound) {
                        matching_bound = mm_bound;
                    }
                }
            }

            stats.notify_bound_delta(clique_sz, matching_bound);
            if (cur_lb < matching_bound) {
                cur_lb = matching_bound;
            }
        }
        stats.notify_matching_time(minicsp::cpuTime() - tbefore);

        tbefore = minicsp::cpuTime();
        brelaz.clear();
        nub = brelaz.brelaz_color_score(
            *this, ub - 1, global_criterion, size(), 12345);
				

        clique_sz = 0;
        for (auto v : brelaz.order) {
            if (brelaz.color[v] < clique_sz) {
                if (brelaz.color[brelaz.order[brelaz.color[v]]]
                    != brelaz.color[v]) {
                    std::cout << "col[" << v << "] = " << brelaz.color[v]
                              << " < " << clique_sz << std::endl;
                }

                assert(brelaz.color[brelaz.order[brelaz.color[v]]]
                    == brelaz.color[v]);
                pedge = arc{v, brelaz.order[brelaz.color[v]]};
                break;
            }
            if (++clique_sz == ub)
                break;
        }
        stats.notify_dsatur_time(minicsp::cpuTime() - tbefore);

        if (clique_sz > cur_lb)
            cur_lb = clique_sz;

        if (depth == 0 and cur_lb > lb) {
            lb = cur_lb;
            stats.notify_lb(lb);
            while (degeneracy[df.order[lb_frontier]] < lb)
                ++lb_frontier;
        }

        if (ub > nub) {
            ub = nub;
            stats.notify_ub(ub);						
            while (ub_frontier > 0 and degeneracy[df.order[ub_frontier]] >= ub)
                --ub_frontier;
        }

        if (options.verbosity > 4) {
            if (options.verbosity > 5)
                std::cout << std::endl << *this;

            std::cout //<< std::endl << *this
                << nodes.size() << "/" << num_edges << ": [" << lb << ".."
                << cur_lb << ".." << ub << "]" << std::endl;
        }

        if (ub > cur_lb) {

            ++depth;

            if (options.verbosity > 3) {
                std::cout << std::setw(3) << std::left << depth << std::right;
                for (int i = 0; i < depth; i++)
                    std::cout << " ";
                std::cout << "/(" << pedge[0] << "," << pedge[1] << ")\n";
            }

            trail.push_back(limit);
            contract(pedge[0], pedge[1]);
            contraction = true;

        } else if (depth > 0) {

            ++stats.total_conflicts;

            --depth;
            pedge = backtrack(limit);
            addition(pedge[0], pedge[1]);
            contraction = false;

            if (options.verbosity > 3) {
                std::cout << std::setw(3) << std::left << depth << std::right;
                for (int i = 0; i < depth; i++)
                    std::cout << " ";
                std::cout << "+(" << pedge[0] << "," << pedge[1] << ")\n";
            }

            cur_lb = ub - 1;

        } else {
            break;
        }

    }

    stats.custom_force_display(std::cout);
}

std::ostream& ca_graph::describe(std::ostream& os, const int verbosity) const
{
    // void gc::ca_graph::describe(std::ostream& os, const int verbosity) const
    // {
    switch (verbosity) {
    case 0:
        os << "n=" << nodes.size() << " m=" << num_edges << std::endl;
        break;

    case 1:
        os << "V=" << nodes << " m=" << num_edges;
        break;

    case 2:
        for (auto v : nodes) {
            os << v << ":";
            for (auto u : matrix[v]) {
                os << " " << u;
            }
            os << std::endl;
        }
        break;

    case 3:
        for (auto v : nodes) {
            os << v << ":";
            for (auto u : matrix[v]) {
                os << " " << u;
            }
            os << std::endl;
        }
        for (auto t : trail)
            std::cout << " " << t;
        std::cout << std::endl;
    }

    return os;
}

void ca_graph::check_consistency(const char* msg)
{
    auto count{2 * num_edges};
    std::vector<int> f(nodes.capacity(), 0);
    std::vector<int> b(nodes.capacity(), 0);

    if (nodeset.size() != nodes.size()) {
        std::cout << msg << ": |nodes|=" << nodes.size()
                  << " and |nodeset|=" << nodeset.size() << std::endl;
        exit(1);
    }
    // assert(nodeset.size() == nodes.size());

    int total_rank{0};
    for (auto u : nodes) {
        if (!nodeset.contain(u)) {
            std::cout << msg << ": " << u << " in nodes but not in nodeset\n";
            exit(1);
        }
        // assert(nodeset.contain(u));
        total_rank += rank[u];

        auto prev{-1};
        for (auto v : matrix[u]) {

            if (!nodeset.contain(u)) {
                std::cout << msg << ": " << v << " in N(" << u
                          << ") but not in nodes\n"
                          << std::endl;
                exit(1);
            }

            if (v <= prev) {
                std::cout << msg << ": nodes not sorted N(" << u << ") =";
                for (auto x : matrix[u]) {
                    std::cout << " " << x;
                    if (x == v) {
                        std::cout << "...\n";
                        break;
                    }
                }
                exit(1);
            }

            // assert(nodes.contain(v));
            ++f[u];
            ++b[v];
            --count;
        }
    }

    assert(total_rank == nodes.capacity());

    if (count)
        std::cout << count << std::endl;

    assert(count == 0);

    for (auto u : nodes) {

        if (f[u] != b[u]) {
            std::cout << msg << ": " << u << " has " << f[u]
                      << " neighbors but appears in " << b[u]
                      << " neighborhoods\n"
                      << std::endl;
            exit(1);
        }
        // assert(f[u] == b[u]);
    }
}

// Returns size of maximum matching
int bi_graph::hopcroftKarp()
{
    // Initialize result
    int result = 0;

    // Keep updating the result while there is an
    // augmenting path.
    while (bfs()) {

        // std::cout << "\nmatching (" << result << "):\n";
        //
        // for (int u = 0; u < N; ++u) {
        //     if (matching[u] != T)
        //         std::cout << u << " -- " << matching[u] << std::endl;
        // }
        //
        // std::cout << "BFS\n";

        // Find a free vertex
        for (int u = N; u < T; ++u)
            // If current vertex is free and there is
            // an augmenting path from current vertex
            if (matching[u] == T && dfs(u))
                result++;

        // for (int u = 0; u < T; ++u) {
        //     std::cout << std::setw(2) << u << " " << std::setw(3) <<
        //     original[u]
        //               << " ";
        //     if (dist[u] != INFTY)
        //         std::cout << dist[u] << std::endl;
        //     else
        //         std::cout << "inf" << std::endl;
        // }
        // if (dist[T] != INFTY)
        //     std::cout << std::setw(2) << T << "     " << dist[T] <<
        //     std::endl;
        // else
        //     std::cout << std::setw(2) << T << "     inf\n";
        // std::cout << std::endl;
    }
    return result;
}

// Returns true if there is an augmenting path, else returns
// false
bool bi_graph::bfs()
{

    // First layer of vertices (set distance as 0)
    for (int u = N; u < T; ++u) {
        // If this is a free vertex, add it to queue
        if (matching[u] == T) {
            // u is not matched
            dist[u] = 0;
            Q.push_back(u);
            // std::cout << " " << u;
        } else {
            dist[u] = INFTY;
            // std::cout << " -" << u;
        }
    }

    // std::cout << std::endl;

    // Initialize distance to T as infinite
    dist[T] = INFTY;

    // Q is going to contain vertices of left side only.
    while (q < Q.size()) {
        // Dequeue a vertex
        int u = Q[q++];

        // std::cout << u << ":";

        // If this node is not T and can provide a shorter path to T
        if (dist[u] < dist[T]) {
            // Get all adjacent vertices of the dequeued vertex u
            for (auto v : matrix[u]) {
                // If pair of v is not considered so far
                // (v, pairV[V]) is not yet explored edge.
                if (dist[matching[v]] == INFTY) {
                    // Consider the pair and add it to queue
                    dist[matching[v]] = dist[u] + 1;
                    Q.push_back(matching[v]);

                    // std::cout << " " << matching[v];
                }
            }
        }

        // std::cout << std::endl;
    }

    // If we could come back to T using alternating path of distinct
    // vertices then there is an augmenting path
    return (dist[T] != INFTY);
}

// Returns true if there is an augmenting path beginning with free vertex u
bool bi_graph::dfs(const int u)
{
    if (u != T) {
        for (auto v : matrix[u]) {
            if (dist[matching[v]] == dist[u] + 1) {
                if (dfs(matching[v])) {
                    matching[u] = v;
                    matching[v] = u;
                    return true;
                }
            }
        }

        dist[u] = INFTY;
        return false;
    }

    return true;
}

std::ostream& bi_graph::describe(std::ostream& os) const
{

    for (int i = 0; i < T; ++i) {
        if (i == N)
            std::cout << std::endl;

        os << original[i] << ":";
        for (auto j : matrix[i]) {
            os << " " << original[j];
        }
        os << std::endl;
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const gc::ca_graph& g)
{
    return g.describe(os, 3);
}

std::ostream& operator<<(std::ostream& os, const gc::bi_graph& g)
{
    return g.describe(os);
}
}