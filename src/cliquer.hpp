#ifndef __CG_CLIQUER_HH
#define __CG_CLIQUER_HH

#include <vector>

#include "bitset.hpp"


namespace gc
{


template <class graph_struct> struct cliquer {
    const graph_struct& g;
    std::vector<std::vector<int>> cliques;
    std::vector<gc::bitset> candidates;
    std::vector<int> last_clique;
    size_t num_cliques;
    size_t maxcliquesize;

    cliquer(const graph_struct& ig)
        : g(ig)
        , num_cliques(0)
        , maxcliquesize(0)
    {
        last_clique.resize(g.capacity());
    }

    // clear previously cached results
    void clear()
    {
        num_cliques = 0;
        maxcliquesize = 0;
    }
    // initialize a new clique
    void new_clique();

    // insert v into the clq^th clique. assumes it fits
    void insert(int v, int clq);

    // heuristically find a set of cliques and return the size of the
    // largest
    template <class ForwardIterator>
    int find_cliques(ForwardIterator beg_ordering, ForwardIterator end_ordering, const int limit=0xfffffff)
    {
        clear();
        for (auto it{beg_ordering}; it != end_ordering; ++it) {

            // std::cout << *it ;

            bool found{false};
            for (int i = 0; i != num_cliques; ++i)
                if (candidates[i].fast_contain(*it)) {
                    found = true;
                    insert(*it, i);
                    // std::cout << ":" << i << " ";
                }
            if (!found && num_cliques < limit) {
                new_clique();
                insert(*it, num_cliques - 1);
                // std::cout << ":" << (num_cliques - 1) << "* ";
            }
        }

        // std::cout << "| ";

        for (auto it{beg_ordering}; it != end_ordering; ++it) {
            for (int i = last_clique[*it] + 1; i < num_cliques; ++i)
                if (candidates[i].fast_contain(*it)) {
                    insert(*it, i);
                    // std::cout << *it << ":" << (num_cliques - 1) << " ";
                }
        }

        // std::cout << " -> " << maxcliquesize << std::endl;

        return maxcliquesize;
    }

    // heuristically find a set of cliques and return the size of the
    // largest
    template <typename ForwardIterator, typename set>
    int find_cliques_in(
        ForwardIterator beg_ordering, ForwardIterator end_ordering, set nodes)
    {
        clear();
        for (auto it{beg_ordering}; it != end_ordering; ++it) {
            // std::cout << *it ;

            if (!nodes.contain(*it))
                continue;

            bool found{false};
            for (int i = 0; i != num_cliques; ++i)
                if (candidates[i].fast_contain(*it)) {
                    found = true;
                    insert(*it, i);
                    // std::cout << ":" << i << " ";
                }
            if (!found) {
                new_clique();
                insert(*it, num_cliques - 1);
                // std::cout << ":" << (num_cliques - 1) << "* ";
            }
        }

        // std::cout << "| ";

        for (auto it{beg_ordering}; it != end_ordering; ++it) {

            if (!nodes.contain(*it))
                continue;

            for (int i = last_clique[*it] + 1; i < num_cliques; ++i)
                if (candidates[i].fast_contain(*it)) {
                    insert(*it, i);
                    // std::cout << *it << ":" << (num_cliques - 1) << " ";
                }
        }

        // std::cout << " -> " << maxcliquesize << std::endl;

        return maxcliquesize;
    }
};

template <class graph_struct> void cliquer<graph_struct>::new_clique()
{
    if (cliques.size() == num_cliques) {
        cliques.resize(num_cliques + 1);
        candidates.resize(num_cliques + 1);
        candidates.back().initialise(0, g.capacity(), bitset::empt);
    }

    cliques[num_cliques].clear();
    candidates[num_cliques].copy(g.nodeset);
    ++num_cliques;
}

// insert v into the clq^th clique. assumes it fits
template <class graph_struct> void cliquer<graph_struct>::insert(int v, int clq)
{
    cliques[clq].push_back(v);
    candidates[clq].intersect_with(begin(g.matrix[v]), end(g.matrix[v]));
    last_clique[v] = clq;
    maxcliquesize = std::max(maxcliquesize, cliques[clq].size());
}

template <class InputIterator1, class InputIterator2, class OutputIterator,
    class DeltaIterator, class InterIterator>
OutputIterator set_union_delta(InputIterator1 first1, InputIterator1 last1,
    InputIterator2 first2, InputIterator2 last2, OutputIterator result_union,
    DeltaIterator delta, // 2 \ 1
    InterIterator inter)
{
    while (first1 != last1 and first2 != last2) {
        if (*first1 < *first2) {
            *result_union = *first1;
            ++result_union;
            ++first1;
        } else {
            *result_union = *first2;
            ++result_union;
            if (*first2 < *first1) {
                *delta = *first2;
                ++first2;
                ++delta;
            } else {
                *inter = *first1;
                ++inter;
                ++first1;
                ++first2;
            }
        }
    }

    std::copy(first2, last2, delta);
    std::copy(first2, last2, result_union);

    return std::copy(first1, last1, result_union);
}

template <class ForwardIt, class T>
ForwardIt find_in(ForwardIt first, ForwardIt last, const T& value)
{
    ForwardIt it;
    typename std::iterator_traits<ForwardIt>::difference_type count, step;
    count = std::distance(first, last);

    while (count > 0) {
        it = first;
        step = count / 2;
        std::advance(it, step);
        if (*it < value) {
            first = ++it;
            count -= step + 1;
        } else if (*it == value) {
            return it;
        } else
            count = step;
    }

    return first;
}

} // namespace gc

#endif
