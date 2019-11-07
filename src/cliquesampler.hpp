#include <iostream>

#include <algorithm>

#include "graph.hpp"

#ifndef __CLIQUESAMPLER_HPP
#define __CLIQUESAMPLER_HPP

// #define _DEBUG_SAMPLE

namespace gc
{

static unsigned long x = 123456789, y = 362436069, z = 521288629;

unsigned long xorshf96(void)
{ // period 2^96-1
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;
}

template <typename T> class no_weight
{

private:
    T val;

public:
    no_weight(const T v)
        : val(v)
    {
    }

    inline T operator[](const int x) const { return val; }
};

template <typename T> struct clique_sampler {

    // static const no_weight<T> default_w{1};

    std::vector<int> start_set;
    std::vector<int> cand_set;
    std::vector<int> probed;
    std::vector<std::vector<int>> domain;
    std::vector<T> maxes;
    std::vector<int> buffer;

    gc::bitset nodeset;

    std::vector<int> clique;
    std::vector<int> best;
    T clique_weight;
    T max_weight;
    T lb;
    size_t probewidth;

    void set_seed(unsigned long s) { x = s; }

    clique_sampler() {}

    template <class viterator>
    clique_sampler(viterator first, viterator last, const int n)
    {
        set_domain(first, last, n, true);
    }

    template <class viterator>
    void set_domain(
        viterator first, viterator last, const int n, const bool full)
    {

#ifdef _DEBUG_SAMPLE
        std::cout << " sample domain:";
#endif

        nodeset.reinitialise(0, n - 1, 0);
        start_set.clear();

        if (full)
            nodeset.fill();

        for (auto vp{first}; vp != last; ++vp) {
            start_set.push_back(*vp);
            nodeset.fast_add(*vp);

#ifdef _DEBUG_SAMPLE
            std::cout << " " << (*vp);
#endif
        }

#ifdef _DEBUG_SAMPLE
        std::cout << std::endl;
#endif
    }

    template <class graph_struct, class viterator, class map_struct>
    T find_clique(graph_struct& g, const T l, viterator first, viterator last,
        const size_t basewidth, const size_t pw, const map_struct& weight)
    {

#ifdef _DEBUG_SAMPLE
        if (first != last)
            std::cout << " start set (" << start_set.size() << "):";
#endif

        lb = l;
        probewidth = pw;

        domain.resize(probewidth);
        maxes.resize(probewidth);

        // collect possible starting points
        for (auto vp{first}; vp != last; ++vp) {
            start_set.push_back(*vp);

#ifdef _DEBUG_SAMPLE
            std::cout << " " << (*vp);
#endif
        }

#ifdef _DEBUG_SAMPLE
        if (first != last)
            std::cout << std::endl;
#endif

        auto limit{std::min(basewidth, start_set.size())};

        // shuffle
        for (auto i{0}; i < limit; ++i) {
            std::swap(start_set[i],
                start_set[i + (xorshf96() % (start_set.size() - i))]);
        }

        start_set.resize(limit);

        for (auto v : start_set) {
            
						clique.clear();
            clique.push_back(v);
            clique_weight = weight[v];

            cand_set.clear();
            for (auto u : g.matrix[v])
                if (nodeset.fast_contain(u)) {
                    cand_set.push_back(u);
                    max_weight += weight[u];
                }
            probed = cand_set;

            probe(g, weight);

#ifdef _DEBUG_SAMPLE
            std::cout << " find clique: "
                      << print_container<std::vector<int>>(clique) << " ("
                      << clique_weight << ")\n";
#endif

#ifdef _DEBUG_SAMPLE
            T total{0};
            for (auto x : clique)
                total += weight[x];

            assert(clique_weight == total);
            std::vector<int> Nv;
            std::sort(clique.begin(), clique.end());
            for (auto v : clique) {
                std::set_intersection(g.matrix[v].begin(), g.matrix[v].end(),
                    clique.begin(), clique.end(), std::back_inserter(Nv));
                if (Nv.size() != clique.size() - 1) {
                    std::cout << "N(" << v << ") = ";
                    for (auto u : g.matrix[v])
                        std::cout << " " << u;
                    std::cout << "\nclique = ";
                    for (auto u : clique)
                        std::cout << " " << u;
                    std::cout << "\n";
                    std::cout << "intersection = ";
                    for (auto u : Nv)
                        std::cout << " " << u;
                    std::cout << "\n";
                    exit(1);
                }
                Nv.clear();
            }
#endif

            // std::cout << " " << p ;
            // std::cout.flush();
            // lb = std::max(lb, clique_weight);
            if (clique_weight > lb) {
                lb = clique_weight;
                best = clique;
            }
        }
				// std::cout << std::endl;

        return lb;
    }

    template <class graph_struct, class map_struct>
    T probe(graph_struct& g, const map_struct& weight)
    {

#ifdef _DEBUG_SAMPLE
        std::cout << "probe (" << lb
                  << "): " << print_container<std::vector<int>>(clique) << " ("
                  << clique_weight << ") + "
                  << print_container<std::vector<int>>(cand_set) << " ("
                  << max_weight << ")" << std::endl;
#endif

        // if (clique.size() + cand_set.size() <= lb or cand_set.size() == 0)
        //     return clique.size();

        if (clique_weight + max_weight <= lb or cand_set.size() == 0)
            return clique_weight;

        auto w{std::min(probewidth, cand_set.size())};

        int largest{0};
        // select 'probewidth' samples
        for (int i = 0; i < w; ++i) {
            std::swap(
                probed[i], probed[i + (xorshf96() % (probed.size() - i))]);

            auto v{probed[i]};
            domain[i].clear();
            for (auto u : g.matrix[v])
                if (nodeset.fast_contain(u))
                    domain[i].push_back(u);

            // Put intersection in the buffer
            buffer.clear();
            std::set_intersection(domain[i].begin(), domain[i].end(),
                cand_set.begin(), cand_set.end(), std::back_inserter(buffer));

            // Copy the buffer in domain[i]
            std::swap(buffer, domain[i]);

            maxes[i] = 0;
            for (auto u : domain[i])
                maxes[i] += weight[u];

#ifdef _DEBUG_SAMPLE
            std::cout << v << " -> " << maxes[i] << std::endl;
#endif

            if (maxes[largest] < maxes[i])
                largest = i;
            // if (domain[largest].size() < domain[i].size())
            //     largest = i;
        }

        clique.push_back(probed[largest]);
        clique_weight += weight[probed[largest]];
        max_weight = maxes[largest];

        std::swap(cand_set, domain[largest]);
        probed = cand_set;

        return probe(g, weight);
    }
};

} // namespace gc

#endif // __CLIQUESAMPLER_HPP
