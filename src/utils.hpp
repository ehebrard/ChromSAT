#ifndef __GC_UTILS_HPP
#define __GC_UTILS_HPP

#include "minicsp/mtl/Vec.h"
#include <iostream>

struct MarsagliaXORSHF {

    static unsigned long x;
    static unsigned long y;
    static unsigned long z;

    MarsagliaXORSHF() { seed(123456789); }

    void seed(unsigned long s)
    {
        x = s;
        y = 362436069;
        z = 521288629;
    }

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
};

namespace gc
{
template <typename Cont> struct print_container {
    const Cont& cont;
    print_container(const Cont& cont)
        : cont(cont)
    {
    }
};

template <typename Cont>
std::ostream& operator<<(std::ostream& os, print_container<Cont> p)
{
    os << "[";
    for (auto&& e : p.cont)
        os << e << " ";
    os << "]";
    return os;
}

// //--------------------------------------------------
// // erase elements from a vector/vec, subject to a Pred. Returns true
// // if it removed any elements
// template <typename T, typename A, typename Pred>
// bool erase_if(std::vector<T, A>& v, Pred p)
// {
//     auto newend = remove_if(begin(v), end(v), p);
//     if (newend == end(v))
//         return false;
//     v.erase(newend, end(v));
//     return true;
// }
//
// template <typename T, typename Pred> bool erase_if(vec<T>& v, Pred p)
// {
//     int i{0}, j{0};
//     for (; i != v.size(); ++i)
//         if (p(v[i]))
//             v[j++] = v[i];
//     v.shrink(i - j);
//     return i != j;
// }

} // namespace gc

#endif
