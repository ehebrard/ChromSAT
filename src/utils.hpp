#ifndef __GC_UTILS_HPP
#define __GC_UTILS_HPP

#include <iostream>

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
} // namespace gc

#endif
