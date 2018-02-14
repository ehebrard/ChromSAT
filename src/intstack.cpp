#include "constants.hpp"
#include "intstack.hpp"
#include "bitset.hpp"
#include <cstring>
#include <cassert>

void IntStack::initialise(IntStack& shared, const int sz)
{
    index_capacity = shared.index_capacity;

    list_capacity = sz;
    list_ = new int[list_capacity];

    start_ = shared.start_;
    index_ = shared.index_;

    size = 0;
}

void IntStack::initialise(
    const int lb, const int ub, const int sz, const bool full)
{

    index_capacity = ub - lb + 1;

    list_capacity = sz;
    list_ = new int[list_capacity];
    start_ = new unsigned int[index_capacity];
    index_ = start_ - lb;

    for (int i = lb; i <= ub; ++i) {
        index_[i] = i - lb;
        if (i - lb < (int)list_capacity)
            list_[i - lb] = i;
    }

    size = (full ? sz : 0);
}

std::ostream& IntStack::display(std::ostream& os) const
{
    int min = INFTY;
    int max = -INFTY;

    if (size) {
        for (unsigned int i = 0; i < size; ++i) {
            if (list_[i] < min)
                min = list_[i];
            if (list_[i] > max)
                max = list_[i];
        }
    } else {
        min = max = 0;
    }

    BitSet elts(min, max, BitSet::empt);
    for (unsigned int i = 0; i < size; ++i) {
        elts.add(list_[i]);
    }

    os << elts;

    if (elts.size() != size) {
        std::cout << "ERROR " << elts.size() << " / " << size << std::endl;
        assert(0);
    }

    // os << "(";
    // bool not_empty = (size>0);

    // if(not_empty) os << list_[0];
    // for(unsigned int i=1; i<size; ++i)
    //  os << " " << list_[i];
    // os << ")";
    return os;
}

std::string IntStack::to_str() const
{
    int min = INFTY;
    int max = -INFTY;

    if (size) {
        for (unsigned int i = 0; i < size; ++i) {
            if (list_[i] < min)
                min = list_[i];
            if (list_[i] > max)
                max = list_[i];
        }
    } else {
        min = max = 0;
    }

    BitSet elts(min, max, BitSet::empt);
    for (unsigned int i = 0; i < size; ++i) {
        elts.add(list_[i]);
    }

    return elts.to_str();
}

std::ostream& IntStack::display(std::ostream& os, const int rank) const
{
    int min = 0;
    int max = index_capacity - 1;

    BitSet elts(min, max, BitSet::empt);
    for (unsigned int i = 0; i < size; ++i) {
        elts.add(list_[i]);
    }

    os << elts;

    elts.clear();

    for (int i = size; i < rank; ++i) {
        elts.add(list_[i]);
    }

    os << elts;

    elts.clear();

    for (unsigned int i = rank; i < list_capacity; ++i) {
        elts.add(list_[i]);
    }

    os << elts;

    return os;
}

void IntStack::extend_list()
{
    unsigned int increment = ((list_capacity + 1) << 1);
    list_capacity += increment;

    int* new_list = new int[list_capacity];
    memcpy(new_list, list_, (list_capacity - increment) * sizeof(int));

    // for(unsigned int i=0; i<list_capacity-increment; ++i)
    //          new_list[i] = list_[i];

    delete[] list_;
    list_ = new_list;
}

void IntStack::extend(const int new_elt)
{
    int lb = (int)(start_ - index_), new_lb = lb;
    int ub = index_capacity + lb - 1, new_ub = ub;
    if (new_elt < lb) {
        new_lb = new_elt;
    } else if (new_elt > ub) {
        new_ub = new_elt;
    } else {
        return;
    }

    unsigned int new_index_capacity = new_ub - new_lb + 1;
    if (new_index_capacity < index_capacity * 2)
        new_index_capacity = index_capacity * 2;
    if (new_lb < lb) {
        new_lb = ub - new_index_capacity + 1;
    } else {
        new_ub = lb + new_index_capacity - 1;
    }

    unsigned int* aux_start = start_;
    start_ = new unsigned int[new_index_capacity];
    std::fill(start_ + index_capacity, start_ + new_index_capacity, INFTY);
    memcpy(start_ + (lb - new_lb), aux_start,
        index_capacity * sizeof(unsigned int));
    delete[] aux_start;

    index_ = start_ - new_lb;
    int k = 0;
    for (int i = new_lb; i < lb; ++i) {
        index_[i] = size + k++;
        // list_[index_capacity+k++] = i;
    }
    for (int i = ub + 1; i <= new_ub; ++i) {
        index_[i] = size + k++;
        // list_[index_capacity+k++] = i;
    }

    index_capacity = new_index_capacity;
}
//@}

/*!@name Accessors*/
//@{
int IntStack::get_min() const
{
    int the_min = INFTY;
    int offset = (index_ - start_);
    unsigned int explored;
    int stop_crit = size + offset;

    for (explored = 0; explored < size
         //&& size - explored <= (the_min - offset);
         && the_min + (int)explored >= stop_crit;
         ++explored) {
        if (list_[explored] < the_min)
            the_min = list_[explored];
    }
    if (explored < size) {
        while (offset < the_min && index_[offset] >= size)
            ++offset;
        the_min = offset;
    }

    // int the_min = INFTY;
    // if(size) {
    //  the_min = list_[0];
    //  // ratio size / index_capacity
    //  if(4*size < index_capacity) {
    //    for(unsigned int i=1; i<size; ++i)
    //      if(list_[i] < the_min) the_min = list_[i];
    //  } else {
    //    int val=(index_ - start_);
    //    while( val<the_min && index_[val] >= size )
    //      ++val;
    //    the_min = val;
    //  }
    // }

    return the_min;
}

int IntStack::get_max() const
{
    int the_max = -INFTY;
    if (size) {
        the_max = list_[0];
        // ratio size / index_capacity
        if (4 * size < index_capacity) {
            for (unsigned int i = 1; i < size; ++i)
                if (list_[i] > the_max)
                    the_max = list_[i];
        } else {
            int val = (index_capacity + index_ - start_);
            while (val > the_max && index_[val] >= size)
                --val;
            the_max = val;
        }
    }
    return the_max;
}

bool IntStack::safe_contain(const int elt) const
{

    // std::cout << "does " << this << " contain " << elt << " (absolute index:
    // "
    //          << ((unsigned)(elt-(int)(start_-index_))) << "/" <<
    //          index_capacity
    // ;
    // if((unsigned)(elt-(int)(start_-index_))<index_capacity) {
    //  std::cout << "list index: " << index_[elt] ;
    // }
    // std::cout << ") => " <<
    // (((unsigned)(elt-(int)(start_-index_))<index_capacity &&
    // index_[elt]<size) ? "yep" : "no") << "\n";

    return ((unsigned)(elt - (int)(start_ - index_)) < index_capacity
        && index_[elt] < size);
}

bool IntStack::contain(const int elt) const { return index_[elt] < size; }

bool IntStack::contain(
    const int elt, const int min_idx, const int max_idx) const
{
    int idx = index_[elt];
    return idx < max_idx && idx >= min_idx;
}

bool IntStack::empty() const { return !size; }

int IntStack::next(const int elt) const
{
    unsigned int idx = index_[elt] + 1;
    return (idx < size ? list_[idx] : elt);
}
int IntStack::prev(const int elt) const
{
    int idx = index_[elt] - 1;
    return (idx >= 0 ? list_[idx] : elt);
}

int IntStack::operator[](const unsigned int idx) const { return list_[idx]; }

int& IntStack::operator[](const unsigned int idx) { return list_[idx]; }
//@}

/*!@name List Manipulation*/
//@{
int* IntStack::begin() { return list_; }
const int* IntStack::begin() const { return list_; }

int* IntStack::end() { return &(list_[size]); }
const int* IntStack::end() const { return &(list_[size]); }

int* IntStack::end_mem() { return list_ + list_capacity; }

void IntStack::fill() { size = list_capacity; }

void IntStack::clear() { size = 0; }

void IntStack::set_to(const int elt)
{
    size = 1;
    index_[*list_] = index_[elt];
    list_[index_[elt]] = *list_;
    *list_ = elt;
    index_[elt] = 0;
}

void IntStack::remove(const int elt)
{
    --size;
    index_[list_[size]] = index_[elt];
    list_[index_[elt]] = list_[size];
    list_[size] = elt;
    index_[elt] = size;
}

void IntStack::remove(const int elt, int& rank)
{
    index_[list_[rank]] = index_[elt];
    list_[index_[elt]] = list_[rank];
    list_[rank] = elt;
    index_[elt] = rank;
    ++rank;
}

void IntStack::move(const int elt, const int idx)
{
    int o_elt = list_[idx];
    int s_idx = index_[elt];

    list_[s_idx] = o_elt;
    index_[o_elt] = s_idx;

    index_[elt] = idx;
    list_[idx] = elt;
}

int IntStack::next() { return list_[size]; }

int IntStack::pop() { return list_[--size]; }

int IntStack::pop_head()
{
    --size;
    index_[list_[size]] = 0;
    const int elt = *list_;
    *list_ = list_[size];
    list_[size] = elt;
    index_[elt] = size;
    return elt;
}

int IntStack::head() const { return *list_; }

int IntStack::back() const { return list_[size - 1]; }

void IntStack::init_add(const int elt)
{
    if (size == list_capacity)
        extend_list();
    if ((unsigned)(elt - (int)(start_ - index_)) >= index_capacity) {
        extend(elt);
    }
    list_[size] = elt;
    index_[elt] = size;
    ++size;
}

void IntStack::add(const int elt)
{
    index_[list_[size]] = index_[elt];
    list_[index_[elt]] = list_[size];
    list_[size] = elt;
    index_[elt] = size;
    ++size;
}

void IntStack::add(const int elt, int& rank)
{
    --rank;
    index_[list_[rank]] = index_[elt];
    list_[index_[elt]] = list_[rank];
    list_[rank] = elt;
    index_[elt] = rank;
}

void IntStack::set(const int elt, const int idx)
{
    index_[list_[idx]] = index_[elt];
    list_[index_[elt]] = list_[idx];
    list_[idx] = elt;
    index_[elt] = idx;
}

void IntStack::safe_add(const int elt)
{
    if (!safe_contain(elt))
        extend(elt);
    add(elt);
}

// create a new element that can potentially be outside the bounds
void IntStack::create(const int elt)
{
    extend(elt);
    add(elt);
}

void IntStack::ordered_add(const int elt)
{
    // the first non-element goes where elt was
    index_[list_[size]] = index_[elt];
    list_[index_[elt]] = list_[size];

    int idx = size;
    while (idx
        && list_[idx - 1]
            > elt) { // push every values greater than elt above elt
        list_[idx] = list_[idx - 1];
        index_[list_[idx - 1]] = idx;
        --idx;
    }

    list_[idx] = elt;
    index_[elt] = idx;
    ++size;
}

void IntStack::revert_to(const int level) { size = level; }

void IntStack::index()
{
    for (unsigned int i = 0; i < list_capacity; ++i)
        index_[list_[i]] = i;
}

int IntStack::get_index(const int elt) { return index_[elt]; }
//@}

std::ostream& operator<<(std::ostream& os, const IntStack& x)
{
    return x.display(os);
}

std::ostream& operator<<(std::ostream& os, const IntStack* x)
{
    return (x ? x->display(os) : os);
}
