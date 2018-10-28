
#include <assert.h>

#include "sparseset.hpp"

// sparseset::sparseset(std::vector<size_t>& i, const size_t n)
// 	: index_(i)
// {
//     size_ = 0;
// 		index_.resize(n, NOVAL);
// }
//
//
// sparseset::sparseset(sparseset& s)
// 	: index_(s.index_)
// {
//     size_ = 0;
// }

sparseset::sparseset() { size_ = 0; }

void sparseset::binds(std::vector<size_t>* i) { index_ = i; }

void sparseset::save(size_t& stamp) { stamp = size_; }
void sparseset::restore(const size_t stamp) { size_ = stamp; }

//@}

/*!@name Accessors*/
//@{

bool sparseset::safe_contain(const int elt) const
{
    if (elt >= 0 && (size_t)elt < index_->size())
        return contain(elt);
    return false;
}

bool sparseset::contain(const int elt) const { return index(elt) < size_; }

size_t sparseset::size() const { return size_; }

bool sparseset::empty() const { return size_ == 0; }

int sparseset::next(const int elt) const
{
    size_t idx = index(elt) + 1;
    return (idx < size_ ? list_[idx] : elt);
}
int sparseset::prev(const int elt) const
{
    size_t idx = index(elt);
    return (idx > 0 ? list_[idx - 1] : elt);
}

int sparseset::operator[](const size_t idx) const { return list_[idx]; }

int& sparseset::operator[](const size_t idx) { return list_[idx]; }
//@}

/*!@name List Manipulation*/
//@{
std::vector<int>::iterator sparseset::begin() { return list_.begin(); }
std::vector<int>::reverse_iterator sparseset::rbegin()
{
    return list_.rend() + size_;
}

std::vector<int>::iterator sparseset::end() { return list_.begin() + size_; }
std::vector<int>::reverse_iterator sparseset::rend() { return list_.rend(); }

std::vector<int>::const_iterator sparseset::begin() const
{
    return list_.begin();
}
std::vector<int>::const_reverse_iterator sparseset::rbegin() const
{
    return list_.rend() + size_;
}

std::vector<int>::const_iterator sparseset::end() const
{
    return list_.begin() + size_;
}
std::vector<int>::const_reverse_iterator sparseset::rend() const
{
    return list_.rend();
}

void sparseset::fill() { size_ = list_.size(); }

void sparseset::clear() { size_ = 0; }

void sparseset::remove(const int elt)
{
    auto last = list_[--size_];
    (*index_)[last] = index(elt);
    list_[index(elt)] = last;
    list_[size_] = elt;
    (*index_)[elt] = size_;
}

void sparseset::move_up(const int elt, const int idx_to)
{
    auto idx_from = index(elt);

    assert(index(elt) <= static_cast<size_t>(idx_to));

    auto last = list_[idx_to];
    (*index_)[last] = idx_from;
    list_[idx_from] = last;
    list_[idx_to] = elt;
    (*index_)[elt] = idx_to;
}

void sparseset::pop_back() { --size_; }
void sparseset::pop_head() { remove(list_[0]); }

int sparseset::head() const { return list_[0]; }
int sparseset::back() const { return list_[size_ - 1]; }

void sparseset::safe_add(const int elt)
{
    if (elt >= 0) {
        if (static_cast<size_t>(elt) >= index_->size()) {
            index_->resize(elt + 1);
        }
        add(elt);
    }
}

void sparseset::add(const int elt)
{
    if (index(elt) >= size_) {
        if (index(elt) == NOINDEX) {
            (*index_)[elt] = list_.size();
            list_.push_back(elt);
        }
        push(elt);
    }
}

void sparseset::push(const int elt)
{
    auto next = list_[size_];
    (*index_)[next] = index(elt);
    list_[index(elt)] = next;
    (*index_)[elt] = size_;
    list_[size_++] = elt;
}

int sparseset::index(const int elt) const { return *(index_->begin() + elt); }
//@}

std::ostream& sparseset::display(std::ostream& os) const
{
    os << "(";
    for (auto it = begin(); it < end(); ++it) {
        os << " " << *it;
    }
    os << " )";
    return os;
}

std::ostream& operator<<(std::ostream& os, const sparseset& x)
{
    return x.display(os);
}

std::ostream& operator<<(std::ostream& os, const sparseset* x)
{
    return (x ? x->display(os) : os);
}
