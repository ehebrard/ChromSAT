
#include <assert.h>

#include "intstack.hpp"


IntStack::IntStack(const size_t n)
{
		size_ = 0;
    reserve(n);
}

void IntStack::reserve(const size_t n)
{
	while(list_.size() < n) {
		index_.push_back(list_.size());
		list_.push_back(list_.size());
	}
}

//@}

/*!@name Accessors*/
//@{

bool IntStack::safe_contain(const int elt) const
{
	if(elt >= 0 && (size_t)elt < index_.size())
		return contain(elt);
	return false;
}

bool IntStack::contain(const int elt) const { return index_[elt] < size_; }

size_t IntStack::size() const { return size_; }

bool IntStack::empty() const { return size_ == 0; }

int IntStack::next(const int elt) const
{
    size_t idx = index_[elt] + 1;
    return (idx < size_ ? list_[idx] : elt);
}
int IntStack::prev(const int elt) const
{
    size_t idx = index_[elt];
    return (idx > 0 ? list_[idx-1] : elt);
}

int IntStack::operator[](const size_t idx) const { return list_[idx]; }

int& IntStack::operator[](const size_t idx) { return list_[idx]; }
//@}

/*!@name List Manipulation*/
//@{
std::vector< int >::iterator IntStack::begin() { return list_.begin(); }
std::vector< int >::reverse_iterator IntStack::rbegin() { return list_.rend()+size_; }

std::vector< int >::iterator IntStack::end() { return list_.begin()+size_; }
std::vector< int >::reverse_iterator IntStack::rend() { return list_.rend(); }

std::vector< int >::const_iterator IntStack::begin() const { return list_.begin(); }
std::vector< int >::const_reverse_iterator IntStack::rbegin() const { return list_.rend()+size_; }

std::vector< int >::const_iterator IntStack::end() const { return list_.begin()+size_; }
std::vector< int >::const_reverse_iterator IntStack::rend() const { return list_.rend(); }


void IntStack::fill() { size_ = list_.size(); }

void IntStack::clear() { size_ = 0; }

void IntStack::remove(const int elt)
{
		auto last = list_[--size_];
		index_[last] = index_[elt];
    list_[index_[elt]] = last;
    list_[size_] = elt;
    index_[elt] = size_;
}

void IntStack::move_up(const int elt, const int idx_to) {
		auto idx_from =  index_[elt];
		
		assert( index_[elt] <= idx_to );
		
		auto last = list_[idx_to];
		index_[last] = idx_from;
    list_[idx_from] = last;
    list_[idx_to] = elt;
    index_[elt] = idx_to;
}

void IntStack::pop_back() { --size_; }
void IntStack::pop_head() { remove(list_[0]); }

int IntStack::head() const { return list_[0]; }
int IntStack::back() const { return list_[size_-1]; }

void IntStack::safe_add(const int elt)
{
	if(elt >= 0) {
		if(elt >= list_.size()) {
			reserve(elt+1);
		} 
		add(elt);
	}
}

void IntStack::add(const int elt)
{
	if( index_[elt] >= size_ ) push(elt);
}

void IntStack::push(const int elt)
{
	auto next = list_[size_];
  index_[next] = index_[elt];
  list_[index_[elt]] = next;
	index_[elt] = size_;
	list_[size_++] = elt;
}
//@}


std::ostream& IntStack::display(std::ostream& os) const
{
		os << "(";
    for( auto it = begin() ; it < end() ; ++it ) {
        os << " " << *it;
    }
    os << " )";
    return os;
}

std::ostream& operator<<(std::ostream& os, const IntStack& x)
{
    return x.display(os);
}

std::ostream& operator<<(std::ostream& os, const IntStack* x)
{
    return (x ? x->display(os) : os);
}
