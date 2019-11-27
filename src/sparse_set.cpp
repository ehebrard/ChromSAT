
#include <assert.h>

#include "sparse_set.hpp"

sparse_set::sparse_set(const size_t n) {
  end_ = 0;
	start_ = 0;
  reserve(n);
}

void sparse_set::reserve(const size_t n) {
  while (list_.size() < n) {
    index_.push_back(list_.size());
    list_.push_back(list_.size());
  }
}

// void sparse_set::re_back(size_t &stamp) { ++end_; }
void sparse_set::save_back(size_t &stamp) { stamp = end_; }
void sparse_set::restore_back(const size_t stamp) { end_ = stamp; }

// void sparse_set::re_front(size_t &stamp) { --start_; }
void sparse_set::save_front(size_t &stamp) { stamp = start_; }
void sparse_set::restore_front(const size_t stamp) { start_ = stamp; }

// void sparse_set::save(size_t &stamp1, size_t &stamp2) { stamp1 = end_; stamp2 = start_; }
// void sparse_set::restore(const size_t stamp1, const size_t stamp2) { end_ = stamp1; start_ = stamp2; }


//@}

/*!@name Accessors*/
//@{

bool sparse_set::safe_contain(const int elt) const {
  if (elt >= 0 && (size_t)elt < index_.size())
    return contain(elt);
  return false;
}

bool sparse_set::contain(const int elt) const { return index_[elt] < end_ and index_[elt] >= start_; }

size_t sparse_set::size() const { return end_ - start_; }
size_t sparse_set::num_front() const { return start_; }
size_t sparse_set::num_back() const { return capacity() - end_; }
size_t sparse_set::capacity() const { return index_.size(); }

bool sparse_set::empty() const { return end_ == start_; }

int sparse_set::next(const int elt) const {
  size_t idx = index_[elt] + 1;
  return (idx < end_ ? list_[idx] : elt);
}
int sparse_set::prev(const int elt) const {
  size_t idx = index_[elt];
  return (idx > start_ ? list_[idx - 1] : elt);
}

int sparse_set::operator[](const size_t idx) const { return list_[idx+start_]; }

int &sparse_set::operator[](const size_t idx) { return list_[idx+start_]; }
//@}

/*!@name List Manipulation*/
//@{
std::vector<int>::iterator sparse_set::begin() { return list_.begin() + start_; }
std::vector<int>::reverse_iterator sparse_set::rbegin() {
  return list_.rend() - end_;
}

std::vector<int>::iterator sparse_set::end() { return list_.begin() + end_; }
std::vector<int>::reverse_iterator sparse_set::rend() { return list_.rend() + start_; }

std::vector<int>::const_iterator sparse_set::begin() const {
  return list_.begin() + start_;
}
std::vector<int>::const_reverse_iterator sparse_set::rbegin() const {
  return list_.rend() + end_;
}

std::vector<int>::const_iterator sparse_set::end() const {
  return list_.begin() + end_;
}
std::vector<int>::const_reverse_iterator sparse_set::rend() const {
  return list_.rend() + start_;
}

std::vector<int>::iterator sparse_set::begin_back() { return list_.begin() + end_; }
std::vector<int>::reverse_iterator sparse_set::rbegin_back() {
  return list_.rbegin();
}

std::vector<int>::iterator sparse_set::end_back() { return list_.end(); }
std::vector<int>::reverse_iterator sparse_set::rend_back() { return list_.rend() - end_; }

std::vector<int>::const_iterator sparse_set::begin_back() const {
  return list_.begin() + end_;
}
std::vector<int>::const_reverse_iterator sparse_set::rbegin_back() const {
  return list_.rbegin();
}

std::vector<int>::const_iterator sparse_set::end_back() const {
  return list_.end();
}
std::vector<int>::const_reverse_iterator sparse_set::rend_back() const {
  return list_.rend() - end_;
}

std::vector<int>::iterator sparse_set::begin_front() { return list_.begin(); }
std::vector<int>::reverse_iterator sparse_set::rbegin_front() {
  return list_.rend() - start_;
}

std::vector<int>::iterator sparse_set::end_front() { return list_.begin() + start_; }
std::vector<int>::reverse_iterator sparse_set::rend_front() { return list_.rend(); }

std::vector<int>::const_iterator sparse_set::begin_front() const {
  return list_.begin();
}
std::vector<int>::const_reverse_iterator sparse_set::rbegin_front() const {
  return list_.rend() - start_;
}

std::vector<int>::const_iterator sparse_set::end_front() const {
  return list_.begin() + start_;
}
std::vector<int>::const_reverse_iterator sparse_set::rend_front() const {
  return list_.rend();
}

// std::vector<int>::iterator sparse_set::begin_after() { return end(); }
// std::vector<int>::reverse_iterator sparse_set::rbegin_after() {
//   return list_.rend();
// }
//
// std::vector<int>::iterator sparse_set::end_after() { return list_.end(); }
// std::vector<int>::reverse_iterator sparse_set::rend_after() { return rbegin(); }
//
// std::vector<int>::const_iterator sparse_set::begin_after() const {
//   return end();
// }
// std::vector<int>::const_reverse_iterator sparse_set::rbegin_after() const {
//   return list_.rend();
// }
//
// std::vector<int>::const_iterator sparse_set::end_after() const {
//   return list_.end();
// }
// std::vector<int>::const_reverse_iterator sparse_set::rend_after() const {
//   return rend();
// }

void sparse_set::fill() { end_ = list_.size(); start_ = 0; }

void sparse_set::clear() { end_ = 0; start_ = 0; }

// void sparse_set::set_size(const int s) { end_ = s; }

// void sparse_set::safe_remove_back(const int elt) {
//   if (elt >= 0) {
//     if (static_cast<size_t>(elt) >= list_.size()) {
//       reserve(elt + 1);
//     }
//     remove_back(elt);
//   }
// }

void sparse_set::remove_back(const int elt) {
  if (index_[elt] < end_ and index_[elt] >= start_)
    pull_back(elt);
}

void sparse_set::pull_back(const int elt) {
  auto last = list_[--end_];
  index_[last] = index_[elt];
  list_[index_[elt]] = last;
  list_[end_] = elt;
  index_[elt] = end_;
}

// void sparse_set::safe_remove_front(const int elt) {
//   if (elt >= 0) {
//     if (static_cast<size_t>(elt) >= list_.size()) {
//       reserve(elt + 1);
//     }
//     remove_front(elt);
//   }
// }

void sparse_set::remove_front(const int elt) {
  if (index_[elt] < end_ and index_[elt] >= start_)
    pull_front(elt);
}

void sparse_set::pull_front(const int elt) {
  auto first = list_[start_];
  index_[first] = index_[elt];
  list_[index_[elt]] = first;
  list_[start_] = elt;
  index_[elt] = start_++;
}

// void sparse_set::move(const int elt, const int idx_to) {
//   auto idx_from = index_[elt];
//
// 	// assert(idx_from )
//
//   // assert(index_[elt] <= static_cast<size_t>(idx_to));
//
//   auto last = list_[idx_to];
//   index_[last] = idx_from;
//   list_[idx_from] = last;
//   list_[idx_to] = elt;
//   index_[elt] = idx_to;
// }

void sparse_set::pop_back() { --end_; }
void sparse_set::pop_front() { ++start_; }

int sparse_set::front() const { return list_[start_]; }
int sparse_set::back() const { return list_[end_ - 1]; }

void sparse_set::safe_add(const int elt) {
  if (elt >= 0) {
    if (static_cast<size_t>(elt) >= list_.size()) {
      reserve(elt + 1);
    }
    add(elt);
  }
}

void sparse_set::add(const int elt) {
  if (index_[elt] >= end_)
    push_back(elt);
	else if(index_[elt] < start_)
		push_front(elt);
}

void sparse_set::push_back(const int elt) {
  auto next = list_[end_];
  index_[next] = index_[elt];
  list_[index_[elt]] = next;
  index_[elt] = end_;
  list_[end_++] = elt;
}

void sparse_set::push_front(const int elt) {
  auto next = list_[--start_];
  index_[next] = index_[elt];
  list_[index_[elt]] = next;
  index_[elt] = start_;
  list_[start_] = elt;
}

int sparse_set::index(const int elt) const { return index_[elt]; }
//@}

std::ostream &sparse_set::display(std::ostream &os) const {
  os << "(";
  for (auto it = begin(); it < end(); ++it) {
    os << " " << *it;
  }
  os << " )";
  return os;
}

std::ostream &operator<<(std::ostream &os, const sparse_set &x) {
  return x.display(os);
}

std::ostream &operator<<(std::ostream &os, const sparse_set *x) {
  return (x ? x->display(os) : os);
}
