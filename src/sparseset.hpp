
#include <iostream>

#include <vector>

#ifndef __SPARSESET_HPP
#define __SPARSESET_HPP

/**********************************************
 * sparseset
 **********************************************/
/// Sparse set representation

#define NOINDEX 0xfffffff

class sparseset
{

public:
    // static const size_t NOVAL{0xfffffff};

    /*!@name Parameters*/
    //@{
    /// list of values
    std::vector<int> list_;
    size_t size_;

    /// values' indices
    std::vector<size_t>* index_;
    //@}

    /*!@name Constructors*/
    //@{
    explicit sparseset();
    void binds(std::vector<size_t>* index);

    //     explicit sparseset(std::vector<size_t> *i, const size_t n);
    // explicit sparseset(sparseset& s);
    //
    //     void reserve(const size_t n);

    /*!@name Accessors*/
    //@{
    bool safe_contain(const int elt) const;
    bool contain(const int elt) const;

    size_t size() const;
    bool empty() const;

    int next(const int elt) const;

    int prev(const int elt) const;

    int operator[](const size_t idx) const;

    int& operator[](const size_t idx);
    //@}

    /*!@name List Manipulation*/
    //@{
    std::vector<int>::iterator begin();
    std::vector<int>::reverse_iterator rbegin();

    std::vector<int>::iterator end();
    std::vector<int>::reverse_iterator rend();

    std::vector<int>::const_iterator begin() const;
    std::vector<int>::const_reverse_iterator rbegin() const;

    std::vector<int>::const_iterator end() const;
    std::vector<int>::const_reverse_iterator rend() const;

    void fill();

    void clear();

    void resize(const size_t n);

    void remove(const int elt);

    void move_up(const int elt, const int idx);

    void pop_back();

    void pop_head();

    int head() const;

    int back() const;

    void push(const int elt);

    void add(const int elt);

    void safe_add(const int elt);

    int index(const int elt) const;

    void save(size_t& stamp);
    void restore(const size_t stamp);
    //@}

    /*!@name Miscellaneous*/
    //@{
    std::ostream& display(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const sparseset& x);

#endif // __SPARSESET_HPP
