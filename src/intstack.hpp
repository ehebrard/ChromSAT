
#include <iostream>

#ifndef __INTSTACK_HPP
#define __INTSTACK_HPP

/**********************************************
* IntStack
**********************************************/
/// Sparse set representation

class IntStack
{
public:
    /*!@name Parameters*/
    //@{
    /// list of values
    int* list_;
    /// current max capacity
    unsigned int index_capacity;
    unsigned int list_capacity;
    /// current size
    unsigned int size;
    /// values' indices
    unsigned int* index_;
    unsigned int* start_;
    //@}

    /*!@name Constructors*/
    //@{
    IntStack()
    {
        size = 0;
        index_capacity = 0;
        list_capacity = 0;

        list_ = NULL;
        index_ = NULL;
        start_ = NULL;
    }

    IntStack(const int lb, const int ub, const bool full = true)
    {
        initialise(lb, ub, ub - lb + 1, full);
    }

    IntStack(const int lb, const int ub, const int sz, const bool full)
    {
        initialise(lb, ub, sz, full);
    }

    IntStack(IntStack& shared, const int sz) { initialise(shared, sz); }

    void neutralise() { start_ = NULL; }

    virtual ~IntStack()
    {
        delete[] list_;
        delete[] start_;
    }

    IntStack& operator=(const IntStack& other)
    {
        size = other.size;
        index_ = other.index_;
        start_ = other.start_;
        index_capacity = other.index_capacity;

        list_ = new int[list_capacity];
        std::copy(other.list_, other.list_+list_capacity, list_);
        start_ = new unsigned[index_capacity];
        std::copy(other.start_, other.start_+index_capacity, start_);

        return *this;
    }
    IntStack& operator=(IntStack&& other)
    {
        list_ = other.list_;
        index_capacity = other.index_capacity;
        list_capacity = other.list_capacity;
        size = other.size;
        index_ = other.index_;
        start_ = other.start_;

        other.size = 0;
        other.index_capacity = 0;
        other.list_capacity = 0;

        other.list_ = NULL;
        other.index_ = NULL;
        other.start_ = NULL;
        return *this;
    };

    virtual void initialise(IntStack& shared, const int sz);
    virtual void initialise(
        const int lb, const int ub, const int sz, const bool full);

    void extend_list();
    void extend(const int new_elt);

    /*!@name Accessors*/
    //@{
    int get_min() const;

    int get_max() const;

    bool safe_contain(const int elt) const;

    bool contain(const int elt) const;

    bool contain(const int elt, const int min_idx, const int max_idx) const;

    bool empty() const;

    int next(const int elt) const;

    int prev(const int elt) const;

    int operator[](const unsigned int idx) const;

    int& operator[](const unsigned int idx);
    //@}

    /*!@name List Manipulation*/
    //@{
    int* begin();
    const int* begin() const;

    int* end();
    const int* end() const;

    int* end_mem();

    void fill();

    void clear();

    void set_to(const int elt);

    void remove(const int elt);
    void remove(const int elt, int& rk); // put elt at rank rk and increase rk

    int next();

    int pop();

    int pop_head();

    int head() const;

    int back() const;

    void init_add(const int elt);

    void add(const int elt);
    void add(const int elt, int& rk); // put elt at rank rk-1 and decrease rk

    void safe_add(const int elt);

    // create a new element that can potentially be outside the bounds
    void create(const int elt);

    void ordered_add(const int elt);

    void revert_to(const int level);

    void index();

    int get_index(const int elt);

    void move(const int elt, const int idx);

    void set(const int elt, const int idx);
    //@}

    /*!@name Miscellaneous*/
    //@{
    std::ostream& display(std::ostream& os) const;
    std::ostream& display(std::ostream& os, const int rank) const;
    std::string to_str() const;
};

std::ostream& operator<<(std::ostream& os, const IntStack& x);

#endif // __INTSTACK_HPP
