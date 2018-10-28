
#include <iostream>

#include <vector>

#ifndef __PARTITION_HPP
#define __PARTITION_HPP

/**********************************************
 * partition
 **********************************************/
/// Sparse set representation

class partition
{
private:
  /// values' indices
  std::vector<size_t> index_;
  //@}
	
public:
    /*!@name Parameters*/
    //@{
    /// list of values
    std::vector<std::vector<int>> bag;

    /*!@name Constructors*/
    //@{
    explicit partition(const size_t n = 0, const size_t m = 1);

    void resize(const size_t n);

		size_t size();
		
		/*!@name List Manipulation*/
    //@{
    void move(const int elt, const int from, const int to);
		void add(const int elt, const int to);
		void swap(const int a, const int b);
		void remove(const int a);
    //@}

    /*!@name Miscellaneous*/
    //@{
    std::ostream& display(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const partition& x);

#endif // __PARTITION_HPP
