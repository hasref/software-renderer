#ifndef ARRAY_2D_
#define ARRAY_2D_

#include <vector>

#include "intdefines.h"

template <typename T>
class Array2d {
 private:
  const usize num_rows{};
  const usize num_cols{};
  std::vector<T> arr{};

 public:
  Array2d(usize num_rows, usize num_cols)
      : num_rows(num_rows), num_cols(num_cols), arr(num_rows * num_cols, T{}) {}

  explicit Array2d(const T& fill_value)
      : arr(num_rows * num_cols, fill_value) {}

  inline T& operator()(usize row, usize col) {
    return arr[row * num_cols + col];
  }

  inline const T& operator()(usize row, usize col) const {
    return arr[row * num_cols + col];
  }

  usize get_numrows() const { return num_rows; }
  usize get_numcols() const { return num_cols; }

  // all of the iterator methods!
  // we may need some sort of a slice later.
  auto begin() { return arr.begin(); }
  auto end() { return arr.end(); }
  auto cbegin() const { return arr.cbegin(); }
  auto cend() const { return arr.cend(); }
  auto rbegin() { return arr.rbegin(); }
  auto rend() { return arr.rend(); }
  auto crbegin() const { return arr.crbegin(); }
  auto crend() const { return arr.crend(); }
};

#endif
