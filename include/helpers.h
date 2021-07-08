#ifndef HELPERS_H_
#define HELPERS_H_

#include <type_traits>

template <typename T>
bool in_range_inclusive(const T& val, const T& a, const T& b) {
  if constexpr (std::is_arithmetic_v<T>) {
    return (val >= a) && (val <= b);
  } else {
    static_assert(std::is_arithmetic_v<T>);
  }
}

#endif
