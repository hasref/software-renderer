#ifndef INTDEFINES_H_
#define INTDEFINES_H_

#include <cstdint>

/*
 * This is mainly for ergonomics.
 * I just find u8 to be a lot easier to type than uint8_t.
 * I directly borrowed these from rust.
 *
 */

using u8 = std::uint8_t;
using i32 = std::int32_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;
using usize = std::size_t;

#endif
