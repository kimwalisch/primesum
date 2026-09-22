///
/// @file  popcnt.hpp
/// @brief Functions to count the number of 1 bits in a 64-bit
///        variable using compiler intrinsics or a portable fallback.
///        No runtime CPU detection is performed.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef POPCNT_HPP
#define POPCNT_HPP

#include <macros.hpp>
#include <stdint.h>

// GCC & Clang
#if defined(__GNUC__) || \
    __has_builtin(__builtin_popcountl)

namespace {

ALWAYS_INLINE uint64_t popcnt64(uint64_t x)
{
#if __cplusplus >= 201703L
  if constexpr(sizeof(int) >= sizeof(uint64_t))
    return (uint64_t) __builtin_popcount(x);
  else if constexpr(sizeof(long) >= sizeof(uint64_t))
    return (uint64_t) __builtin_popcountl(x);
  else if constexpr(sizeof(long long) >= sizeof(uint64_t))
    return (uint64_t) __builtin_popcountll(x);
#else
    return (uint64_t) __builtin_popcountll(x);
#endif
}

} // namespace

#elif defined(_MSC_VER) && \
      defined(_M_X64) && \
      __has_include(<intrin.h>)

#include <intrin.h>

namespace {

ALWAYS_INLINE uint64_t popcnt64(uint64_t x)
{
  return __popcnt64(x);
}

} // namespace

#elif defined(_MSC_VER) && \
      defined(_M_IX86) && \
      __has_include(<intrin.h>)

#include <intrin.h>

namespace {

ALWAYS_INLINE uint64_t popcnt64(uint64_t x)
{
  return __popcnt(uint32_t(x)) +
         __popcnt(uint32_t(x >> 32));
}

} // namespace

#elif __cplusplus >= 202002L && \
      __has_include(<bit>)

#include <bit>

namespace {

/// We only use the C++ standard library as a fallback if there
/// are no compiler intrinsics available for POPCNT.
/// Compiler intrinsics often generate faster assembly.
ALWAYS_INLINE uint64_t popcnt64(uint64_t x)
{
  return std::popcount(x);
}

} // namespace

#else

namespace {

/// This uses fewer arithmetic operations than any other known
/// implementation on machines with fast multiplication.
/// It uses 12 arithmetic operations, one of which is a multiply.
/// http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
///
ALWAYS_INLINE uint64_t popcnt64(uint64_t x)
{
  uint64_t m1 = 0x5555555555555555;
  uint64_t m2 = 0x3333333333333333;
  uint64_t m4 = 0x0F0F0F0F0F0F0F0F;
  uint64_t h01 = 0x0101010101010101;

  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;

  return (x * h01) >> 56;
}

} // namespace

#endif

#endif // POPCNT_HPP
