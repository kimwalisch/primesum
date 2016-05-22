///
/// @file  bit_scan_forward.hpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BIT_SCAN_FORWARD_HPP
#define BIT_SCAN_FORWARD_HPP

#include <cassert>
#include <stdint.h>

namespace primesum {

/// @return Index of the first set bit
inline uint64_t bit_scan_forward(uint64_t x)
{
  assert (x != 0);

#if defined(__GNUC__) && \
    defined(__x86_64__)
  asm ("bsfq %0, %0" : "=r" (x) : "0" (x));
  return x;
#else
  // GCC emits the bsfq instruction on x86-64 when using
  // __builtin_ctzll(x) but the assembly code above is still faster
  return __builtin_ctzll(x);
#endif
}

} // namespace

#endif
