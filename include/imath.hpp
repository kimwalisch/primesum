///
/// @file  imath.hpp
/// @brief Integer math functions used in primesum.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef IMATH_HPP
#define IMATH_HPP

#include <isqrt.hpp>

#include <stdint.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace primesum {

inline int64_t isquare(int64_t x)
{
  return x * x;
}

template <typename A, typename B>
inline A ceil_div(A a, B b)
{
  assert(b > 0);
  return (A) ((a + b - 1) / b);
}

template <typename T>
inline T number_of_bits(T)
{
  return (T) std::numeric_limits<T>::digits;
}

template <typename T>
inline T next_power_of_2(T x)
{
  if (x == 0)
    return 1;

  x--;
  for (T i = 1; i < number_of_bits(x); i += i)
    x |= (x >> i);

  return ++x;
}

template <typename T>
inline T prev_power_of_2(T x)
{
  for (T i = 1; i < number_of_bits(x); i += i)
    x |= (x >> i);

  return x - (x >> 1);
}

template <typename T>
inline int ilog(T x)
{
  return (int) std::log((double) x);
}

template <typename T>
inline T ipow(T x, int n)
{
  T r = 1;
  for (int i = 0; i < n; i++)
    r *= x;

  return r;
}

/// Integer nth root
template <int N, typename T>
inline T iroot(T x)
{
  if (N == 0)
    return 0;

  T r = (T) std::pow((double) x, 1.0 / N);

  // fix root too large
  for (; r > 0; r--)
  {
    if (ipow(r, N - 1) <= x / r)
      break;
  }

  // fix root too small
  while (ipow(r + 1, N - 1) <= x / (r + 1))
    r += 1;

  return r;
}

template <typename T1, typename T2, typename T3>
inline T2 in_between(T1 min, T2 x, T3 max)
{
  if (x < min)
    return (T2) min;
  if (x > max)
    return (T2) max;

  return x;
}

/// Count the number of primes <= x using binary search.
/// @pre primes[1] = 2, primes[3] = 3, ...
/// @pre x <= primes.back()
///
template <typename T1, typename T2>
inline T2 pi_bsearch(const std::vector<T1>& primes, T2 x)
{
  assert(primes.size() < 2 || primes[1] == 2);
  auto start = primes.begin() + 1;
  return (T2) (std::upper_bound(start, primes.end(), x) - start);
}

} // namespace

#endif
