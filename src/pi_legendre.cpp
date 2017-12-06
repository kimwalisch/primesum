///
/// @file  pi_legendre.cpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum.hpp>
#include <primesum-internal.hpp>
#include <imath.hpp>

#include <stdint.h>

namespace primesum {

/// Count the number of primes <= x using Legendre's formula.
/// Run time: O(x)
/// Memory usage: O(x^(1/2))
///
int64_t pi_legendre(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  // temporarily disable printing
  bool print = is_print();
  set_print(false);

  int64_t a = pi_primesieve(isqrt(x));
  int64_t sum = phi(x, a, threads) + a - 1;

  set_print(print);
  return sum;
}

} // namespace
