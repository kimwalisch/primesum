///
/// @file  pi_lmo1.cpp
/// @brief Simple demonstration implementation of the
///        Lagarias-Miller-Odlyzko prime summing algorithm.
///        Usually in the Lagarias-Miller-Odlyzko algorithm phi(x, a)
///        is calculated using a prime sieve but this simple
///        implementation calculates phi(x, a) using the recursive
///        formula with caching.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum.hpp>
#include <primesum-internal.hpp>
#include <generate.hpp>
#include <PhiTiny.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <int256_t.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

using namespace std;

namespace primesum {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3)) space.
///
int256_t pi_lmo1(int128_t x)
{
  if (x < 2)
    return 0;

  int64_t y = iroot<3>(x); 
  int64_t pi_y = pi_legendre(y, 1);
  int64_t c = PhiTiny::get_c(y);
  int256_t S1 = 0;
  int256_t S2 = 0;

  vector<int32_t> primes = generate_primes(y);
  vector<int32_t> lpf = generate_lpf(y);
  vector<int32_t> mu = generate_moebius(y);

  // Calculate the contribution of the ordinary leaves
  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > primes[c])
      S1 += phi_sum(x / n, c) * (mu[n] * n);
  
  // Calculate the contribution of the special leaves
  for (int64_t b = c + 1; b < pi_y; b++)
    for (int128_t m = (y / primes[b]) + 1; m <= y; m++)
      if (lpf[m] > primes[b])
        S2 -= phi_sum(x / (primes[b] * m), b - 1) * (mu[m] * m * primes[b]);

  int256_t phi = S1 + S2;
  int256_t p2 = P2(x, y, 1);
  int256_t sum = phi + prime_sum_tiny(y) - 1 - p2;

  return sum;
}

} // namespace
