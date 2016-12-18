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

#include <stdint.h>
#include <algorithm>
#include <vector>

using namespace std;

namespace primesum {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3)) space.
///
maxint_t pi_lmo1(maxint_t x)
{
  if (x < 2)
    return 0;

  int64_t y = iroot<3>(x); 
  int64_t pi_y = pi_legendre(y, 1);
  int64_t c = PhiTiny::get_c(y);
  maxint_t S1 = 0;
  maxint_t S2 = 0;

  vector<int32_t> primes = generate_primes(y);
  vector<int32_t> lpf = generate_least_prime_factors(y);
  vector<int32_t> mu = generate_moebius(y);

  // Calculate the contribution of the ordinary leaves
  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > primes[c])
      S1 += mu[n] * n * phi_sum(x / n, c);

  // Calculate the contribution of the special leaves
  for (int64_t b = c + 1; b < pi_y; b++)
    for (maxint_t m = (y / primes[b]) + 1; m <= y; m++)
      if (lpf[m] > primes[b])
        S2 -= mu[m] * m * primes[b] * phi_sum(x / (primes[b] * m), b - 1);

  maxint_t phi = S1 + S2;
  maxint_t p2 = P2(x, y, 1);
  maxint_t sum = phi + prime_sum_tiny(y) - 1 - p2;

  return sum;
}

} // namespace
