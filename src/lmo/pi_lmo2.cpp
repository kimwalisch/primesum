///
/// @file  pi_lmo2.cpp
/// @brief Simple implementation of the Lagarias-Miller-Odlyzko prime
///        counting algorithm. This implementation uses the sieve
///        of Eratosthenes (without segmentation) to calculate S2(x).
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum-internal.hpp>
#include <pmath.hpp>
#include <generate.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

using namespace std;
using namespace primesum;

namespace {

/// Calculate the contribution of the special leaves.
/// This implementation uses the sieve of Eratosthenes (without
/// segmentation). Space complexity: O(x^(2/3) * (log x)^2).
/// @pre y > 0 && c > 1
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t c,
           vector<int32_t>& primes,
           vector<int32_t>& lpf,
           vector<int32_t>& mu)
{
  int64_t limit = x / y + 1;
  int64_t pi_y = pi_bsearch(primes, y);
  int64_t S2_result = 0;
  int64_t b = 1;
  vector<char> sieve(limit, 1);

  // phi(y, b) nodes with b <= c do not contribute to S2, so we
  // simply sieve out the multiples of the first c primes
  for (; b <= c; b++)
  {
    int64_t prime = primes[b];
    for (int64_t k = prime; k < limit; k += prime)
      sieve[k] = 0;
  }

  for (; b < pi_y; b++)
  {
    int64_t prime = primes[b];
    int64_t i = 1;
    int64_t phi = 0;

    for (int64_t m = y; m > y / prime; m--)
    {
      if (mu[m] != 0 && prime < lpf[m])
      {
        // We have found a special leaf, compute it's contribution
        // phi(x / (primes[b] * m), b - 1) by counting the number
        // of unsieved elements <= x / (primes[b] * m) after having
        // removed the multiples of the first b - 1 primes
        //
        for (int64_t xn = x / (prime * m); i <= xn; i++)
          phi += i * sieve[i];

        S2_result -= mu[m] * m * prime * phi;
      }
    }

    // Remove the multiples of (b)th prime
    for (int64_t k = prime; k < limit; k += prime * 2)
      sieve[k] = 0;
  }

  return S2_result;
}

} // namespace

namespace primesum {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(2/3) / (log x)^2) space.
///
int64_t pi_lmo2(int64_t x)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_lmo(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t c = PhiTiny::get_c(y);
  maxint_t p2 = P2(x, y, 1);

  vector<int32_t> primes = generate_primes(y);
  vector<int32_t> lpf = generate_least_prime_factors(y);
  vector<int32_t> mu = generate_moebius(y);

  maxint_t s1 = S1(x, y, c, 1);
  maxint_t s2 = S2(x, y, c, primes, lpf, mu);
  maxint_t phi = s1 + s2;
  maxint_t sum = phi + prime_sum_tiny(y) - 1 - p2;

  return sum.convert_to<int64_t>();
}

} // namespace primesum
