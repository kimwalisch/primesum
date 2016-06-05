///
/// @file  pi_lmo3.cpp
/// @brief Implementation of the Lagarias-Miller-Odlyzko prime
///        counting algorithm. This implementation uses the segmented
///        sieve of Eratosthenes and a special tree data structure
///        for faster counting in S2(x).
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

/// Cross-off the multiples of prime in the sieve array.
/// For each element that is unmarked the first time update
/// the special counters tree data structure.
///
template <typename T>
void cross_off(int64_t prime,
               int64_t low,
               int64_t high,
               int64_t& next_multiple,
               T& sieve)
{
  int64_t m = next_multiple;

  for (; m < high; m += prime * 2)
    sieve[m - low] = 0;

  next_multiple = m;
}

/// Calculate the contribution of the special leaves.
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
  int64_t segment_size = next_power_of_2(isqrt(limit));
  int64_t pi_y = pi_bsearch(primes, y);
  int64_t S2_result = 0;

  vector<char> sieve(segment_size);
  vector<int64_t> next(primes.begin(), primes.end());
  vector<int64_t> phi(primes.size(), 0);

  // Segmented sieve of Eratosthenes
  for (int64_t low = 1; low < limit; low += segment_size)
  {
    fill(sieve.begin(), sieve.end(), 1);

    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);

    // phi(y, b) nodes with b <= c do not contribute to S2, so we
    // simply sieve out the multiples of the first c primes
    for (int64_t b = 1; b <= c; b++)
    {
      int64_t k = next[b];
      for (int64_t prime = primes[b]; k < high; k += prime)
        sieve[k - low] = 0;
      next[b] = k;
    }

    for (int64_t b = c + 1; b < pi_y; b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t max_m = min(x / (prime * low), y);
      int64_t i = 0;

      // Obviously if (prime >= max_m) then (prime >= lpf[max_m])
      // if so then (prime < lpf[m]) will always evaluate to
      // false and no special leaves are possible
      if (prime >= max_m)
        break;

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (mu[m] != 0 && prime < lpf[m])
        {
          int64_t xn = x / (prime * m);
          int64_t stop = xn - low;

          for (; i <= stop; i++)
            phi[b] += (low + i) * sieve[i];

          S2_result -= mu[m] * m * prime * phi[b];
        }
      }

      for (; i < high - low; i++)
        phi[b] += (low + i) * sieve[i];

      // Remove the multiples of (b)th prime
      cross_off(prime, low, high, next[b], sieve);
    }
  }

  return S2_result;
}

} // namespace

namespace primesum {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo3(int64_t x)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_lmo(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t c = PhiTiny::get_c(y);
  maxint_t p2 = P2(x, y, 1);

  vector<int32_t> mu = generate_moebius(y);
  vector<int32_t> lpf = generate_least_prime_factors(y);
  vector<int32_t> primes = generate_primes(y);

  maxint_t s1 = S1(x, y, c, 1);
  maxint_t s2 = S2(x, y, c, primes, lpf, mu);
  maxint_t phi = s1 + s2;
  maxint_t sum = phi + prime_sum_tiny(y) - 1 - p2;

  return sum.convert_to<int64_t>();
}

} // namespace primesum
