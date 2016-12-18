///
/// @file  phi.cpp
/// @brief The PhiCache class calculates the partial sieve function
///        (a.k.a. Legendre-sum) using the recursive formula:
///        phi(x, a) = phi(x, a - 1) - phi(x / primes[a], a - 1).
///        phi(x, a) counts the numbers <= x that are not divisible by
///        any of the first a primes. The algorithm used is an
///        optimized version of the algorithm described in Tomás
///        Oliveira e Silva's paper [1]. I have added 5 optimizations
///        to my implementation which significantly speed up the
///        calculation:
///
///        * Cache results of phi(x, a)
///        * Calculate phi(x, a) using formula [2] if a <= 6
///        * Calculate phi(x, a) using pi(x) lookup table
///        * Calculate all phi(x, a) = 1 upfront
///        * Stop recursion at c instead of 1
///
///       [1] Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///           method, Revista do DETUA, vol. 4, no. 6, March 2006, p. 761.
///           http://sweet.ua.pt/tos/bib/5.4.pdf
///       [2] phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a)
///           with pp = 2 * 3 * ... * prime[a] 
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primesum-internal.hpp>
#include <primesum.hpp>
#include <generate.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <fast_div.hpp>
#include <min_max.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>
#include <limits>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primesum;

namespace {

/// Cache phi(x, a) results if a <= MAX_A
const int MAX_A = 500;

/// Keep the cache size below MAX_BYTES per thread
const int MAX_BYTES = 16 << 20;

class PhiCache
{
public:
  PhiCache(vector<int32_t>& primes,
           PiTable& pi) :
    primes_(primes),
    pi_(pi),
    bytes_(0)
  {
    size_t max_size = MAX_A + 1;
    size_t size = min(primes.size(), max_size);
    cache_.resize(size);
  }

  /// Calculate phi(x, a) using the recursive formula:
  /// phi(x, a) = phi(x, a - 1) - phi(x / primes_[a], a - 1)
  ///
  template <int SIGN>
  int64_t phi(int64_t x, int64_t a)
  {
    int64_t sum = 0;

    if (x <= primes_[a])
      sum = SIGN;
    else if (is_phi_tiny(a))
      sum = phi_tiny(x, a) * SIGN;
    else if (is_phi_by_pix(x, a))
      sum = (pi_[x] - a + 1) * SIGN;
    else
    {
      int64_t sqrtx = isqrt(x);
      int64_t pi_sqrtx = a;

      if (sqrtx < pi_.size() && sqrtx < primes_[a])
        pi_sqrtx = pi_[sqrtx];

      // Move out of the loop the calculations where phi(x2, a2) = 1
      // phi(x, a) = 1 if primes_[a] >= x
      // x2 = x / primes_[a2 + 1]
      // phi(x2, a2) = 1 if primes_[a2] >= x / primes_[a2 + 1]
      // phi(x2, a2) = 1 if primes_[a2] >= sqrt(x)
      // phi(x2, a2) = 1 if a2 >= pi(sqrt(x))
      // \sum_{a2 = pi(sqrt(x))}^{a - 1} phi(x2, a2) = a - pi(sqrt(x))
      //
      sum = (a - pi_sqrtx) * -SIGN;

      // phi(x, c) = phi(x, 1) - \sum_{a2 = 1}^{c - 1} phi(x / primes_[a2 + 1], a2)
      int64_t c = min(PhiTiny::max_a(), pi_sqrtx);
      sum += phi_tiny(x, c) * SIGN;

      for (int64_t a2 = c; a2 < pi_sqrtx; a2++)
      {
        int64_t x2 = fast_div(x, primes_[a2 + 1]);
        if (is_cached(x2, a2))
          sum += cache_[a2][x2] * -SIGN;
        else
          sum += phi<-SIGN>(x2, a2);
      }
    }

    if (write_to_cache(x, a))
      cache_[a][x] = (uint16_t) (sum * SIGN);

    return sum;
  }

private:
  vector<vector<uint16_t> > cache_;
  vector<int32_t>& primes_;
  PiTable& pi_;
  int64_t bytes_;

  int64_t cache_size(int64_t a) const
  {
    return (int64_t) cache_[a].size();
  }

  bool is_phi_by_pix(int64_t x, int64_t a) const
  {
    return x < pi_.size() &&
           x < isquare(primes_[a + 1]);
  }

  bool is_cached(int64_t x, int64_t a) const
  {
    return a <= MAX_A && 
           x < cache_size(a) && 
           cache_[a][x] != 0;
  }

  bool write_to_cache(int64_t x, int64_t a)
  {
    if (a > MAX_A || 
        x > numeric_limits<uint16_t>::max())
      return false;

    // we need to increase cache size
    if (x >= cache_size(a))
    {
      if (bytes_ > MAX_BYTES)
        return false;
      bytes_ += (x + 1 - cache_size(a)) * 2;
      cache_[a].resize(x + 1, 0);
    }

    return true;
  }
};

} // namespace

namespace primesum {

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, int threads)
{
  if (x < 1) return 0;
  if (a > x) return 1;
  if (a < 1) return x;

  print("");
  print("=== phi(x, a) ===");
  print("Count the numbers <= x coprime to the first a primes");

  double time = get_wtime();
  int64_t sum = 0;

  if (is_phi_tiny(a))
    sum = phi_tiny(x, a);
  else
  {
    vector<int32_t> primes = generate_n_primes(a);

    if (primes.at(a) >= x)
      sum = 1;
    else
    {
      // use a large pi(x) lookup table for speed
      int64_t sqrtx = isqrt(x);
      PiTable pi(max(sqrtx, primes[a]));
      PhiCache cache(primes, pi);

      int64_t pi_sqrtx = min(pi[sqrtx], a); 
      sum = x - a + pi_sqrtx;

      int64_t p14 = ipow((int64_t) 10, 14);
      int64_t thread_threshold = p14 / primes[a];
      threads = ideal_num_threads(threads, x, thread_threshold);

      // this loop scales only up to about 8 CPU cores
      threads = min(8, threads);

      #pragma omp parallel for schedule(dynamic, 16) \
          num_threads(threads) firstprivate(cache) reduction(+: sum)
      for (int64_t a2 = 0; a2 < pi_sqrtx; a2++)
        sum += cache.phi<-1>(x / primes[a2 + 1], a2);
    }
  }

  print("phi", sum, time);
  return sum;
}

} // namespace
