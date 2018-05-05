///
/// @file  S2_easy.cpp
/// @brief Calculate the contribution of the clustered easy leaves
///        and the sparse easy leaves in parallel using OpenMP
///        (Deleglise-Rivat algorithm).
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primesum-internal.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <int256_t.hpp>
#include <min_max.hpp>
#include <imath.hpp>
#include <S2Status.hpp>
#include <S2.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primesum;

namespace {

/// Calculate the contribution of the clustered easy leaves
/// and the sparse easy leaves.
/// @param T  either int64_t or uint128_t.
///
template <typename res_t, typename Primes, typename PrimeSums>
res_t S2_easy_OpenMP(uint128_t x,
                     int64_t y,
                     int64_t z,
                     int64_t c,
                     Primes& primes,
                     PrimeSums& prime_sums,
                     int threads)
{
  res_t s2_easy = 0;
  int64_t x13 = iroot<3>(x);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);
  using PS = typename PrimeSums::value_type;

  PiTable pi(y);
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_x13 = pi[x13];
  S2Status status(x);

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: s2_easy)
  for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
  {
    int64_t prime = primes[b];
    uint128_t x2 = x / prime;
    int64_t min_trivial = min(x2 / prime, y);
    int64_t min_clustered = (int64_t) isqrt(x2);
    int64_t min_sparse = z / prime;
    int64_t min_hard = max(y / prime, prime);

    min_clustered = in_between(min_hard, min_clustered, y);
    min_sparse = in_between(min_hard, min_sparse, y);

    int64_t l = pi[min_trivial];
    int64_t pi_min_clustered = pi[min_clustered];
    int64_t pi_min_sparse = pi[min_sparse];

    // Find all clustered easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) == phi(x / m, b - 1)
    // where phi(x / n, b - 1) = pi(x / n) - b + 2
    while (l > pi_min_clustered)
    {
      int64_t xn = (int64_t) fast_div(x2, primes[l]);
      int64_t phi_xn = pi[xn] - b + 2;
      res_t phi_xn_sum = prime_sums[pi[xn]] + 1 - prime_sums[b - 1];
      int64_t xm = (int64_t) fast_div(x2, primes[b + phi_xn - 1]);
      xm = max(xm, min_clustered);
      int64_t l2 = pi[xm];
      s2_easy += (phi_xn_sum * prime) * (prime_sums[l] - prime_sums[l2]);
      l = l2;
    }

    // Find all sparse easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
    for (; l > pi_min_sparse; l--)
    {
      int64_t xn = (int64_t) fast_div(x2, primes[l]);
      res_t phi = prime_sums[pi[xn]] + 1 - prime_sums[b - 1];
      s2_easy += phi * ((PS) prime * primes[l]);
    }

    if (is_print())
      status.print(b, pi_x13);
  }

  return s2_easy;
}

} // namespace

namespace primesum {

int256_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads)
{
  print("");
  print("=== S2_easy(x, y) ===");
  print("Computation of the easy special leaves");
  print(x, y, c, threads);

  double time = get_time();
  int256_t s2_easy;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(y);
    auto prime_sums = generate_prime_sums<uint64_t>(y);

    if (x <= numeric_limits<uint64_t>::max())
      s2_easy = S2_easy_OpenMP<int128_t>(x, y, z, c, primes, prime_sums, threads);
    else
      s2_easy = S2_easy_OpenMP<int256_t>(x, y, z, c, primes, prime_sums, threads);
  }
  else
  {
    auto primes = generate_primes<int64_t>(y);
    auto prime_sums = generate_prime_sums<int128_t>(y);

    s2_easy = S2_easy_OpenMP<int256_t>(x, y, z, c, primes, prime_sums, threads);
  }

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

} // namespace
