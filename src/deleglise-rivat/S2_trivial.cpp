///
/// @file  S2_trivial.cpp
/// @brief Calculate the contribution of the trivial special leaves
///        in parallel using OpenMP.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primesum.hpp>
#include <primesum-internal.hpp>
#include <primesieve.hpp>
#include <generate.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primesum;

namespace {

template <typename T>
T S2_trivial_OpenMP(T x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int threads)
{
  int64_t thread_threshold = ipow(10, 7);
  threads = validate_threads(threads, y, thread_threshold);

  PiTable pi(y);
  vector<int64_t> prime_sums = generate_prime_sums(y);
  int64_t sqrtz = isqrt(z);
  int64_t prime_c = nth_prime(c);

  T s2_trivial = 0;

  // Find all trivial leaves: n = primes[b] * primes[l]
  // which satisfy phi(x / n), b - 1) = 1
  #pragma omp parallel for num_threads(threads) reduction(+: s2_trivial)
  for (int64_t i = 0; i < threads; i++)
  {
    int64_t start = max(prime_c, sqrtz) + 1;
    int64_t thread_interval = ceil_div(y - start, threads);
    start += thread_interval * i;
    int64_t stop = min(start + thread_interval, y);
    primesieve::iterator iter(start - 1, stop);
    T prime;

    while ((prime = iter.next_prime()) < stop)
    {
      int64_t xn = (int64_t) max(x / (prime * prime), prime);
      T diff = prime_sums[pi[y]] - prime_sums[pi[xn]];
      s2_trivial += prime * diff;
    }
  }

  return s2_trivial;
}

} // namespace

namespace primesum {

maxint_t S2_trivial(maxint_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int threads)
{
  print("");
  print("=== S2_trivial(x, y) ===");
  print("Computation of the trivial special leaves");
  print(x, y, c, threads);

  double time = get_wtime();
  maxint_t s2_trivial = S2_trivial_OpenMP(x, y, z, c, threads);

  print("S2_trivial", s2_trivial, time);
  return s2_trivial;
}

} // namespace
