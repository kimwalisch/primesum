///
/// @file  S1.cpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S1.hpp>
#include <primesum-internal.hpp>
#include <generate.hpp>
#include <imath.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primesum;

namespace {

/// Recursively iterate over the square free numbers coprime to the
/// first b primes and calculate the sum of the ordinary leaves.
/// This algorithm is based on section 2.2 of the paper:
/// Douglas Staple, "The Combinatorial Algorithm For Computing pi(x)",
/// arXiv:1503.01839, 6 March 2015.
///
template <int MU, typename T, typename P>
T S1_OpenMP_thread(T x,
                   int64_t y,
                   int64_t b,
                   int64_t c,
                   T square_free,
                   vector<P>& primes)
{
  T s1_sum = 0;

  for (b += 1; b < (int64_t) primes.size(); b++)
  {
    T next = square_free * primes[b];
    if (next > y) break;
    s1_sum += MU * next * phi_sum(x / next, c);
    s1_sum += S1_OpenMP_thread<-MU>(x, y, b, c, next, primes);
  }

  return s1_sum;
}

/// Calculate the contribution of the ordinary leaves in parallel.
/// Run time: O(y * log(log(y))) operations.
/// Space complexity: O(y / log(y)).
///
template <typename X, typename Y>
X S1_OpenMP_master(X x,
                   Y y,
                   int64_t c,
                   int threads)
{
  int64_t thread_threshold = ipow(10, 6);
  threads = validate_threads(threads, y, thread_threshold);
  vector<Y> primes = generate_primes<Y>(y);
  X s1_sum = phi_sum(x, c);

  #pragma omp parallel for schedule(static, 1) num_threads(threads) reduction (+: s1_sum)
  for (int64_t b = c + 1; b < (int64_t) primes.size(); b++)
  {
    s1_sum -= primes[b] * phi_sum(x / primes[b], c);
    s1_sum += S1_OpenMP_thread<1>(x, y, b, c, (X) primes[b], primes);
  }

  return s1_sum;
}

} // namespace

namespace primesum {

maxint_t S1(maxint_t x,
            int64_t y,
            int64_t c,
            int threads)
{
  print("");
  print("=== S1(x, y) ===");
  print("Computation of the ordinary leaves");
  print(x, y, c, threads);

  double time = get_wtime();
  maxint_t s1_sum;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
    s1_sum = S1_OpenMP_master(x, (uint32_t) y, c, threads);
  else
    s1_sum = S1_OpenMP_master(x, y, c, threads);

  print("S1_sum", s1_sum, time);
  return s1_sum;
}

} // namespace primesum
