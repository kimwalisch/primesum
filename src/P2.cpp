///
/// @file  P2.cpp
/// @brief 2nd partial sieve function.
///        P2(x, y) sums the numbers <= x that have exactly 2 prime
///        factors each exceeding the a-th prime.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum-internal.hpp>
#include <primesieve.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <min_max.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primesum;

namespace {

/// Calculate the segment_size per thread.
/// The idea is to gradually increase the segment_size per
/// thread in order to keep all CPU cores busy.
///
int64_t balanceLoad(int64_t segment_size, double start_time)
{
  double seconds = get_wtime() - start_time;
  int min_segment_size = 1 << 20;

  if (seconds < 30)
    segment_size *= 2;
  else if (segment_size > min_segment_size)
    segment_size -= segment_size / 4;

  return segment_size;
}

template <typename T>
T P2_OpenMP_thread(T x,
                   int64_t y,
                   int64_t z,
                   int64_t segment_size,
                   int64_t thread_num,
                   int64_t low,
                   T& prime_sum,
                   T& correct)
{
  prime_sum = 0;
  correct = 0;
  low += thread_num * segment_size;
  z = min(low + segment_size, z);
  T P2_thread = 0;

  int64_t sqrtx = isqrt(x);
  int64_t start = (int64_t) max(x / z, y);
  int64_t stop  = (int64_t) min(x / low, sqrtx);
  int64_t x_div_prime = 0;

  primesieve::iterator nit(low - 1, z);
  primesieve::iterator pit(stop + 1, start);

  int64_t next_prime = nit.next_prime();
  int64_t prev_prime = pit.previous_prime();

  while (prev_prime > start &&
         (x_div_prime = (int64_t) (x / prev_prime)) < z)
  {
    // Compute the sum of the primes <= x / prev_prime
    while (next_prime <= x_div_prime)
    {
      if (next_prime > y && 
          next_prime <= sqrtx)
      {
        P2_thread -= next_prime * prime_sum;
        correct -= next_prime;
      }

      prime_sum += next_prime;
      next_prime = nit.next_prime();
    }

    P2_thread += prev_prime * prime_sum;
    correct += prev_prime;
    prev_prime = pit.previous_prime();
  }

  // Compute the sum of the primes < z
  while (next_prime < z)
  {
    if (next_prime > y && 
        next_prime <= sqrtx)
    {
      P2_thread -= next_prime * prime_sum;
      correct -= next_prime;
    }

    prime_sum += next_prime;
    next_prime = nit.next_prime();
  }

  return P2_thread;
}

/// P2(x, y) sums the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime, a = pi(y).
/// Space complexity: O((x / y)^(1/2)).
///
template <typename T>
T P2_OpenMP_master(T x, int64_t y, int threads)
{
#if __cplusplus >= 201103L
  static_assert(prt::is_signed<T>::value,
                "P2(T x, ...): T must be signed integer type");
#endif

  if (x < 4)
    return 0;

  T a = pi_legendre(y, threads);
  T b = pi_legendre((int64_t) isqrt(x), threads);

  if (a >= b)
    return 0;

  int64_t low = 2;
  int64_t z = (int64_t)(x / max(y, 1));
  int64_t segment_size = 1 << 20;
  threads = validate_threads(threads, z);

  aligned_vector<T> prime_sums(threads);
  aligned_vector<T> correct(threads);

  T p2 = 0;
  T prime_sum = prime_sum_tiny(y);

  while (low < z)
  {
    int64_t segments = ceil_div(z - low, segment_size);
    threads = in_between(1, threads, segments);
    double time = get_wtime();

    #pragma omp parallel for \
        num_threads(threads) reduction(+: p2)
    for (int i = 0; i < threads; i++)
      p2 += P2_OpenMP_thread(x, y, z, segment_size, 
         i, low, prime_sums[i], correct[i]);

    for (int i = 0; i < threads; i++)
    {
      p2 += prime_sum * correct[i];
      prime_sum += prime_sums[i];
    }

    low += segment_size * threads;
    segment_size = balanceLoad(segment_size, time);

    if (print_status())
    {
      double percent = get_percent((double) low, (double) z);
      cout << "\rStatus: " << fixed << setprecision(get_status_precision(x))
           << percent << '%' << flush;
    }
  }

  return p2;
}

} // namespace

namespace primesum {

maxint_t P2(maxint_t x, int64_t y, int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return P2_mpi(x, y, threads);
#endif

  print("");
  print("=== P2(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_wtime();
  maxint_t p2 = P2_OpenMP_master(x, y, threads);

  print("P2", p2, time);
  return p2;
}

} // namespace
