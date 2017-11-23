///
/// @file  P2.cpp
/// @brief 2nd partial sieve function.
///        P2(x, y) sums the numbers <= x that have exactly 2 prime
///        factors each exceeding the a-th prime.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum-internal.hpp>
#include <primesieve.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <int256_t.hpp>
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

/// Calculate the thread sieving distance. The idea is to
/// gradually increase the thread_distance in order to
/// keep all CPU cores busy.
///
void balanceLoad(int64_t* thread_distance, 
                 int64_t low,
                 int64_t z,
                 int threads,
                 double start_time)
{
  double seconds = get_wtime() - start_time;

  int64_t min_distance = 1 << 23;
  int64_t max_distance = ceil_div(z - low, threads);

  if (seconds < 60)
    *thread_distance *= 2;
  if (seconds > 60)
    *thread_distance /= 2;

  *thread_distance = in_between(min_distance, *thread_distance, max_distance);
}

template <typename T>
T P2_OpenMP_thread(int128_t x,
                   int64_t y,
                   int64_t z,
                   int64_t thread_distance,
                   int64_t thread_num,
                   int64_t low,
                   T& prime_sum,
                   T& correct)
{
  prime_sum = 0;
  correct = 0;
  low += thread_distance * thread_num;
  z = min(low + thread_distance, z);

  int64_t sqrtx = isqrt(x);
  int64_t start = (int64_t) max(x / z, y);
  int64_t stop  = (int64_t) min(x / low, sqrtx);
  int64_t x_div_prime = 0;

  primesieve::iterator rit(stop + 1, start);
  primesieve::iterator it(low - 1, z);

  int64_t next = it.next_prime();
  int64_t prime = rit.prev_prime();
  T P2_thread = 0;

  while (prime > start && (x_div_prime = (int64_t) (x / prime)) < z)
  {
    // Sum the primes <= x / prime
    while (next <= x_div_prime)
    {
      if (next > y && 
          next <= sqrtx)
      {
        P2_thread -= prime_sum * next;
        correct -= next;
      }

      prime_sum += next;
      next = it.next_prime();
    }

    P2_thread += prime_sum * prime;
    correct += prime;
    prime = rit.prev_prime();
  }

  // Sum the primes < z
  while (next < z)
  {
    if (next > y && 
        next <= sqrtx)
    {
      P2_thread -= prime_sum * next;
      correct -= next;
    }

    prime_sum += next;
    next = it.next_prime();
  }

  return P2_thread;
}

/// P2(x, y) sums the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime, a = pi(y).
/// Space complexity: O((x / y)^(1/2)).
///
int256_t P2_OpenMP_master(int128_t x,
                          int64_t y,
                          int threads)
{
  if (x < 4)
    return 0;

  int64_t a = pi_legendre(y, threads);
  int64_t b = pi_legendre(isqrt(x), threads);

  if (a >= b)
    return 0;

  int64_t low = 2;
  int64_t z = (int64_t)(x / max(y, 1));
  int64_t min_distance = 1 << 23;
  int64_t thread_distance = min_distance;

  aligned_vector<int256_t> prime_sums(threads);
  aligned_vector<int256_t> correct(threads);

  int256_t p2 = 0;
  int256_t prime_sum = prime_sum_tiny(y);

  while (low < z)
  {
    int64_t segments = ceil_div(z - low, thread_distance);
    threads = in_between(1, threads, segments);
    double time = get_wtime();

    #pragma omp parallel for \
        num_threads(threads) reduction(+: p2)
    for (int i = 0; i < threads; i++)
      p2 += P2_OpenMP_thread(x, y, z, thread_distance, 
         i, low, prime_sums[i], correct[i]);

    for (int i = 0; i < threads; i++)
    {
      p2 += prime_sum * correct[i];
      prime_sum += prime_sums[i];
    }

    low += thread_distance * threads;
    balanceLoad(&thread_distance, low, z, threads, time);

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

int256_t P2(int128_t x, int64_t y, int threads)
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
  int256_t p2 = P2_OpenMP_master(x, y, threads);

  print("P2", p2, time);
  return p2;
}

} // namespace
