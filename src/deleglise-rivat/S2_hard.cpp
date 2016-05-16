///
/// @file  S2_hard.cpp
/// @brief Calculate the contribution of the hard special leaves which
///        require use of a sieve (Deleglise-Rivat algorithm).
///        This is a parallel implementation which uses compression
///        (PiTable & FactorTable) to reduce the memory usage by
///        about 10x.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <FactorTable.hpp>
#include <primesum-internal.hpp>
#include <BitSieve.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <int128.hpp>
#include <min_max.hpp>
#include <pmath.hpp>
#include <S2.hpp>
#include <S2LoadBalancer.hpp>
#include <Wheel.hpp>

#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primesum;

namespace {

/// Cross-off the multiples of prime in the sieve array.
void cross_off(BitSieve& sieve,
               int64_t low,
               int64_t high,
               int64_t prime,
               WheelItem& w)
{
  int64_t m = w.next_multiple;
  int64_t wheel_index = w.wheel_index;

  for (; m < high; m += prime * Wheel::next_multiple_factor(&wheel_index))
    sieve.unset(m - low);

  w.set(m, wheel_index);
}

/// Compute the S2 contribution of the hard special leaves which
/// require use of a sieve. Each thread processes the interval
/// [low_thread, low_thread + segments * segment_size[
/// and the missing special leaf contributions for the interval
/// [1, low_process[ are later reconstructed and added in
/// the parent S2_hard_OpenMP_master() function.
///
template <typename T, typename FactorTable, typename Primes>
T S2_hard_OpenMP_thread(T x,
                        int64_t y,
                        int64_t z,
                        int64_t c,
                        int64_t segment_size,
                        int64_t segments_per_thread,
                        int64_t thread_num,
                        int64_t low,
                        int64_t limit,
                        FactorTable& factors,
                        PiTable& pi,
                        Primes& primes,
                        vector<T>& mu_sum,
                        vector<T>& phi)
{
  low += segment_size * segments_per_thread * thread_num;
  limit = min(low + segment_size * segments_per_thread, limit);
  int64_t max_b = pi[min3(isqrt(x / low), isqrt(z), y)];
  int64_t pi_sqrty = pi[isqrt(y)];
  T s2_hard = 0;

  if (c > max_b)
    return s2_hard;

  BitSieve sieve(segment_size);
  Wheel wheel(primes, max_b + 1, low);
  vector<int32_t> counters;
  phi.resize(max_b + 1, 0);
  mu_sum.resize(max_b + 1, 0);

  // Segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = c + 1;

    // pre-sieve the multiples of the first c primes
    sieve.pre_sieve(c, low);

    // For c + 1 <= b <= pi_sqrty
    // Find all special leaves: n = primes[b] * m
    // which satisfy: mu[m] != 0 && primes[b] < lpf[m] && low <= (x / n) < high
    for (int64_t end = min(pi_sqrty, max_b); b <= end; b++)
    {
      int64_t prime = primes[b];
      T x2 = x / prime;
      int64_t x2_div_low = min(fast_div(x2, low), y);
      int64_t x2_div_high = min(fast_div(x2, high), y);
      int64_t min_m = max(x2_div_high, y / prime);
      int64_t max_m = x2_div_low;
      int64_t start = 0;

      if (prime >= max_m)
        goto next_segment;

      factors.to_index(&min_m);
      factors.to_index(&max_m);

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (prime < factors.lpf(m))
        {
          int64_t fm = factors.get_number(m);
          int64_t xn = (int64_t) fast_div(x2, fm);
          int64_t stop = xn - low;
          phi[b] += sieve.prime_sum(start, low, stop);
          start = stop + 1;
          int64_t mu_m = factors.mu(m);
          int64_t pmul = mu_m * fm * prime;
          s2_hard -= pmul * phi[b];
          mu_sum[b] -= pmul;
        }
      }

      phi[b] += sieve.prime_sum(start, low, (high - low) - 1);
      cross_off(sieve, low, high, prime, wheel[b]);
    }

    // For pi_sqrty <= b <= pi_sqrtz
    // Find all hard special leaves: n = primes[b] * primes[l]
    // which satisfy: low <= (x / n) < high
    for (; b <= max_b; b++)
    {
      int64_t prime = primes[b];
      T x2 = x / prime;
      int64_t x2_div_low = min(fast_div(x2, low), y);
      int64_t x2_div_high = min(fast_div(x2, high), y);
      int64_t l = pi[min(x2_div_low, z / prime)];
      int64_t min_hard_leaf = max3(x2_div_high, y / prime, prime);
      int64_t start = 0;

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_hard_leaf; l--)
      {
        int64_t xn = (int64_t) fast_div(x2, primes[l]);
        int64_t stop = xn - low;
        phi[b] += sieve.prime_sum(start, low, stop);
        start = stop + 1;
        int64_t pmul = primes[l] * prime;
        s2_hard += pmul * phi[b];
        mu_sum[b] += pmul;
      }

      phi[b] += sieve.prime_sum(start, low, (high - low) - 1);
      cross_off(sieve, low, high, prime, wheel[b]);
    }

    next_segment:;
  }

  return s2_hard;
}

/// Calculate the contribution of the hard special leaves which
/// require use of a sieve (to reduce the memory usage).
/// This is a parallel implementation with advanced load balancing.
/// As most special leaves tend to be in the first segments we
/// start off with a small segment size and few segments
/// per thread, after each iteration we dynamically increase
/// the segment size and the segments per thread.
///
template <typename T, typename FactorTable, typename Primes>
T S2_hard_OpenMP_master(T x,
                        int64_t y,
                        int64_t z,
                        int64_t c,
                        Primes& primes,
                        FactorTable& factors,
                        int threads)
{
  threads = validate_threads(threads, z);

  T s2_hard = 0;
  int64_t low = 1;
  int64_t limit = z + 1;
  int64_t max_prime = z / isqrt(y);

  S2LoadBalancer loadBalancer(x, y, z, threads);
  int64_t min_segment_size = loadBalancer.get_min_segment_size();
  int64_t segment_size = min_segment_size;
  int64_t segments_per_thread = 1;

  PiTable pi(max_prime);
  vector<T> phi_total(pi[isqrt(z)] + 1, 0);

  while (low < limit)
  {
    int64_t segments = ceil_div(limit - low, segment_size);
    threads = in_between(1, threads, segments);
    segments_per_thread = in_between(1, segments_per_thread, ceil_div(segments, threads));

    aligned_vector<vector<T> > phi(threads);
    aligned_vector<vector<T> > mu_sum(threads);
    aligned_vector<double> timings(threads);

    #pragma omp parallel for num_threads(threads) reduction(+: s2_hard)
    for (int i = 0; i < threads; i++)
    {
      timings[i] = get_wtime();
      s2_hard += S2_hard_OpenMP_thread(x, y, z, c, segment_size, segments_per_thread,
          i, low, limit, factors, pi, primes, mu_sum[i], phi[i]);
      timings[i] = get_wtime() - timings[i];
    }

    // Once all threads have finished reconstruct and add the 
    // missing contribution of all special leaves. This must
    // be done in order as each thread (i) requires the sum of
    // the phi values from the previous threads.
    //
    for (int i = 0; i < threads; i++)
    {
      for (size_t j = 1; j < phi[i].size(); j++)
      {
        s2_hard += phi_total[j] * mu_sum[i][j];
        phi_total[j] += phi[i][j];
      }
    }

    low += segments_per_thread * threads * segment_size;
    loadBalancer.update(low, threads, &segment_size, &segments_per_thread, timings);
  }

  return s2_hard;
}

} // namespace

namespace primesum {

maxint_t S2_hard(maxint_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return S2_hard_mpi(x, y, z, c, threads);
#endif

  print("");
  print("=== S2_hard(x, y) ===");
  print("Computation of the hard special leaves");
  print(x, y, c, threads);

  double time = get_wtime();
  maxint_t s2_hard;

  // uses less memory
  if (y <= FactorTable<uint16_t>::max())
  {
    FactorTable<uint16_t> factors(y, threads);
    int64_t max_prime = z / isqrt(y);
    vector<uint32_t> primes = generate_primes<uint32_t>(max_prime);

    s2_hard = S2_hard_OpenMP_master((maxuint_t) x, y, z, c, primes, factors, threads);
  }
  else
  {
    FactorTable<uint32_t> factors(y, threads);
    int64_t max_prime = z / isqrt(y);
    vector<int64_t> primes = generate_primes<int64_t>(max_prime);

    s2_hard = S2_hard_OpenMP_master((maxuint_t) x, y, z, c, primes, factors, threads);
  }

  print("S2_hard", s2_hard, time);
  return s2_hard;
}

} // namespace primesum
