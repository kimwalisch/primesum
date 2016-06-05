///
/// @file  S2_easy_libdivide.cpp
/// @brief This is an optimized version of S2_easy which uses
///        libdivide. libdivide allows to replace expensive integer
///        divides with comparatively cheap multiplication and
///        bitshifts.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primesum-internal.hpp>
#include <generate.hpp>
#include <int128.hpp>
#include <min_max.hpp>
#include <pmath.hpp>
#include <S2Status.hpp>
#include <S2.hpp>

#include <libdivide.h>
#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primesum;

namespace {

template <typename T>
bool is_libdivide(T x)
{
  return x <= numeric_limits<uint64_t>::max();
}

typedef libdivide::divider<uint64_t, libdivide::BRANCHFREE> libdivide_u64_t;

template <typename Primes>
vector<libdivide_u64_t>
libdivide_divisors(Primes& primes)
{
  return vector<libdivide_u64_t>(primes.begin(), primes.end());
}

/// Calculate the contribution of the clustered easy
/// leaves and the sparse easy leaves.
///
template <typename Primes>
maxint_t S2_easy_OpenMP(uint128_t x,
                        int64_t y,
                        int64_t z,
                        int64_t c,
                        Primes& primes,
                        int threads)
{
  maxint_t s2_easy = 0;
  int64_t x13 = iroot<3>(x);
  threads = validate_threads(threads, x13, 1000);
  vector<libdivide_u64_t> fastdiv = libdivide_divisors(primes);

  PiTable pi(y);
  vector<int64_t> prime_sums = generate_prime_sums(y);
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

    if (is_libdivide(x2))
    {
      // Find all clustered easy leaves:
      // n = primes[b] * primes[l]
      // x / n <= y && phi(x / n, b - 1) == phi(x / m, b - 1)
      // where phi(x / n, b - 1) = pi(x / n) - b + 2
      while (l > pi_min_clustered)
      {
        int64_t xn = (uint64_t) x2 / fastdiv[l];
        int64_t phi_xn = pi[xn] - b + 2;
        maxint_t phi_xn_sum = 1 + prime_sums[pi[xn]] - prime_sums[b - 1];
        int64_t xm = (uint64_t) x2 / fastdiv[b + phi_xn - 1];
        xm = max(xm, min_clustered);
        int64_t l2 = pi[xm];
        s2_easy += prime * phi_xn_sum * (prime_sums[l] - prime_sums[l2]);
        l = l2;
      }

      // Find all sparse easy leaves:
      // n = primes[b] * primes[l]
      // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
      for (; l > pi_min_sparse; l--)
      {
        int64_t xn = (uint64_t) x2 / fastdiv[l];
        maxint_t phi = 1 + prime_sums[pi[xn]] - prime_sums[b - 1];
        s2_easy += prime * primes[l] * phi;
      }
    }
    else
    {
      // Find all clustered easy leaves:
      // n = primes[b] * primes[l]
      // x / n <= y && phi(x / n, b - 1) == phi(x / m, b - 1)
      // where phi(x / n, b - 1) = pi(x / n) - b + 2
      while (l > pi_min_clustered)
      {
        int64_t xn = (int64_t) (x2 / primes[l]);
        int64_t phi_xn = pi[xn] - b + 2;
        maxint_t phi_xn_sum = 1 + prime_sums[pi[xn]] - prime_sums[b - 1];
        int64_t xm = (int64_t) (x2 / primes[b + phi_xn - 1]);
        xm = max(xm, min_clustered);
        int64_t l2 = pi[xm];
        s2_easy += prime * phi_xn_sum * (prime_sums[l] - prime_sums[l2]);
        l = l2;
      }

      // Find all sparse easy leaves:
      // n = primes[b] * primes[l]
      // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
      for (; l > pi_min_sparse; l--)
      {
        int64_t xn = (int64_t) (x2 / primes[l]);
        maxint_t phi = 1 + prime_sums[pi[xn]] - prime_sums[b - 1];
        s2_easy += prime * primes[l] * phi;
      }
    }

    if (print_status())
      status.print(b, pi_x13);
  }

  return s2_easy;
}

} // namespace

namespace primesum {

maxint_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return S2_easy_mpi(x, y, z, c, threads);
#endif

  print("");
  print("=== S2_easy(x, y) ===");
  print("Computation of the easy special leaves");
  print(x, y, c, threads);

  double time = get_wtime();
  maxint_t s2_easy;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
  {
    vector<uint32_t> primes = generate_primes<uint32_t>(y);
    s2_easy = S2_easy_OpenMP(x, y, z, c, primes, threads);
  }
  else
  {
    vector<int64_t> primes = generate_primes<int64_t>(y);
    s2_easy = S2_easy_OpenMP(x, y, z, c, primes, threads);
  }

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

} // namespace primesum
