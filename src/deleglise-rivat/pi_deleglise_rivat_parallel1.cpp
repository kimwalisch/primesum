///
/// @file  pi_deleglise_rivat_parallel1.cpp
/// @brief Parallel implementation of the Deleglise-Rivat prime
///        counting algorithm. This implementation is identical to
///        pi_deleglise_rivat_parallel2(x) but uses 128-bit integers.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum.hpp>
#include <primesum-internal.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <int128_t.hpp>
#include <S1.hpp>
#include <S2.hpp>

#include <stdint.h>
#include <string>

using namespace std;
using namespace primesum;

namespace {

/// Calculate the contribution of the special leaves.
/// @pre y > 0 && c > 1
///
maxint_t S2(maxint_t x,
            int64_t y,
            int64_t z,
            int64_t c,
            double alpha,
            int threads)
{
  maxint_t s2_trivial = S2_trivial(x, y, z, c, threads);
  maxint_t s2_easy = S2_easy(x, y, z, c, threads);
  maxint_t s2_hard = S2_hard(x, y, z, c, alpha, threads);
  maxint_t s2 = s2_trivial + s2_easy + s2_hard;

  return s2;
}

} // namespace

namespace primesum {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
maxint_t pi_deleglise_rivat_parallel1(maxint_t x, int threads)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primesum_error("pi(x): x must be <= " + limit);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  print("");
  print("=== pi_deleglise_rivat_parallel1(x) ===");
  print("pi(x) = S1 + S2 + pi(y) - 1 - P2");
  print(x, y, z, c, alpha, threads);

  maxint_t p2 = P2(x, y, threads);
  maxint_t s1 = S1(x, y, c, threads);
  maxint_t s2 = S2(x, y, z, c, alpha, threads);
  maxint_t phi = s1 + s2;
  maxint_t sum = phi + prime_sum_tiny(y) - 1 - p2;

  return sum;
}

} // namespace
