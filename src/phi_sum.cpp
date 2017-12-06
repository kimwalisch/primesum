///
/// @file  phi_sum.cpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum-internal.hpp>
#include <primesum.hpp>
#include <generate.hpp>
#include <fast_div.hpp>
#include <int128_t.hpp>
#include <int256_t.hpp>

#include <stdint.h>
#include <array>

using namespace std;
using namespace primesum;

namespace {

const array<int, 10> small_primes_ = { 0, 2, 3, 5, 7, 11, 13, 17, 19, 23 };

template <int SIGN, typename Primes>
int128_t phi_sum128(int64_t x,
                    int64_t a,
                    Primes&& primes)
{
  int128_t sum = 0;

  for (; a > 0; a--)
  {
    if (x <= primes[a])
      return sum + SIGN;

    int64_t x2 = fast_div(x, primes[a]);
    sum += phi_sum128<-SIGN>(x2, a - 1, primes) * primes[a];
  }

  int128_t n = x;
  int128_t fx = (n * (n + 1)) >> 1;
  sum += fx * SIGN;

  return sum;
}

template <int SIGN, typename Primes>
int256_t phi_sum256(int128_t x,
                    int64_t a,
                    Primes&& primes)
{
  int256_t sum = 0;

  for (; a > 0; a--)
  {
    if (x <= primes[a])
      return sum + SIGN;

    int128_t x2 = fast_div(x, primes[a]);
    int256_t phi_sum;

    if (x2 <= numeric_limits<int64_t>::max())
      phi_sum = phi_sum128<-SIGN>((int64_t) x2, a - 1, primes);
    else
      phi_sum = phi_sum256<-SIGN>(x2, a - 1, primes);

    sum += phi_sum * primes[a];
  }

  int256_t n = x;
  int256_t fx = (n * (n + 1)) >> 1;
  sum += fx * SIGN;

  return sum;
}

} // namespace

namespace primesum {

int128_t phi_sum(int64_t x, int64_t a)
{
  if (x < 1)
    return 0;

  if (a < 10)
    return phi_sum128<1>(x, a, small_primes_);
  else
    return phi_sum128<1>(x, a, generate_n_primes(a));
}

int256_t phi_sum(int128_t x, int64_t a)
{
  if (x <= numeric_limits<int64_t>::max())
    return phi_sum((int64_t) x, a);

  if (a < 10)
    return phi_sum256<1>(x, a, small_primes_);
  else
    return phi_sum256<1>(x, a, generate_n_primes(a));
}

} // namespace
