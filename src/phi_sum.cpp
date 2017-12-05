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

using namespace std;
using namespace primesum;

namespace {

const int small_primes_[] = { 0, 2, 3, 5, 7, 11, 13, 17, 19, 23 };

template <int SIGN, typename T, typename Primes>
typename next_larger_type<T>::type
phi_sum(T x,
        int64_t a,
        Primes& primes)
{
  using res_t = typename next_larger_type<T>::type;

  res_t sum = 0;

  for (; a > 0; a--)
  {
    if (x <= primes[a])
      return sum + SIGN;

    T x2 = fast_div(x, primes[a]);
    sum += phi_sum<-SIGN>(x2, a - 1, primes) * primes[a];
  }

  res_t n = x;
  res_t fx = (n * (n + 1)) >> 1;
  sum += fx * SIGN;

  return sum;
}

}

namespace primesum {

int128_t phi_sum(int64_t x, int64_t a)
{
  if (x < 1)
    return 0;

  if (a < 10)
    return ::phi_sum<1>(x, a, small_primes_);
  else
  {
    auto primes = generate_n_primes(a);
    return ::phi_sum<1>(x, a, primes);
  }
}

int256_t phi_sum(int128_t x, int64_t a)
{
  if (x < 1)
    return 0;

  if (a < 10)
    return ::phi_sum<1>(x, a, small_primes_);
  else
  {
    auto primes = generate_n_primes(a);
    return ::phi_sum<1>(x, a, primes);
  }
}

} // namespace
