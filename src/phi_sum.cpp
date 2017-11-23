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

const int primes_[] = { 0, 2, 3, 5, 7, 11, 13, 17, 19, 23 };

int256_t F(int256_t u)
{
  return ((u * u) + u) >> 1;
}

template <int SIGN>
int256_t phi_sum_tiny(int128_t x, int64_t a)
{
  int256_t sum = 0;

  for (; a > 0; a--)
  {
    if (x <= primes_[a])
      return sum + SIGN;

    int128_t x2 = fast_div(x, primes_[a]);
    sum += phi_sum_tiny<-SIGN>(x2, a - 1) * primes_[a];
  }

  sum += F(x) * SIGN;

  return sum;
}

template <int SIGN>
int256_t phi_sum(int128_t x,
                 int64_t a,
                 vector<int>& primes)
{
  int256_t sum = 0;

  for (; a > 0; a--)
  {
    if (x <= primes[a])
      return sum + SIGN;

    int128_t x2 = fast_div(x, primes[a]);
    sum += phi_sum<-SIGN>(x2, a - 1, primes) * primes[a];
  }

  sum += F(x) * SIGN;

  return sum;
}

}

namespace primesum {

int256_t phi_sum(int128_t x, int64_t a)
{
  if (x < 1)
    return 0;

  if (a < 10)
    return phi_sum_tiny<1>(x, a);

  vector<int> primes = generate_n_primes(a);
  return ::phi_sum<1>(x, a, primes);
}

} // namespace
