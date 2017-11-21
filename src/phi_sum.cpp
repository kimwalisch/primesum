///
/// @file  phi_sum.cpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum-internal.hpp>
#include <primesum.hpp>
#include <generate.hpp>
#include <fast_div.hpp>

#include <stdint.h>

using namespace std;
using namespace primesum;

namespace {

const int primes_[] = { 0, 2, 3, 5, 7, 11, 13, 17, 19, 23 };

template <int SIGN>
maxint_t F(maxint_t u)
{
  maxint_t u1 = u;
  u1 += 1;
  u1 *= u;
  u1 >>= 1;
  return u1;
}

template <>
maxint_t F<-1>(maxint_t u)
{
  maxint_t u1 = u;
  u1 += 1;
  u1 *= u;
  u1 >>= 1;
  u1 = -u1;
  return u1;
}

template <int SIGN>
maxint_t phi_sum_tiny(int128_t x, int64_t a)
{
  maxint_t phi;
  maxint_t sum = 0;

  for (; a > 0; a--)
  {
    if (x <= primes_[a])
      return sum + SIGN;

    int128_t x2 = fast_div(x, primes_[a]);
    phi = phi_sum_tiny<-SIGN>(x2, a - 1);
    phi *= primes_[a];
    sum += phi;
  }

  sum += F<SIGN>(x);

  return sum;
}

template <int SIGN>
maxint_t phi_sum(int128_t x,
                 int64_t a,
                 vector<int>& primes)
{
  maxint_t phi;
  maxint_t sum = 0;

  for (; a > 0; a--)
  {
    if (x <= primes[a])
      return sum + SIGN;

    int128_t x2 = fast_div(x, primes[a]);
    phi = phi_sum<-SIGN>(x2, a - 1, primes);
    phi *= primes[a];
    sum += phi;
  }

  sum += F<SIGN>(x);

  return sum;
}

}

namespace primesum {

maxint_t phi_sum(int128_t x, int64_t a)
{
  if (x < 1)
    return 0;

  if (a < 10)
    return phi_sum_tiny<1>(x, a);

  vector<int> primes = generate_n_primes(a);
  return ::phi_sum<1>(x, a, primes);
}

} // namespace
