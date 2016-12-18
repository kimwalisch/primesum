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

maxint_t F(maxint_t u)
{
  return (u * (u + 1)) / 2;
}

template <int SIGN>
maxint_t phi_sum_tiny(maxint_t x, int64_t a)
{
  maxint_t sum = 0;

  for (; a > 0; a--)
  {
    if (x <= primes_[a])
      return sum + SIGN;
    maxint_t x2 = fast_div(x, primes_[a]);
    sum += primes_[a] * phi_sum_tiny<-SIGN>(x2, a - 1);
  }

  sum += SIGN * F(x);

  return sum;
}

template <int SIGN>
maxint_t phi_sum(maxint_t x,
                 int64_t a,
                 vector<int>& primes)
{
  maxint_t sum = 0;

  for (; a > 0; a--)
  {
    if (x <= primes[a])
      return sum + SIGN;
    maxint_t x2 = fast_div(x, primes[a]);
    sum += primes[a] * phi_sum<-SIGN>(x2, a - 1, primes);
  }

  sum += SIGN * F(x);

  return sum;
}

}

namespace primesum {

maxint_t phi_sum(maxint_t x, int64_t a)
{
  if (x < 1)
    return 0;

  if (a < 10)
    return phi_sum_tiny<1>(x, a);

  vector<int> primes = generate_n_primes(a);
  return ::phi_sum<1>(x, a, primes);
}

} // namespace
