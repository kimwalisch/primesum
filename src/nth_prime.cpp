///
/// @file  nth_prime.cpp
/// @brief Find the nth prime
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum.hpp>
#include <primesum-internal.hpp>

#include <stdint.h>
#include <string>

namespace {

// primes[1] = 2, primes[2] = 3, ...
const int primes[] = { 0, 2, 3, 5, 7, 11, 13, 17, 19, 23 };

}

namespace primesum {

int64_t nth_prime(int64_t n, int threads)
{
  if (n >= 10)
    throw primesum_error("nth_prime(n): n must be < " + 10);

  if (n < 2)
    return primes[1];

  return primes[n];
}

} // namespace primesum
