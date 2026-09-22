///
/// @file  generate.hpp
///
/// Copyright (C) 2018-2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef GENERATE_HPP
#define GENERATE_HPP

#include <primesieve.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace primesum {

/// Generate a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
template <typename T>
Vector<T> generate_primes(int64_t max)
{
  Vector<T> primes;
  primes.push_back(0);
  primesieve::generate_primes(max, &primes);
  return primes;
}

/// Generate a vector with the prime sums <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
template <typename T>
Vector<T> generate_prime_sums(int64_t max)
{
  Vector<T> primes;
  primes.push_back(0);
  primesieve::generate_primes(max, &primes);

  for (uint64_t i = 1; i < primes.size(); i++)
    primes[i] += primes[i - 1];

  return primes;
}

/// Generate a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
//
Vector<int32_t> generate_primes(int64_t max);

/// Generate a vector with the first n primes.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
//
Vector<int32_t> generate_n_primes(int64_t n);

/// Generate a vector with Möbius function values.
Vector<int32_t> generate_moebius(int64_t max);

/// Generate a vector with the least prime
/// factors of the integers <= max.
///
Vector<int32_t> generate_lpf(int64_t max);

/// Generate a vector with the prime counts below max
/// using the sieve of Eratosthenes.
///
Vector<int32_t> generate_pi(int64_t max);

} // namespace

#endif
