///
/// @file   primesum-internal.hpp
/// @brief  primesum internal function definitions.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMESUM_INTERNAL_HPP
#define PRIMESUM_INTERNAL_HPP

#include <int128_t.hpp>
#include <aligned_vector.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>
#include <string>
#include <vector>

namespace primesum {

#pragma omp declare reduction(+ : int256_t : omp_out += omp_in)
#pragma omp declare reduction(+ : uint256_t : omp_out += omp_in)

enum {
  /// Uses all CPU cores.
  MAX_THREADS = -1
};

/// Silence unused parameter compiler warning
template<class T>
void unused_param(const T&)
{ }

std::string pi(const std::string& x, int threads);

/// Calculate the number of primes below x using Legendre's formula.
/// Run time: O(x) operations, O(x^(1/2)) space.
///
int64_t pi_legendre(int64_t x);

/// Calculate the number of primes below x using an optimized 
/// segmented sieve of Eratosthenes implementation.
/// Run time: O(x log log x) operations, O(x^(1/2)) space.
///
int64_t pi_primesieve(int64_t x);

/// Calculate the nth prime using a combination of the counting
/// function and the sieve of Eratosthenes.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/2)) space.
/// @pre n <= 216289611853439384
///
int64_t nth_prime(int64_t n);

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a);

maxint_t pi(int128_t x, int threads);

maxint_t pi(int128_t x);

maxint_t pi_deleglise_rivat(int128_t x, int threads);

maxint_t pi_deleglise_rivat_parallel1(int128_t x, int threads);

int64_t pi_legendre(int64_t x, int threads);

int64_t pi_lehmer(int64_t x, int threads);

int64_t pi_lehmer2(int64_t x, int threads);

maxint_t pi_lmo(int128_t x, int threads);

maxint_t pi_lmo1(int128_t x);

int64_t pi_lmo2(int64_t x);

int64_t pi_lmo3(int64_t x);

int64_t pi_lmo4(int64_t x);

int64_t pi_lmo5(int64_t x);

maxint_t pi_lmo_parallel1(int128_t x, int threads);

maxint_t pi_lmo_parallel1(int128_t x, int threads);

int64_t pi_meissel(int64_t x, int threads);

int64_t pi_primesieve(int64_t x, int threads);

int64_t phi(int64_t x, int64_t a, int threads);

maxint_t phi_sum(int128_t x, int64_t a);

int64_t nth_prime(int64_t n, int threads);

maxint_t P2(int128_t x, int64_t y, int threads);

int128_t prime_sum_tiny(int64_t x);

void set_status_precision(int precision);

int get_status_precision(int128_t x);

void set_alpha(double alpha);

double get_alpha();

double get_alpha(int128_t x, int64_t y);

double get_alpha_lmo(int128_t x);

double get_alpha_deleglise_rivat(int128_t x);

double get_wtime();

int ideal_num_threads(int threads, int64_t sieve_limit, int64_t thread_threshold = 100000);

int128_t to_int128(const std::string& expr);

template <typename T>
T get_percent(T low, T limit)
{
  T percent = (T) (100.0 * low / std::max<T>(1, limit));
  return in_between(0, percent, 100);
}

bool test();

} // namespace

#endif
