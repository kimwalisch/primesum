///
/// @file   primesum.cpp
/// @brief  primesum C++ API
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum-internal.hpp>
#include <primesum.hpp>
#include <primesieve.hpp>
#include <calculator.hpp>
#include <int128_t.hpp>
#include <int256_t.hpp>
#include <imath.hpp>

#include <algorithm>
#include <ctime>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <stdint.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace {

#ifdef _OPENMP
  int threads_ = 0;
#endif

int status_precision_ = -1;

double alpha_ = -1;

}

namespace primesum {

int256_t pi(int128_t x, int threads)
{
  return pi_deleglise_rivat(x, threads);
}

int256_t pi(int128_t x)
{
  return pi(x, get_num_threads());
}

/// Alias for the fastest prime summing function in primesum.
/// @param x  integer arithmetic expression e.g. "10^12".
/// @pre   x  <= get_max_x().
///
string pi(const string& x)
{
  return pi(x, get_num_threads());
}

/// Alias for the fastest prime summing function in primesum.
/// @param x  integer arithmetic expression e.g. "10^12".
/// @pre   x  <= get_max_x().
///
string pi(const string& x, int threads)
{
  int256_t pi_x = pi(to_int128(x), threads);
  ostringstream oss;
  oss << pi_x;
  return oss.str();
}

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int256_t pi_deleglise_rivat(int128_t x, int threads)
{
  return pi_deleglise_rivat_parallel1(x, threads);
}

/// Calculate the number of primes below x using Legendre's formula.
/// Run time: O(x) operations, O(x^(1/2)) space.
///
int64_t pi_legendre(int64_t x)
{
  return pi_legendre(x, get_num_threads());
}

/// Parallel implementation of the Lagarias-Miller-Odlyzko
/// prime summing algorithm using OpenMP.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * (log x)^2) space.
///
int256_t pi_lmo(int128_t x, int threads)
{
  return pi_lmo_parallel1(x, threads);
}

/// Calculate the number of primes below x using an optimized 
/// segmented sieve of Eratosthenes implementation.
/// Run time: O(x log log x) operations, O(x^(1/2)) space.
///
int64_t pi_primesieve(int64_t x)
{
  return pi_primesieve(x, get_num_threads());
}

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a)
{
  return phi(x, a, get_num_threads());
}

int128_t prime_sum_tiny(int64_t x)
{
  int64_t prime = 0;
  int128_t prime_sum = 0;
  primesieve::iterator iter(0, x);

  while ((prime = iter.next_prime()) <= x)
    prime_sum += prime;

  return prime_sum;
}

/// Returns the largest integer that can be used with
/// pi(string x). The return type is a string as max can be a 128-bit
/// integer which is not supported by all compilers.
///
string get_max_x(double alpha)
{
  ostringstream oss;

  // primesum is limited by:
  // z < 2^62, with z = x^(2/3) / alpha
  // x^(2/3) / alpha < 2^62
  // x < (2^62 * alpha)^(3/2)

  // safety buffer: use 61 instead of 62 
  double max_x = pow(pow(2.0, 61.0) * alpha, 3.0 / 2.0);
  oss << (int128_t) max_x; 

  return oss.str();
}

/// Get the wall time in seconds.
double get_wtime()
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
#endif
}

int ideal_num_threads(int threads, int64_t sieve_limit, int64_t thread_threshold)
{
  thread_threshold = max((int64_t) 1, thread_threshold);
  threads = (int) min((int64_t) threads, sieve_limit / thread_threshold);
  threads = max(1, threads);
  return threads;
}

void set_alpha(double alpha)
{
  alpha_ = alpha;
}

double get_alpha()
{
  return alpha_;
}

double get_alpha(int128_t x, int64_t y)
{
  // y = x13 * alpha, thus alpha = y / x13
  double x13 = (double) iroot<3>(x);
  return (double) y / x13;
}

/// Calculate the Lagarias-Miller-Odlyzko alpha tuning factor.
/// alpha = a log(x)^2 + b log(x) + c
/// a, b and c are constants that should be determined empirically.
/// @see ../doc/alpha-factor-tuning.pdf
///
double get_alpha_lmo(int128_t x)
{
  double alpha = get_alpha();

  // use default alpha if no command-line alpha provided
  if (alpha < 1)
  {
    double a = 0.00156512;
    double b = -0.0261411;
    double c = 0.990948;
    double logx = log((double) x);

    alpha = a * pow(logx, 2) + b * logx + c;
  }

  return in_between(1, alpha, iroot<6>(x));
}

/// Calculate the Deleglise-Rivat alpha tuning factor.
/// alpha = a log(x)^3 + b log(x)^2 + c log(x) + d
/// a, b, c and d are constants that should be determined empirically.
/// @see ../doc/alpha-tuning-factor.pdf
///
double get_alpha_deleglise_rivat(int128_t x)
{
  double alpha = get_alpha();
  double x2 = (double) x;

  // use default alpha if no command-line alpha provided
  if (alpha < 1)
  {
    double a = 0.000356618;
    double b = 0.00263762;
    double c = -0.125227;
    double d = 1.39952;
    double logx = log(x2);

    alpha = a * pow(logx, 3) + b * pow(logx, 2) + c * logx + d;
  }

  return in_between(1, alpha, iroot<6>(x));
}

void set_num_threads(int threads)
{
#ifdef _OPENMP
  threads_ = in_between(1, threads, omp_get_max_threads());
#else
  unused_param(threads);
#endif
}

int get_num_threads()
{
#ifdef _OPENMP
  if (threads_)
    return threads_;
  else
    return max(1, omp_get_max_threads());
#else
  return 1;
#endif
}

void set_status_precision(int precision)
{
  status_precision_ = in_between(0, precision, 5);
}

int get_status_precision(int128_t x)
{
  // use default precision when no command-line precision provided
  if (status_precision_ < 0)
  {
    if (x >= 1e23)
      return 2;
    if (x >= 1e21)
      return 1;
  }

  return (status_precision_ > 0) ? status_precision_ : 0;
}

int128_t to_int128(const string& expr)
{
  int128_t n = calculator::eval<int128_t>(expr);
  return n;
}

/// Get the primesum version number, in the form “i.j”.
string primesum_version()
{
  return PRIMESUM_VERSION;
}

} // namespace
