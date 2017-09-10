///
/// @file  primesum.hpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License.
///

#ifndef PRIMESUM_HPP
#define PRIMESUM_HPP

#include <stdexcept>
#include <string>
#include <stdint.h>

#define PRIMESUM_VERSION "1.2"
#define PRIMESUM_VERSION_MAJOR 1
#define PRIMESUM_VERSION_MINOR 2

namespace primesum {

class primesum_error : public std::runtime_error
{
public:
  primesum_error(const std::string& msg)
    : std::runtime_error(msg)
  { }
};

/// 128-bit prime summing function.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
/// @param expr  Integer arithmetic expression e.g. "1000", "10^22"
/// @pre   expr  <= get_max_x()
///
std::string pi(const std::string& x);

/// Enable/disable printing status information during computation.
void set_print_status(bool print_status);

/// Set the number of threads.
void set_num_threads(int num_threads);

/// Get the currently set number of threads.
int get_num_threads();

/// Largest integer supported by pi(const std::string& x).
/// The return type is a string as max may be a 128-bit integer
/// which is not supported by some compilers.
/// @param alpha Tuning factor used in LMO type algorithms.
/// @return for 32-bit CPUs: 2^63-1, 
///         for 64-bit CPUs: max >= 10^27
///
std::string get_max_x(double alpha = 1.0);

/// Get the primesum version number, in the form “i.j”.
std::string primesum_version();

} // namespace

#endif
