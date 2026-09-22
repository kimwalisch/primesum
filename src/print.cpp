///
/// @file  print.cpp
///
/// Copyright (C) 2017-2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <print.hpp>
#include <primesum-internal.hpp>
#include <int128_t.hpp>
#include <int256_t.hpp>
#include <stdint.h>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace {

bool print_ = false;

bool print_variables_ = false;

}

namespace primesum {

/// The compiler supports the non standard __int128_t
/// type, but the int128_t type is missing in <stdint.h>.
/// Hence, we need to define a few int128_t functions
/// that are not supported by the C++ STL.
///
#if defined(ENABLE_INT128_TO_STRING)

std::string to_string(uint128_t n)
{
  std::string str;

  while (n > 0)
  {
    str += '0' + n % 10;
    n /= 10;
  }

  if (str.empty())
    str = "0";

  std::reverse(str.begin(), str.end());

  return str;
}

std::string to_string(int128_t n)
{
  if (n >= 0)
    return to_string((uint128_t) n);
  else
  {
    // -n causes undefined behavior for n = INT128_MIN.
    // Hence we use the defined two's complement negation: ~n + 1.
    // Casting ~n to unsigned ensures the result of the addition
    // (2^127 for INT128_MIN) is safely stored in a uint128_t
    // without signed overflow.
    uint128_t abs_n = uint128_t(~n) + 1;
    return "-" + to_string(abs_n);
  }
}

std::ostream& operator<<(std::ostream& stream, int128_t n)
{
  stream << to_string(n);
  return stream;
}

std::ostream& operator<<(std::ostream& stream, uint128_t n)
{
  stream << to_string(n);
  return stream;
}

#endif

void set_print(bool is_print)
{
  print_ = is_print;
}

void set_print_variables(bool print_variables)
{
  print_variables_ = print_variables;
}

bool print_result()
{
  return !print_variables();
}

bool is_print()
{
  return print_;
}

bool print_variables()
{
  return print_variables_;
}

void print(const string& str)
{
  if (is_print())
    cout << str << endl;
}

void print(int128_t x, int64_t y, int64_t z, int64_t c, double alpha, int threads)
{
  if (is_print())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << c << endl;
    cout << "alpha = " << fixed << setprecision(3) << alpha << endl;
    cout << "threads = " << threads << endl;
  }
}

void print(int128_t x, int64_t y, int threads)
{
  if (print_variables())
  {
    int128_t z = x / y;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
    cout << "threads = " << threads << endl;
    cout << endl;
  }
}

void print(int128_t x, int64_t y, int64_t c, int threads)
{
  if (print_variables())
  {
    int128_t z = x / y;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << c << endl;
    cout << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
    cout << "threads = " << threads << endl;
    cout << endl;
  }
}

void print(const string& res_str, int256_t res, double time)
{
  if (is_print())
  {
    cout << "\r" << string(50,' ') << "\r";
    cout << "Status: 100%" << endl;
    cout << res_str << " = " << res << endl;
    print_seconds(get_time() - time);
  }
}

void print_seconds(double seconds)
{
  cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;
}

} // namespace
