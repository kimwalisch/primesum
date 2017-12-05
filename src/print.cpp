///
/// @file  print.cpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <print.hpp>
#include <primesum-internal.hpp>
#include <int128_t.hpp>
#include <int256_t.hpp>
#include <stdint.h>

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace {

bool print_ = false;

bool print_variables_ = false;

}

namespace primesum {

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
    print_seconds(get_wtime() - time);
  }
}

void print_seconds(double seconds)
{
  cout << "Seconds: " << fixed << setprecision(3) << seconds << endl;
}

} // namespace
