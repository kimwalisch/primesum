///
/// @file  print.cpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <print.hpp>
#include <primesum-internal.hpp>
#include <int128_t.hpp>
#include <stdint.h>

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace {

bool print_status_ = false;

bool print_variables_ = false;

}

namespace primesum {

void set_print_status(bool print_status)
{
#ifdef HAVE_MPI
  print_status_ = print_status && is_mpi_master_proc();
#else
  print_status_ = print_status;
#endif
}

void set_print_variables(bool print_variables)
{
#ifdef HAVE_MPI
  print_variables_ = print_variables && is_mpi_master_proc();
#else
  print_variables_ = print_variables;
#endif
}

bool print_result()
{
#ifdef HAVE_MPI
  return !print_variables() && is_mpi_master_proc();
#else
  return !print_variables();
#endif
}

void print_threads(int threads)
{
#ifdef HAVE_MPI
  cout << "processes = " << mpi_num_procs() << endl;
  cout << "threads = " << mpi_num_procs() << " * " << validate_threads(threads) << endl;
#else
  cout << "threads = " << validate_threads(threads) << endl;
#endif
}

bool print_status()
{
  return print_status_;
}

bool print_variables()
{
  return print_variables_;
}

void print(const string& str)
{
  if (print_status())
    cout << str << endl;
}

void print(maxint_t x, int64_t y, int64_t z, int64_t c, double alpha, int threads)
{
  if (print_status())
  {
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << c << endl;
    cout << "alpha = " << fixed << setprecision(3) << alpha << endl;
    print_threads(threads);
  }
}

void print(maxint_t x, int64_t y, int threads)
{
  if (print_variables())
  {
    maxint_t z = x / y;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
    print_threads(threads);
    cout << endl;
  }
}

void print(maxint_t x, int64_t y, int64_t c, int threads)
{
  if (print_variables())
  {
    maxint_t z = x / y;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << c << endl;
    cout << "alpha = " << fixed << setprecision(3) << get_alpha(x, y) << endl;
    print_threads(threads);
    cout << endl;
  }
}

void print(const string& res_str, maxint_t res, double time)
{
  if (print_status())
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
