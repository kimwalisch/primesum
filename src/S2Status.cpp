///
/// @file  S2Status.cpp
/// @brief Print the status of S2(x, y) in percent.
///        Requires use of --status[=N] command-line flag.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2Status.hpp>
#include <primesum-internal.hpp>
#include <pmath.hpp>
#include <int128_t.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace primesum {

S2Status::S2Status(maxint_t x) :
  old_percent_(0),
  old_time_(0),
  print_threshold_(1.0 / 20),
  precision_(get_status_precision(x))
{
  precision_factor_ = ipow(10, precision_);
}

double S2Status::skewed_percent(maxint_t n, maxint_t limit) const
{
  double exp = 0.96;
  double percent = get_percent((double) n, (double) limit);
  double base = exp + percent / (101 / (1 - exp));
  double low = pow(base, 100.0);
  percent = 100 - in_between(0, 100 * (pow(base, percent) - low) / (1 - low), 100);

  return max(old_percent_, percent);
}

bool S2Status::is_print(double time) const
{
  return (time - old_time_) >= print_threshold_;
}

bool S2Status::is_print(double time, double percent) const
{
  if (old_time_ == 0)
    return true;

  if (!is_print(time))
    return false;

  int new_val = (int) (precision_factor_ * percent);
  int old_val = (int) (precision_factor_ * old_percent_);

  return new_val > old_val;
}

void S2Status::print(maxint_t n, maxint_t limit)
{
  double time = get_wtime();
  double percent = skewed_percent(n, limit);

  if (is_print(time, percent))
  {
    ostringstream status;
    ostringstream out;

    status << "Status: " << fixed << setprecision(precision_) << percent << "%";
    size_t spaces = status.str().length();
    string reset_line = "\r" + string(spaces,' ') + "\r";
    out << reset_line << status.str();
    cout << out.str() << flush;

    #pragma omp critical (s2_status)
    {
      old_percent_ = percent;
      old_time_ = time;
    }
  }
}

void S2Status::print(maxint_t n, maxint_t limit, double rsd)
{
  double time = get_wtime();

  if (is_print(time))
  {
    double percent = skewed_percent(n, limit);
    int load_balance = (int) in_between(0, 100 - rsd + 0.5, 100);

    ostringstream status;
    ostringstream out;

    status << "Status: " << fixed << setprecision(precision_) << percent << "%, ";
    status << "Load balance: " << load_balance << "%";
    size_t spaces = status.str().length() + 2;
    string reset_line = "\r" + string(spaces,' ') + "\r";
    out << reset_line << status.str();
    cout << out.str() << flush;

    old_percent_ = percent;
    old_time_ = time;
  }
}

} //namespace
