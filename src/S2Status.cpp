///
/// @file  S2Status.cpp
/// @brief Print the status of S2(x, y) in percent.
///        Requires use of --status[=N] command-line flag.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2Status.hpp>
#include <primesum-internal.hpp>
#include <imath.hpp>
#include <int128_t.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

namespace primesum {

S2Status::S2Status(int128_t x)
{
  precision_ = get_status_precision(x);
  int q = ipow(10, precision_);
  epsilon_ = 1.0 / q;
}

/// Dirty hack!
double S2Status::skewed_percent(int128_t x, int128_t y)
{
  double exp = 0.96;
  double percent = get_percent(x, y);
  double base = exp + percent / (101 / (1 - exp));
  double low = pow(base, 100.0);
  double dividend = pow(base, percent) - low;
  percent = 100 - (100 * dividend / (1 - low));

  return percent;
}

#if defined(_OPENMP)

bool S2Status::is_print(double time)
{
  TryLock lock(lock_);
  if (lock.ownsLock())
  {
    double old = time_;
    return old == 0 ||
          (time - old) >= is_print_;
  }

  return false;
}

#else

bool S2Status::is_print(double time)
{
  double old = time_;
  return old == 0 ||
        (time - old) >= is_print_;
}

#endif

void S2Status::print(int128_t n, int128_t limit)
{
  double time = get_time();

  if (is_print(time))
  {
    time_ = time;

    double percent = skewed_percent(n, limit);
    double old = percent_;

    if ((percent - old) >= epsilon_)
    {
      percent_ = percent;
      ostringstream status;
      ostringstream out;

      status << "Status: " << fixed << setprecision(precision_) << percent << "%";
      size_t spaces = status.str().length();
      string reset_line = "\r" + string(spaces,' ') + "\r";
      out << reset_line << status.str();
      cout << out.str() << flush;
    }
  }
}

} //namespace
