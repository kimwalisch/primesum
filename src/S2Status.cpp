///
/// @file  S2Status.cpp
/// @brief Print the status of S2(x, y) in percent.
///        Requires use of --status[=N] command-line flag.
///
/// Copyright (C) 2018-2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S2Status.hpp>
#include <primesum-internal.hpp>
#include <int128_t.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

namespace {

double log_percent(double ratio, double factor)
{
  double percent = 100.0 * log1p(factor * ratio) / log1p(factor);
  return primesum::in_between(0.0, percent, 100.0);
}

double smoothstep(double x)
{
  x = primesum::in_between(0.0, x, 1.0);
  return x * x * (3.0 - 2.0 * x);
}

double log_percent(double ratio,
                   double early_factor,
                   double base_factor,
                   double delay,
                   double cap,
                   double cutoff)
{
  double base = log_percent(ratio, base_factor);
  double boost = log_percent(ratio, early_factor);
  boost -= delay * (1.0 - smoothstep(ratio / cutoff));
  boost = primesum::in_between(0.0, boost, 100.0);

  // Dampen the early boost after its cap without creating a plateau.
  if (boost > cap)
    boost = (boost + cap) / 2.0;

  double floor = min(500.0 * ratio, 0.5);
  return max({base, boost, floor});
}

} // namespace

namespace primesum {

S2Status::S2Status(int128_t x, int64_t y)
{
  precision_ = get_status_precision(x);
  epsilon_ = 1.0;
  for (int i = 0; i < precision_; i++)
    epsilon_ /= 10.0;

  if (y > 0)
  {
    y_log_y_ = int64_t(y * log(double(y)));
    x_tune_ = in_between(0.0, (log10(double(x)) - 20.0) / 2.0, 1.0);
  }
}

double S2Status::get_percent_hard(int64_t low, int64_t limit) const
{
  // The sieve position is most useful near completion.
  double percent1 = get_percent(low, limit);

  // The first y log(y) values cover much of the early work.
  double percent2 = get_percent(low, y_log_y_);
  percent2 = min(percent2, 30.0);

  // Primecount's estimate for the uneven distribution of hard leaves.
  double ratio = percent1 / 100.0;
  double small = log_percent(ratio, 2643.010656, 21.015052,
                             11.846115, 56.508811, 0.000203036);
  double large = log_percent(ratio, 25589.45108, 15.357592,
                             42.898382, 54.704957, 0.000411627);
  double percent3 = small * (1.0 - x_tune_) + large * x_tune_;

  return max({percent1, percent2, percent3});
}

void S2Status::print_S2_hard(int64_t low, int64_t limit)
{
  // print(res, time) writes the final 100%
  // when the formula completes.
  if (low >= limit)
    return;

  double time = get_time();

  if (is_print(time))
  {
    time_ = time;
    double percent = get_percent_hard(low, limit);

    percent = min(percent, 100.0 - epsilon_);

    if ((percent - percent_) >= epsilon_)
    {
      percent_ = percent;
      ostringstream status;
      status << "Status: " << fixed << setprecision(precision_)
             << percent << '%';
      cout << '\r' << string(status.str().length(), ' ') << '\r'
           << status.str() << flush;
    }
  }
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

bool S2Status::is_print(double time)
{
  double old = time_;
  return old == 0 ||
        (time - old) >= is_print_;
}

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
