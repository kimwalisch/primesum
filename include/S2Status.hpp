///
/// @file  S2Status.hpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2STATUS_HPP
#define S2STATUS_HPP

#include <int128_t.hpp>

namespace primesum {

class S2Status
{
public:
  S2Status(maxint_t x);
  void print(maxint_t n, maxint_t limit);
  void print(maxint_t n, maxint_t limit, double rsd);
  double skewed_percent(maxint_t n, maxint_t limit) const;
private:
  bool is_print(double time) const;
  bool is_print(double time, double percent) const;
  double old_percent_;
  double old_time_;
  double print_threshold_;
  int precision_;
  int precision_factor_;
};

} // namespace

#endif
