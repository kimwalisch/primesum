///
/// @file  S2Status.hpp
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
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
  S2Status(int128_t x);
  void print(int128_t n, int128_t limit);
  double skewed_percent(int128_t n, int128_t limit) const;
private:
  bool is_print(double time, double percent) const;
  double old_percent_;
  double old_time_;
  double print_threshold_;
  int precision_;
  int precision_factor_;
};

} // namespace

#endif
