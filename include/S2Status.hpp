///
/// @file  S2Status.hpp
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2STATUS_HPP
#define S2STATUS_HPP

#include <int128_t.hpp>
#include <stdint.h>

namespace primesum {

class S2Status
{
public:
  S2Status(int128_t x, int64_t y = 0);
  void print(int128_t n, int128_t limit);
  void print_S2_hard(int64_t low, int64_t limit);
  static double skewed_percent(int128_t x, int128_t y);
private:
  double get_percent_hard(int64_t low, int64_t limit) const;
  bool is_print(double time);
  double epsilon_;
  double percent_ = -1;
  double time_ = 0;
  double is_print_ = 1.0 / 20;
  double x_tune_ = 0;
  int64_t y_log_y_ = 0;
  int precision_;
};

} // namespace

#endif
