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

#if defined(_OPENMP)
  #include <OmpLock.hpp>
#endif

namespace primesum {

class S2Status
{
public:
  S2Status(int128_t x);
  void print(int128_t n, int128_t limit);
  static double skewed_percent(int128_t n, int128_t limit);
private:
  bool is_print(double time);
  double epsilon_;
  double percent_ = -1;
  double time_ = 0;
  double is_print_ = 1.0 / 20;
  int precision_;

#if defined(_OPENMP)
  OmpLock lock_;
#endif
};

} // namespace

#endif
