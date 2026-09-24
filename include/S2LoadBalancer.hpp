///
/// @file  S2LoadBalancer.hpp
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2LOADBALANCER_HPP
#define S2LOADBALANCER_HPP

#include <aligned_vector.hpp>
#include <int128_t.hpp>
#include <S2Status.hpp>

#include <stdint.h>

namespace primesum {

class S2LoadBalancer
{
public:
  S2LoadBalancer(int128_t x, int64_t y, int64_t z, int64_t threads,
                 bool is_print = false);
  int64_t get_min_segment_size() const;
  double get_rsd() const;
  void update(int64_t low,
              int64_t threads,
              int64_t* segment_size,
              int64_t* segments_per_thread,
              const aligned_vector<double>& timings);
private:
  void init(int128_t x, int64_t y, int64_t threads);
  void set_min_size(int64_t z);
  void update(int64_t* segments_per_thread, double seconds, double pivot);
  void update_min_size(double divisor);
  double get_avg_seconds() const;
  double get_pivot(double seconds) const;
  bool is_increase(double seconds, double pivot) const;
  bool is_decrease(double seconds, double pivot) const;
  int64_t y_;
  int64_t z_;
  double rsd_;
  double count_;
  double total_seconds_;
  double min_seconds_;
  double decrease_dividend_;
  int64_t min_size_;
  int64_t sqrtz_;
  int64_t smallest_hard_leaf_;
  bool is_print_;
  S2Status status_;
};

} // namespace

#endif
