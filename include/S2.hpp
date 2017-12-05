///
/// @file  S2.hpp.
/// @brief S2 function declarations.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2_HPP
#define S2_HPP

#include <int128_t.hpp>
#include <int256_t.hpp>
#include <stdint.h>

namespace primesum {

/// ------------------------ S2_trivial() ----------------------------

int256_t S2_trivial(int128_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int threads);

/// ------------------------ S2_easy() -------------------------------

int256_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads);

/// ------------------------ S2_hard() -------------------------------

int256_t S2_hard(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads);

} // namespace

#endif
