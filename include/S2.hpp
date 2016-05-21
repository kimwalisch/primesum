///
/// @file  S2.hpp.
/// @brief S2 function declarations.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2_HPP
#define S2_HPP

#include <int128.hpp>
#include <stdint.h>

namespace primesum {

/// ------------------------ S2_trivial() ----------------------------

maxint_t S2_trivial(maxint_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int threads);

/// ------------------------ S2_easy() -------------------------------

maxint_t S2_easy(maxint_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads);

#ifdef HAVE_MPI

maxint_t S2_easy_mpi(maxint_t x,
                     int64_t y,
                     int64_t z,
                     int64_t c,
                     int threads);

#endif

/// ------------------------ S2_hard() -------------------------------

maxint_t S2_hard(maxint_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 double alpha,
                 int threads);

#ifdef HAVE_MPI

maxint_t S2_hard_mpi(maxint_t x,
                     int64_t y,
                     int64_t z,
                     int64_t c,
                     double alpha,
                     int threads);

#endif

} // namespace primesum

#endif
