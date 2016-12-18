///
/// @file  S1.hpp
/// @brief S1 function declarations.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <int128_t.hpp>
#include <stdint.h>

namespace primesum {

maxint_t S1(maxint_t x,
            int64_t y,
            int64_t c,
            int threads);

} // namespace primesum
