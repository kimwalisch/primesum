///
/// @file  primesum-config.hpp
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMESUM_CONFIG_HPP
#define PRIMESUM_CONFIG_HPP

#ifndef MAX_CACHE_LINE_SIZE
  /// Maximum CPU cache line size in bytes (of all CPU types that
  /// will be produced over the next few decades).
  /// In order to prevent false sharing when using a mutex (or atomic
  /// variable) we need to ensure that this mutex is stored on a
  /// cache line where no other data is stored. We achieve this by
  /// adding MAX_CACHE_LINE_SIZE bytes before and after the mutex.
  #define MAX_CACHE_LINE_SIZE 512
#endif

#endif
