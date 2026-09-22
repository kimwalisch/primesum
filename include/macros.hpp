///
/// @file  macros.hpp
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MACROS_HPP
#define MACROS_HPP

#ifndef __has_attribute
  #define __has_attribute(x) 0
#endif

#ifndef __has_builtin
  #define __has_builtin(x) 0
#endif

#ifndef __has_include
  #define __has_include(x) 0
#endif

#if __has_attribute(always_inline)
  #define ALWAYS_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
  #define ALWAYS_INLINE inline __forceinline
#else
  #define ALWAYS_INLINE inline
#endif

#endif
