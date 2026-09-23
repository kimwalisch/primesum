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

#ifndef __has_cpp_attribute
  #define __has_cpp_attribute(x) 0
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

#if __cplusplus >= 202002L && \
    __has_cpp_attribute(likely)
  #define if_likely(x) if (x) [[likely]]
#elif defined(__GNUC__) || \
      __has_builtin(__builtin_expect)
  #define if_likely(x) if (__builtin_expect(!!(x), 1))
#else
  #define if_likely(x) if (x)
#endif

#if __cplusplus >= 202002L && \
    __has_cpp_attribute(unlikely)
  #define if_unlikely(x) if (x) [[unlikely]]
#elif defined(__GNUC__) || \
      __has_builtin(__builtin_expect)
  #define if_unlikely(x) if (__builtin_expect(!!(x), 0))
#else
  #define if_unlikely(x) if (x)
#endif

#if __cplusplus >= 201703L && \
    __has_cpp_attribute(maybe_unused)
  #define MAYBE_UNUSED [[maybe_unused]]
#elif __has_attribute(unused)
  #define MAYBE_UNUSED __attribute__((unused))
#else
  #define MAYBE_UNUSED
#endif

// Silence GCC < 12 warning:
// warning: 'unused' attribute ignored [-Wattributes]
#if defined(__GNUC__) && \
   !defined(__clang__)
  #if __GNUC__ < 12
    #undef MAYBE_UNUSED
    #define MAYBE_UNUSED
  #endif
#endif

/// Unrolling loops that execute very few iterations on average
/// tends to deteriorate performance due to increased branch
/// mispredictions. Using the NO_UNROLL_LOOP macro we can disable
/// loop unrolling for such loops.
#if defined(__clang__)
  #define NO_UNROLL_LOOP _Pragma("nounroll")
#elif defined(__GNUC__) && __GNUC__ >= 8
  #define NO_UNROLL_LOOP _Pragma("GCC unroll 0")
#else
  #define NO_UNROLL_LOOP
#endif

/// Enable expensive debugging assertions.
/// These assertions enable e.g. bounds checks for the
/// Vector and Array types.
///
#if defined(ENABLE_ASSERT)
  namespace primesum {
  [[noreturn]]
  void assert_failed(const char* assertion,
                     const char* file,
                     unsigned int line,
                     const char* function);
  } // namespace

  #if defined(_MSC_VER)
    #define ASSERT_FUNCTION __FUNCSIG__
  #elif defined(__GNUC__) || defined(__clang__)
    #define ASSERT_FUNCTION __PRETTY_FUNCTION__
  #else
    #define ASSERT_FUNCTION __func__
  #endif

  #define ASSERT(x) \
    do { \
      if(!(x)) \
        primesum::assert_failed(#x, __FILE__, __LINE__, ASSERT_FUNCTION); \
    } while (0)
#else
  #define ASSERT(x) ((void) 0)
#endif

#endif
