///
/// @file  RelaxedAtomic.hpp
///        LLVM/OpenMP dynamic scheduling caused a severe scaling issue
///        in primesum: https://bugs.llvm.org/show_bug.cgi?id=49588
///        By default OpenMP dynamic scheduling may process iterations
///        in random order, it is likely that this causes many cache
///        misses when computing the easy special leaves in primesum.
///        It is possible to use schedule(monotonic:dynamic) to ensure
///        that iterations are processed in sequential order. This
///        seems to fix primesum's scaling issue.
///
///        However as a precaution measure we try to avoid OpenMP
///        dynamic scheduling in primesum and instead implement
///        schedule(monotonic:dynamic, 1) ourselves using relaxed
///        atomics.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef RELAXEDATOMIC_HPP
#define RELAXEDATOMIC_HPP

#include <macros.hpp>
#include <primesum-config.hpp>

#include <atomic>

namespace {

template <typename T>
class RelaxedAtomic
{
public:
  RelaxedAtomic(T n) : atomic_(n) { }
  // Postfix Increment
  T operator++(int)
  {
    return atomic_.fetch_add(1, std::memory_order_relaxed);
  }
private:
  // Use padding to avoid CPU false sharing
  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  std::atomic<T> atomic_;
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
};

} // namespace

#endif
