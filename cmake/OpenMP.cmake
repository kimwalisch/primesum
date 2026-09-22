# Check if OpenMP supports 128-bit integers and 256-bit reductions
# out of the box or if we have to link against libatomic.

set(PRIMESUM_WITH_OPENMP OFF)

if(NOT WITH_OPENMP)
    return()
endif()

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

find_package(OpenMP QUIET)

if(TARGET OpenMP::OpenMP_CXX)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX")
    set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

    if(NOT compiler_supports_cpp11)
        if(CMAKE_CXX11_EXTENSION_COMPILE_OPTION)
            set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_EXTENSION_COMPILE_OPTION}")
        elseif(CMAKE_CXX11_STANDARD_COMPILE_OPTION)
            set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_STANDARD_COMPILE_OPTION}")
        endif()
    endif()

    set(OpenMP_TEST_SOURCE "
        #include <int256_t.hpp>
        #include <omp.h>
        #include <stdint.h>
        #include <iostream>
        int main(int, char** argv) {
            using namespace primesum;
            uintptr_t n = (uintptr_t) argv;
            int128_t sum128 = (int128_t) n;
            int256_t sum256 = n;
            int iters = (int) n;
            #pragma omp parallel for reduction(+: sum128, sum256)
            for (int i = 0; i < iters; i++) {
                sum128 += (i / 3) * omp_get_thread_num();
                sum256 += (i / 3) * omp_get_thread_num();
            }
            std::cout << sum128 << sum256;
            return 0;
        }")

    check_cxx_source_compiles("${OpenMP_TEST_SOURCE}" OpenMP_int256)

    if(NOT OpenMP_int256)
        # First try -latomic as the compiler may find libraries
        # in directories which CMake does not search.
        set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX" "-latomic")
        check_cxx_source_compiles("${OpenMP_TEST_SOURCE}" OpenMP_int256_with_latomic)

        if(OpenMP_int256_with_latomic)
            list(APPEND PRIMESUM_LINK_LIBRARIES "-latomic")
        else()
            find_library(LIB_ATOMIC NAMES atomic atomic.so.1 libatomic.so.1)

            if(LIB_ATOMIC)
                set(CMAKE_REQUIRED_LIBRARIES "OpenMP::OpenMP_CXX" "${LIB_ATOMIC}")
                check_cxx_source_compiles("${OpenMP_TEST_SOURCE}" OpenMP_int256_with_libatomic_path)

                if(OpenMP_int256_with_libatomic_path)
                    list(APPEND PRIMESUM_LINK_LIBRARIES "${LIB_ATOMIC}")
                endif()
            endif()
        endif()
    endif()

    cmake_pop_check_state()

    if(OpenMP_int256 OR
       OpenMP_int256_with_latomic OR
       OpenMP_int256_with_libatomic_path)
        set(PRIMESUM_WITH_OPENMP ON)
        list(APPEND PRIMESUM_LINK_LIBRARIES "OpenMP::OpenMP_CXX")
    endif()
endif()

if(NOT PRIMESUM_WITH_OPENMP)
    message(WARNING "OpenMP with 128-bit and 256-bit reductions is required for multithreading in primesum!")
endif()
