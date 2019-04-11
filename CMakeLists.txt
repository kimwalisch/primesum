cmake_minimum_required(VERSION 3.4)
project(primesum CXX)
set(PRIMESUM_VERSION_MAJOR 1)
set(PRIMESUM_VERSION_MINOR 7)
set(PRIMESUM_VERSION "${PRIMESUM_VERSION_MAJOR}.${PRIMESUM_VERSION_MINOR}")

# Build options ######################################################

option(WITH_POPCNT       "Enable x86 POPCNT instruction"  ON)
option(WITH_LIBDIVIDE    "Use libdivide.h"                ON)
option(WITH_OPENMP       "Enable OpenMP support"          ON)
option(BUILD_STATIC_LIBS "Build static libprimesieve"     ON)
option(BUILD_SHARED_LIBS "Build shared libprimesieve"     OFF)
option(BUILD_PRIMESIEVE  "Build primesieve binary"        OFF)

# Set default build type to Release ##################################

if(NOT CMAKE_VERSION VERSION_LESS 3.9)
    get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
elseif(CMAKE_CONFIGURATION_TYPES)
    set(isMultiConfig TRUE)
endif()

if(NOT isMultiConfig AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# primesum binary source files #######################################

set(SRC src/app/cmdoptions.cpp
        src/app/help.cpp
        src/app/main.cpp
        src/BitSieve.cpp
        src/FactorTable.cpp
        src/generate.cpp
        src/nth_prime.cpp
        src/P2.cpp
        src/PhiTiny.cpp
        src/PiTable.cpp
        src/phi.cpp
        src/phi_sum.cpp
        src/pi_legendre.cpp
        src/primesum.cpp
        src/print.cpp
        src/S1.cpp
        src/S2LoadBalancer.cpp
        src/S2Status.cpp
        src/test.cpp
        src/Wheel.cpp
        src/deleglise-rivat/S2_trivial.cpp
        src/deleglise-rivat/S2_hard.cpp
        src/deleglise-rivat/pi_deleglise_rivat_parallel1.cpp
        src/lmo/pi_lmo1.cpp
        src/lmo/pi_lmo2.cpp
        src/lmo/pi_lmo3.cpp
        src/lmo/pi_lmo4.cpp
        src/lmo/pi_lmo5.cpp
        src/lmo/pi_lmo_parallel1.cpp)

# Include Check* #####################################################

include(CheckTypeSize)
include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
include(GNUInstallDirs)

# libprimesieve ######################################################

add_subdirectory(lib/primesieve)

# Check if int128_t type exists ######################################

check_type_size("__int128_t" __INT128_T LANGUAGE CXX)
if(NOT HAVE___INT128_T)
    check_type_size("int128_t" INT128_T LANGUAGE CXX)
    if(NOT HAVE_INT128_T)
        message(FATAL_ERROR "Compiler does not support int128_t!")
    endif()
endif()

# Check if libdivide.h compiles ######################################

if(WITH_LIBDIVIDE)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include")

    check_cxx_source_compiles("
        #include <libdivide.h>
        #include <stdint.h>
        int main() {
            libdivide::divider<uint64_t, libdivide::BRANCHFREE> d(3);
            uint64_t n = 1000000000;
            return (n / d == n / 3) ? 0 : 1;
        }" libdivide_branchfree)

    cmake_pop_check_state()
endif()

if(libdivide_branchfree)
    set(SRC ${SRC} src/deleglise-rivat/S2_easy_libdivide.cpp)
else()
    set(SRC ${SRC} src/deleglise-rivat/S2_easy.cpp)
endif()

# Check for OpenMP ###################################################

if(WITH_OPENMP)
    find_package(OpenMP QUIET)

    if(OpenMP_FOUND OR OpenMP_CXX_FOUND)
        cmake_push_check_state()
        set(CMAKE_REQUIRED_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_CXX11_STANDARD_COMPILE_OPTION}")
        set(CMAKE_REQUIRED_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include")

        check_cxx_source_compiles("
            #include <int128_t.hpp>
            int main() {
                primesum::int128_t sum = 0;
                #pragma omp parallel for reduction(+: sum)
                for (int i = 0; i < 100; i++)
                    sum += i;
                return 0;
            }" OpenMP)

        if(NOT OpenMP)
            find_library(LIB_ATOMIC NAMES atomic libatomic.so.1)

            if(NOT LIB_ATOMIC)
                set(LIB_ATOMIC "")
            else()
                set(CMAKE_REQUIRED_LIBRARIES "${LIB_ATOMIC}")

                check_cxx_source_compiles("
                    #include <int128_t.hpp>
                    int main() {
                        primesum::int128_t sum = 0;
                        #pragma omp parallel for reduction(+: sum)
                        for (int i = 0; i < 100; i++)
                            sum += i;
                        return 0;
                    }" OpenMP_with_libatomic)

                if(NOT OpenMP_with_libatomic)
                    set(LIB_ATOMIC "")
                endif()
            endif()
        endif()

        if(OpenMP OR OpenMP_with_libatomic)

            check_cxx_source_compiles("
                #include <int256_t.hpp>
                using namespace primesum;
                int main() {
                    int256_t sum = 0;
                    #pragma omp parallel for reduction(+: sum)
                    for (int i = 0; i < 100; i++)
                        sum += i;
                    return 0;
                }" OpenMP_declare_reduction)

            if(OpenMP_declare_reduction)
                if(TARGET OpenMP::OpenMP_CXX)
                    set(LIB_OPENMP "OpenMP::OpenMP_CXX")
                else()
                    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
                endif()
            endif()
        endif()

        cmake_pop_check_state()
    else()
        message(STATUS "Performing Test OpenMP")
        message(STATUS "Performing Test OpenMP - Failed")
    endif()
endif()

# Build binary #######################################################

add_executable(primesum "${SRC}")
target_link_libraries(primesum PRIVATE primesieve::primesieve "${LIB_OPENMP}" "${LIB_ATOMIC}")
target_include_directories(primesum PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
install(TARGETS primesum DESTINATION "${CMAKE_INSTALL_BINDIR}")

target_compile_features(primesum
PRIVATE
    cxx_constexpr
    cxx_lambdas
    cxx_range_for)

# Find POPCNT compiler flag ##########################################

if(WITH_POPCNT)
    if("${CMAKE_CXX_COMPILER}" MATCHES "icpc")
        cmake_push_check_state()
        set(CMAKE_REQUIRED_FLAGS -Werror)
        check_cxx_compiler_flag(-msse4.2 msse42)
        cmake_pop_check_state()

        if(msse42)
            set(POPCNT_FLAG "-msse4.2")
        endif()
    endif()

    if(NOT POPCNT_FLAG)
        cmake_push_check_state()
        set(CMAKE_REQUIRED_FLAGS -Werror)
        check_cxx_compiler_flag(-mpopcnt mpopcnt)
        cmake_pop_check_state()

        if(mpopcnt)
            set(POPCNT_FLAG "-mpopcnt")
        endif()
    endif()
endif()

if(POPCNT_FLAG)
    target_compile_options(primesum PRIVATE "${POPCNT_FLAG}")
endif()
