cmake_minimum_required(VERSION 3.9...4.2)
project(primesum CXX)
set(PRIMESUM_VERSION_MAJOR 1)
set(PRIMESUM_VERSION_MINOR 8)
set(PRIMESUM_VERSION "${PRIMESUM_VERSION_MAJOR}.${PRIMESUM_VERSION_MINOR}")

# Build options ######################################################

option(WITH_POPCNT       "Enable x86 POPCNT instruction"  ON)
option(WITH_LIBDIVIDE    "Use libdivide.h"                ON)
option(WITH_OPENMP       "Enable OpenMP support"          ON)
option(BUILD_STATIC_LIBS "Build static libprimesieve"     ON)
option(BUILD_SHARED_LIBS "Build shared libprimesieve"     OFF)
option(BUILD_PRIMESIEVE  "Build primesieve binary"        OFF)
option(BUILD_TESTS       "Build test programs"            OFF)

# Set default build type to Release ##################################

get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

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
        src/util.cpp
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

include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
include(GNUInstallDirs)

# libprimesieve ######################################################

set(COPY_BUILD_TESTS "${BUILD_TESTS}")
set(BUILD_TESTS OFF CACHE BOOL "Build primesieve tests" FORCE)
add_subdirectory(lib/primesieve)
set(BUILD_TESTS "${COPY_BUILD_TESTS}" CACHE BOOL "Build test programs" FORCE)

# Check if compiler supports C++11 or later ##########################

include("${PROJECT_SOURCE_DIR}/cmake/compiler_supports_cpp11.cmake")

# Check if int128_t type exists ######################################

include("${PROJECT_SOURCE_DIR}/cmake/int128_t.cmake")

# Enable/Disable libdivide ###########################################

# Historically, integer division has been one of the slowest
# instructions on most CPU architectures. Hence by default we
# enable libdivide which replaces expensive integer division
# instructions by a sequence of shift, add and multiply
# instructions which is often much faster.

set(DEFAULT_LIBDIVIDE ON)

# Apple Silicon CPUs first released in 2020 have very fast integer
# division instructions. For Apple Silicon CPUs on Apple OSes
# we disable libdivide to get the best performance.

if(APPLE)
    include("${PROJECT_SOURCE_DIR}/cmake/Apple_ARM64.cmake")
    if(Apple_ARM64)
        set(DEFAULT_LIBDIVIDE OFF)
    endif()
endif()

option(WITH_LIBDIVIDE "Use libdivide.h" ${DEFAULT_LIBDIVIDE})

# Check for OpenMP ###################################################

include("${PROJECT_SOURCE_DIR}/cmake/OpenMP.cmake")

# Build binary #######################################################

add_executable(primesum "${SRC}")
target_link_libraries(primesum PRIVATE primesieve::primesieve ${PRIMESUM_LINK_LIBRARIES})
target_include_directories(primesum PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_definitions(primesum PRIVATE $<$<CONFIG:Debug>:ENABLE_ASSERT>)
install(TARGETS primesum DESTINATION "${CMAKE_INSTALL_BINDIR}")

get_target_property(PRIMESIEVE_LIBRARY_TYPE primesieve::primesieve TYPE)
if(WIN32 AND PRIMESIEVE_LIBRARY_TYPE STREQUAL "SHARED_LIBRARY")
    add_custom_command(TARGET primesum POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:primesieve::primesieve>
            $<TARGET_FILE_DIR:primesum>)
endif()

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

if(WITH_POPCNT AND POPCNT_FLAG)
    target_compile_options(primesum PRIVATE "${POPCNT_FLAG}")
endif()

# Testing ############################################################

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
