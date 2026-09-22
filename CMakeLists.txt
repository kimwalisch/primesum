cmake_minimum_required(VERSION 3.9...4.2)

project(primesum
        VERSION 1.8
        LANGUAGES CXX)

# Build options ######################################################

option(BUILD_PRIMESUM    "Build the primesum binary"           ON)
option(BUILD_SHARED_LIBS "Build the shared libprimesum"        OFF)
option(BUILD_STATIC_LIBS "Build the static libprimesum"        ON)
option(BUILD_PRIMESIEVE  "Build the primesieve binary"         OFF)
option(BUILD_TESTS       "Build the test programs"             OFF)

option(WITH_POPCNT       "Enable x86 POPCNT instruction"       ON)
option(WITH_OPENMP       "Enable OpenMP multi-threading"       ON)

# Option sanity checks ###############################################

if(NOT BUILD_SHARED_LIBS AND NOT BUILD_STATIC_LIBS)
    message(FATAL_ERROR "One or both of BUILD_SHARED_LIBS or BUILD_STATIC_LIBS must be set to ON!")
endif()

if(WIN32 AND MSVC AND BUILD_SHARED_LIBS AND NOT BUILD_STATIC_LIBS AND
   (BUILD_PRIMESUM OR BUILD_TESTS))
    message(FATAL_ERROR
        "On Windows with MSVC, BUILD_STATIC_LIBS must be ON when "
        "BUILD_PRIMESUM or BUILD_TESTS is ON. "
        "Either enable BUILD_STATIC_LIBS or disable BUILD_PRIMESUM and "
        "BUILD_TESTS.")
endif()

# Set default build type to Release ##################################

get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

if(NOT isMultiConfig AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

list(APPEND PRIMESUM_COMPILE_DEFINITIONS "$<$<CONFIG:Debug>:ENABLE_ASSERT>")

# Required includes ##################################################

include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
include(GNUInstallDirs)

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

# primesum binary source files #######################################

set(BIN_SRC src/app/cmdoptions.cpp
            src/app/help.cpp
            src/app/main.cpp
            src/app/test.cpp)

# primesum library source files ######################################

set(LIB_SRC src/BitSieve.cpp
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

if(WITH_LIBDIVIDE)
    list(APPEND LIB_SRC src/deleglise-rivat/S2_easy_libdivide.cpp)
else()
    list(APPEND LIB_SRC src/deleglise-rivat/S2_easy.cpp)
endif()

# libprimesieve ######################################################

set(COPY_BUILD_TESTS "${BUILD_TESTS}")
set(BUILD_TESTS OFF CACHE BOOL "Build primesieve tests" FORCE)
add_subdirectory(lib/primesieve)
set(BUILD_TESTS "${COPY_BUILD_TESTS}" CACHE BOOL "Build test programs" FORCE)

# Check if compiler supports C++11 or later ##########################

include("${PROJECT_SOURCE_DIR}/cmake/compiler_supports_cpp11.cmake")

# Check if int128_t type exists ######################################

include("${PROJECT_SOURCE_DIR}/cmake/int128_t.cmake")

# Check for OpenMP ###################################################

include("${PROJECT_SOURCE_DIR}/cmake/OpenMP.cmake")

# Find POPCNT compiler flag ##########################################

if(WITH_POPCNT)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_FLAGS -Werror)
    check_cxx_compiler_flag(-mpopcnt mpopcnt)
    cmake_pop_check_state()

    if(mpopcnt)
        set(POPCNT_FLAG "-mpopcnt")
    endif()
endif()

# libprimesum ########################################################

if(BUILD_SHARED_LIBS)
    add_library(libprimesum SHARED ${LIB_SRC})
    set_target_properties(libprimesum PROPERTIES OUTPUT_NAME primesum)
    set_target_properties(libprimesum PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
    set_target_properties(libprimesum PROPERTIES VERSION ${PROJECT_VERSION})
    set_target_properties(libprimesum PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    target_compile_definitions(libprimesum PRIVATE ${PRIMESUM_COMPILE_DEFINITIONS})
    target_link_libraries(libprimesum PRIVATE primesieve::primesieve ${PRIMESUM_LINK_LIBRARIES})

    target_compile_features(libprimesum
    PRIVATE
        cxx_constexpr
        cxx_lambdas
        cxx_range_for)

    target_include_directories(libprimesum
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

    if(WITH_POPCNT AND POPCNT_FLAG)
        target_compile_options(libprimesum PRIVATE "${POPCNT_FLAG}")
    endif()

    install(TARGETS libprimesum
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# libprimesum-static #################################################

if(BUILD_STATIC_LIBS)
    add_library(libprimesum-static STATIC ${LIB_SRC})
    set_target_properties(libprimesum-static PROPERTIES OUTPUT_NAME primesum)
    target_compile_definitions(libprimesum-static PRIVATE ${PRIMESUM_COMPILE_DEFINITIONS})
    target_link_libraries(libprimesum-static PRIVATE primesieve::primesieve ${PRIMESUM_LINK_LIBRARIES})

    if(CMAKE_IMPORT_LIBRARY_SUFFIX AND
       "${CMAKE_IMPORT_LIBRARY_SUFFIX}" STREQUAL "${CMAKE_STATIC_LIBRARY_SUFFIX}")
        set_target_properties(libprimesum-static PROPERTIES ARCHIVE_OUTPUT_NAME primesum_static)
    endif()

    target_compile_features(libprimesum-static
    PRIVATE
        cxx_constexpr
        cxx_lambdas
        cxx_range_for)

    target_include_directories(libprimesum-static
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

    if(WITH_POPCNT AND POPCNT_FLAG)
        target_compile_options(libprimesum-static PRIVATE "${POPCNT_FLAG}")
    endif()

    install(TARGETS libprimesum-static
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Shared vs. static linking ##########################################

# On Unix-like OSes we use shared linking against libprimesum by
# default. On Windows we use static linking by default because the
# library's path is not encoded into the binary.

if(WIN32 AND TARGET libprimesum-static)
    set(STATICALLY_LINK_LIBPRIMESUM ON)
elseif(NOT TARGET libprimesum)
    set(STATICALLY_LINK_LIBPRIMESUM ON)
endif()

if(STATICALLY_LINK_LIBPRIMESUM)
    add_library(primesum::primesum ALIAS libprimesum-static)
else()
    add_library(primesum::primesum ALIAS libprimesum)
endif()

function(copy_primesum_runtime_libraries target_name)
    if(WIN32 AND NOT STATICALLY_LINK_LIBPRIMESUM)
        add_custom_command(TARGET ${target_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:libprimesum>
                $<TARGET_FILE_DIR:${target_name}>)

        add_custom_command(TARGET ${target_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:libprimesieve>
                $<TARGET_FILE_DIR:${target_name}>)
    endif()
endfunction()

# primesum binary ####################################################

if(BUILD_PRIMESUM)
    add_executable(primesum ${BIN_SRC})
    target_link_libraries(primesum PRIVATE primesum::primesum primesieve::primesieve ${PRIMESUM_LINK_LIBRARIES})
    target_compile_definitions(primesum PRIVATE ${PRIMESUM_COMPILE_DEFINITIONS})
    target_compile_features(primesum PRIVATE cxx_auto_type)
    install(TARGETS primesum DESTINATION ${CMAKE_INSTALL_BINDIR})
    copy_primesum_runtime_libraries(primesum)
endif()

# Install header #####################################################

install(FILES include/primesum.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Testing ############################################################

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
