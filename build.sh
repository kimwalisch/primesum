#!/bin/sh

# Usage: ./build.sh
# Script which automates building primesum and libprimesum.
#
# What this script does:
#
# 1) Download primesieve library
# 2) Build primesieve library using: ./configure && make
# 3) Build primesum library using: ./configure && make
#
# Lots of hacks are needed because:
#
# 1) We want to build a static primesum binary.
# 2) We build primesum without first installing libprimesieve.

CONFIGURE_OPTIONS="$@"
CPU_CORES=$(nproc --all 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)

# Exit on any error
set -e

# Download primesieve
if [ ! -f ./primesieve-latest.tar.gz ]
then
    wget      https://dl.bintray.com/kimwalisch/primesieve/primesieve-latest.tar.gz || \
    curl -LkO https://dl.bintray.com/kimwalisch/primesieve/primesieve-latest.tar.gz
fi

# Build libprimesieve
if [ ! -f primesieve*/.libs/libprimesieve.a ]
then
    tar xvf primesieve-latest.tar.gz
    cd primesieve-*
    ./configure
    make -j$CPU_CORES
    cd ..
fi

# Generate ./configure script, requires Autotools
if [ ! -f ./configure ]
then
    ./autogen.sh
fi

# Patch ./configure script to continue even
# if libprimesieve is not installed
if [ "$(grep 'libprimesieve is missing' configure)" != "" ]
then
    sed '/libprimesieve is missing/c\
    true;
    ' configure > configure.tmp
    mv -f configure.tmp configure
    chmod +x configure
fi

# Generate Makefile using ./configure
if [ ! -f ./Makefile ]
then
    ./configure $CONFIGURE_OPTIONS CXXFLAGS="-O2 -I$(echo primesieve-*/include)"
fi

# Patch Makefile to build primesum binary which links
# statically against libprimesum and libprimesieve
if [ "$(grep libprimesieve.a Makefile)" = "" ]
then
    sed 's/-lprimesieve//g' Makefile > Makefile.tmp
    mv -f Makefile.tmp Makefile

    sed '/primesum_DEPENDENCIES = libprimesum.la/c\
    primesum_DEPENDENCIES = libprimesum.la primesieve-*/.libs/libprimesieve.a
    ' Makefile > Makefile.tmp
    mv -f Makefile.tmp Makefile

    sed '/primesum_LDADD = libprimesum.la/c\
    primesum_LDADD = $(OPENMP_CXXFLAGS) .libs/libprimesum.a primesieve-*/.libs/libprimesieve.a
    ' Makefile > Makefile.tmp
    mv -f Makefile.tmp Makefile

    chmod +x Makefile
fi

if [ "$(uname 2>/dev/null | egrep -i 'windows|cygwin|mingw|msys')" != "" ]
then
    # Windows: build only static library
    make libprimesum.la LDFLAGS="-static" -j$CPU_CORES
else
    # Other OSes: build static and shared library
    make libprimesum.la -j$CPU_CORES
fi

# Build statically linked primesum binary
make primesum$(grep 'EXEEXT =' Makefile | cut -f3 -d' ') LDFLAGS="-static" -j$CPU_CORES
