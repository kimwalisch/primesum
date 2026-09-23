#!/bin/bash

# Usage: scripts/build_mingw64_arm64.sh
# Builds a primesum release binary that is statically linked
# and ready for distribution.

# === Prerequisites arm64 ===
# 1) Install a trial version of both Parallels & Windows on a MacBook ARM64.
# 2) No need to purchase/register Parallels & Windows, keep using the trial version.
# 3) Install MSYS2 x64 (or arm64 if available)
# 4) Open C:/msys64/clangarm64.exe
# 5) pacman -Syu (exit then run it again)
# 6) pacman -S mingw-w64-clang-aarch64-clang mingw-w64-clang-aarch64-llvm-openmp make git zip unzip wget
# 7) git clone https://github.com/kimwalisch/primesum.git
# 8) scripts/build_mingw64_arm64.sh

# Exit if any error occurs
set -e

rm -rf build*

####################################################################

FULL_DATE=$(date +'%B %d, %Y')
YEAR=$(date +'%Y')

cd include
VERSION=$(grep "PRIMESUM_VERSION " primesum.hpp | cut -f2 -d'"')
cd ..

####################################################################

handle_error() {
    echo ""
    echo "Error: $1"
    exit 1
}

# Build primesum binary ############################################

mkdir build-release
cd build-release

mkdir build_primesieve
cd build_primesieve
clang++ -c -I../../lib/primesieve/include -I../../lib/primesieve/src \
  -O3 -flto -static -Wall -Wextra -pedantic \
  -DENABLE_MULTIARCH_ARM_SVE -DNDEBUG -D_WIN32_WINNT=0x0A00 \
  ../../lib/primesieve/src/*.cpp ../../lib/primesieve/src/arch/arm/sve.cpp

cd ..
mkdir build_primesum
cd build_primesum
clang++ -c -I../../include -I../../src -I../../lib/primesieve/include \
  -O3 -flto -fopenmp -static -Wall -Wextra -pedantic \
  -DNDEBUG -D_WIN32_WINNT=0x0A00 \
  ../../src/*.cpp ../../src/lmo/*.cpp \
  ../../src/deleglise-rivat/S2_easy_libdivide.cpp \
  ../../src/deleglise-rivat/S2_hard.cpp \
  ../../src/deleglise-rivat/S2_trivial.cpp \
  ../../src/deleglise-rivat/pi_deleglise_rivat_parallel1.cpp \
  ../../src/app/*.cpp

cd ..
clang++ -O3 -flto -fopenmp -static -Wall -Wextra -pedantic -DENABLE_MULTIARCH_ARM_SVE -DNDEBUG -D_WIN32_WINNT=0x0A00 \
  build_primesieve/*.o build_primesum/*.o -o primesum -lPsapi

strip primesum.exe

# Download the latest x64 release archive as a packaging template.
LATEST_X64_URL=$(wget -qO- https://api.github.com/repos/kimwalisch/primesum/releases/latest |
                 grep -oE 'https://[^" ]+/primesum-[^" ]+-win(-x64|64)\.zip' |
                 head -n 1)
[ -n "$LATEST_X64_URL" ] || handle_error "failed finding the latest x64 release archive"

LATEST_X64_ARCHIVE=${LATEST_X64_URL##*/}
wget "$LATEST_X64_URL"
unzip "$LATEST_X64_ARCHIVE" -d primesum-$VERSION-win-arm64-tmp
rm "$LATEST_X64_ARCHIVE"

echo ""
echo ""
echo "Old file size: $(ls -l --block-size=K primesum-$VERSION-win-arm64-tmp/primesum.exe)"
echo "New file size: $(ls -l --block-size=K primesum.exe)"
echo ""
echo ""

mv -f primesum.exe primesum-$VERSION-win-arm64-tmp
cd primesum-$VERSION-win-arm64-tmp
sed -i "1 s/.*/primesum $VERSION/" README.txt
sed -i "2 s/.*/$FULL_DATE/" README.txt
sed -i "3 s/.*/Copyright \(c\) 2016 - $YEAR, Kim Walisch\./" COPYING

# Verify sed has worked correctly
[ "$(sed -n '1p' < README.txt)" = "primesum $VERSION" ] || handle_error "failed updating README.txt"
[ "$(sed -n '2p' < README.txt)" = "$FULL_DATE" ] || handle_error "failed updating README.txt"
[ "$(sed -n '3p' < COPYING)" = "Copyright (c) 2016 - $YEAR, Kim Walisch." ] || handle_error "failed updating COPYING"

./primesum --test
echo ""
echo ""
./primesum 1e15

# Build release zip archive ########################################

zip -r ../primesum-$VERSION-win-arm64.zip *
cd ..
mv primesum-$VERSION-win-arm64-tmp primesum-$VERSION-win-arm64

####################################################################

echo ""
echo "Release binary built successfully!"
