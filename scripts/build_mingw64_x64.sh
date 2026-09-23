#!/bin/bash

# Usage: scripts/build_mingw64_x64.sh
# Builds a primesum release binary that is statically linked
# and ready for distribution.

# === Prerequisites x64 ===
# 1) Install MSYS2 x64
# 2) pacman -Syu (exit then run it again)
# 3) pacman -S --needed base-devel mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake zip unzip git wget
# 4) git clone https://github.com/kimwalisch/primesum.git
# 5) scripts/build_mingw64_x64.sh

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
cmake .. -G "Unix Makefiles" -DCMAKE_CXX_FLAGS="-ffunction-sections -fdata-sections -mpopcnt -flto -static -static-libgcc -static-libstdc++ -Wall -Wextra -pedantic -D_WIN32_WINNT=0x601" -DCMAKE_EXE_LINKER_FLAGS="-Wl,--gc-sections" -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON
make -j8
rm primesum.exe

# Remove unnecessary libraries for linking,
# keep only GCC libraries + kernel32.
sed -i 's/-lkernel32.*/-lkernel32/g' CMakeFiles/primesum.dir/linklibs.rsp
sed -i 's/\.dll\.a/\.a/g' CMakeFiles/primesum.dir/linklibs.rsp

# Verify that sed has worked correctly,
# last word should be -lkernel32.
[ "$(grep -o '[^ ]\+$' CMakeFiles/primesum.dir/linklibs.rsp)" = "-lkernel32" ] || handle_error "failed updating linklibs.rsp"

make
strip primesum.exe

# Download the latest x64 release archive as a packaging template.
LATEST_X64_URL=$(wget -qO- https://api.github.com/repos/kimwalisch/primesum/releases/latest |
                 grep -oE 'https://[^" ]+/primesum-[^" ]+-win(-x64|64)\.zip' |
                 head -n 1)
[ -n "$LATEST_X64_URL" ] || handle_error "failed finding the latest x64 release archive"

LATEST_X64_ARCHIVE=${LATEST_X64_URL##*/}
wget "$LATEST_X64_URL"
unzip "$LATEST_X64_ARCHIVE" -d primesum-$VERSION-win-x64-tmp
rm "$LATEST_X64_ARCHIVE"

echo ""
echo ""
echo "Old file size: $(ls -l --block-size=K primesum-$VERSION-win-x64-tmp/primesum.exe)"
echo "New file size: $(ls -l --block-size=K primesum.exe)"
echo ""
echo ""

mv -f primesum.exe primesum-$VERSION-win-x64-tmp
cd primesum-$VERSION-win-x64-tmp
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

zip -r ../primesum-$VERSION-win-x64.zip *
cd ..
mv primesum-$VERSION-win-x64-tmp primesum-$VERSION-win-x64

####################################################################

echo ""
echo "Release binary built successfully!"
