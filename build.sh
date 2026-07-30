#!/usr/bin/env bash
# Configure, build, and install ibsimu. Replaces the old:
#   autoreconf -fi && CXXFLAGS="-O3 -fopenmp" ./configure --with-mkl --prefix=... && make -j && sudo make install
#
# Override the install location with: IBSIMU_PREFIX=/some/path ./build.sh
set -euo pipefail

PREFIX="${IBSIMU_PREFIX:-/opt/ibsimu/1.0.6}"

cmake -S . -B build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DIBSIMU_WITH_MKL=ON \
      -DCMAKE_CXX_FLAGS="-O3 -fopenmp"
cmake --build build -j
sudo cmake --install build
