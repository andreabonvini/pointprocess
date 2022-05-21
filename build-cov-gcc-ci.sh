#!/bin/bash

# Rationale: https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euxo pipefail

# BASE_DIR is the project's directory, containing the src/ and tests/ folders.
BASE_DIR=$PWD
COVERAGE_FILE=coverage.info


rm -rf build
mkdir build && cd build

# Configure
cmake -DCODE_COVERAGE=ON -DCMAKE_BUILD_TYPE=Release ..

# Build (for Make on Unix equivalent to `make -j $(nproc)`)
cmake --build . --config Release

# Clean-up for any previous run.
rm -f $COVERAGE_FILE
lcov --zerocounters --directory .
# Run tests
./tests/RunTests
# Create coverage report by taking into account only the files contained in src/
lcov --capture --directory tests/ -o $COVERAGE_FILE --include "$BASE_DIR/src/pointprocess/*"
