#!/bin/bash

# Rationale: https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euxo pipefail

# Set up defaults for CC, CXX, GCOV_PATH,SHOW_COVERAGE
export CC="${CC:-gcc}"
export CXX="${CXX:-g++}"
export GCOV="${GCOV:-gcov}"
export SHOW_COVERAGE="${SHOW_COVERAGE:-true}"

# Record the base directory
BASE_DIR=$PWD

# Configure
cmake --preset gcc-coverage

# Build
cmake --build --preset gcc-coverage

# Enter build directory
cd build-ci

# Clean-up counters for any previous run.
lcov --zerocounters --directory .

# Run tests
./tests/RunTests

# Create coverage report by taking into account only the files contained in src/
lcov --capture --directory tests/ -o coverage.info --include "$BASE_DIR/src/*" --gcov-tool "$GCOV"


if [ "$SHOW_COVERAGE" = "true" ]; then
    # Create HTML report in the out/ directory
    genhtml coverage.info --output-directory out

    # Open HTML
    open out/index.html
fi