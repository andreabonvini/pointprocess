#!/bin/bash

# Rationale: https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euxo pipefail

# Set up defaults for CC, CXX, GCOV_PATH,SHOW_COVERAGE
export CC="${CC:-gcc}"
export CXX="${CXX:-g++}"

# Configure
cmake --preset build-release

# Build
cmake --build --preset build-release --clean-first