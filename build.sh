#!/bin/bash

# Rationale: https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euxo pipefail

# The -p option allows the 'mkdir' command not to return any error in case that the directory already exists.
mkdir -p build && cd build

# Configure
cmake -DCODE_COVERAGE=ON -DCMAKE_BUILD_TYPE=Release ..
# Build (for Make on Unix equivalent to `make -j $(nproc)`)
cmake --build . --config Release -- -j $(nproc)
