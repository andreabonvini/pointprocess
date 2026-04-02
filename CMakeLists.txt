cmake_minimum_required(VERSION 3.25)
project(pointprocess LANGUAGES CXX)

# On macOS, ensure the SDK sysroot is found (needed for non-Apple clang, e.g. Homebrew LLVM)
if(APPLE AND NOT CMAKE_OSX_SYSROOT)
    execute_process(
        COMMAND xcrun --show-sdk-path
        OUTPUT_VARIABLE CMAKE_OSX_SYSROOT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable compilation database for clang-tidy / LSP
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Compiler warnings ---
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wunused)
endif()

# --- Optimization flags (non-MSVC) ---
if(NOT MSVC)
    add_compile_options(
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-fno-math-errno>
        $<$<CONFIG:Release>:-ffast-math>
        $<$<CONFIG:RelWithDebInfo>:-Ofast>
    )
endif()
add_compile_definitions(
    $<$<CONFIG:Release>:EIGEN_NO_DEBUG>
    $<$<CONFIG:RelWithDebInfo>:EIGEN_NO_DEBUG>
)

# --- Optional coverage instrumentation ---
option(ENABLE_COVERAGE "Enable code coverage instrumentation (GCC/Clang)" OFF)
if(ENABLE_COVERAGE AND NOT MSVC)
    add_compile_options(-fprofile-arcs -ftest-coverage)
    add_link_options(-fprofile-arcs -ftest-coverage)
endif()

# --- Build type selection ---
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
    add_subdirectory(src)
elseif(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    add_subdirectory(tests)
endif()



