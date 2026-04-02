# Set minimum CMake's version required in order to build the project.
cmake_minimum_required(VERSION 3.22)

# Define Project name.
project(pointprocess LANGUAGES CXX)

# Set C++ standard to use.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable compilation database for clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add compiler flags for warnings and best practices
if(MSVC)
    # Windows (MSVC)
    add_compile_options(/W4 /WX)
else()
    # GCC and Clang
    add_compile_options(-Wall -Wextra -Wpedantic -Wunused)
    add_compile_options($<$<CONFIG:Release>:-O3>)
endif()

# Build type selection
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_subdirectory(src)
elseif(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  add_subdirectory(tests)
  # add_subdirectory(benchmarking) # TODO: revisit, it takes too much time to compile the google::benchmark target
endif()



