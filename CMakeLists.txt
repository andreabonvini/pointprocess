# Set minimum CMake's version required in order to build the project.
cmake_minimum_required(VERSION 3.19)

# Define Project name.
project(pointprocess)
# Set C++ standard to use.
set(CMAKE_CXX_STANDARD 17)  # 20 Doesn't work on Windows

add_subdirectory(src)

if(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  add_subdirectory(tests)
  add_subdirectory(benchmarking)
endif()



