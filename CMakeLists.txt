# Set minimum CMake's version required in order to build the project.
cmake_minimum_required(VERSION 3.19)

# Define Project name.
project(pointprocess)
# Set C++ standard to use.
set(CMAKE_CXX_STANDARD 17)  # 20 Doesn't work on Windows

if(CODE_COVERAGE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage" )
endif()


add_subdirectory(src)
add_subdirectory(tests)

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_subdirectory(benchmarking)
endif()



