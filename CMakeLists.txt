# Set minimum CMake's version required in order to build the project.
cmake_minimum_required(VERSION 3.22)

# Define Project name.
project(pointprocess)
# Set C++ standard to use.
set(CMAKE_CXX_STANDARD 17)  # 20 Doesn't work on Windows

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_subdirectory(src)
elseif(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  add_subdirectory(tests)
  # add_subdirectory(benchmarking) # TODO: revisit, it takes too much time to compile the google::benchmark target
endif()



