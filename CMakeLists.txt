cmake_minimum_required(VERSION 3.19)
project(pointprocess)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_NO_DEBUG")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-math-errno")
set(CMAKE_CXX_LAGS "${CMAKE_CXX_FLAGS} -fno-math-errno")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -openmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")

# Enable vectorization for Eigen library
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" _march_native_works)
check_cxx_compiler_flag("-xHost" _xhost_works)

set(_CXX_FLAGS)

if(_march_native_works)
    message(STATUS "Using processor's vector instructions (-march=native compiler flag set)")
    set(_CXX_FLAGS "-march=native")
elseif(_xhost_works)
    message(STATUS "Using processor's vector instructions (-xHost compiler flag set)")
    set(_CXX_FLAGS "-xHost")
else()
    message(STATUS "No suitable compiler flag found for vectorization")
endif()

add_library(
        pointprocess SHARED src/library.cpp
        pointprocess/InterEventDistributions.h
        pointprocess/WeightsProducer.cpp
        pointprocess/WeightsProducer.h
        pointprocess/PointProcessDataset.cpp
        pointprocess/PointProcessDataset.h
        pointprocess/RegressionPipeline.cpp
        pointprocess/RegressionPipeline.h
        pointprocess/optimizers/BaseOptimizer.cpp
        pointprocess/optimizers/BaseOptimizer.h
        pointprocess/optimizers/InverseGaussianOptimizer.cpp
        pointprocess/optimizers/InverseGaussianOptimizer.H
        ext/csv.h
        pointprocess/optimizers/GaussianOptimizer.cpp
        pointprocess/optimizers/GaussianOptimizer.h
        pointprocess/optimizers/LogNormalOptimizer.cpp
        pointprocess/optimizers/LogNormalOptimizer.h
        pointprocess/serialize.h
        src/examples/singleRegression.cpp)

# ------------------------------------------- Setup dependencies---------------------------------------------
find_package(Boost REQUIRED)
message(STATUS "Boost Found: ${Boost_FOUND}")
target_link_libraries(pointprocess Boost::boost)

find_package(OpenMP REQUIRED)
message(STATUS "OpenMP Found: ${OpenMP_FOUND}")
target_link_libraries(pointprocess OpenMP::OpenMP_CXX)

find_package (Eigen3 3.3 REQUIRED)
message(STATUS "Eigen3 Found: ${Eigen3_FOUND}")
target_link_libraries (pointprocess Eigen3::Eigen)

# --------------------------------------------- Setup testing ---------------------------------------------
# I followed the instructions at:
# https://www.jetbrains.com/help/clion/creating-google-test-run-debug-configuration-for-test.html#add-google-tests

add_subdirectory(ext/Google_tests)
