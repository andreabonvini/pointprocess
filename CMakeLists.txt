cmake_minimum_required(VERSION 3.19)
project(pointprocess)

set(CMAKE_CXX_STANDARD 20)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_NO_DEBUG")

add_executable(pointprocess main.cpp
        pointprocess/InterEventDistributions.h
        pointprocess/WeightsProducer.cpp
        pointprocess/WeightsProducer.h
        pointprocess/PointProcessDataset.cpp
        pointprocess/PointProcessDataset.h
        pointprocess/RegressionPipeline.cpp pointprocess/RegressionPipeline.h pointprocess/optimizers/BaseOptimizer.cpp pointprocess/optimizers/BaseOptimizer.h pointprocess/optimizers/InverseGaussianOptimizer.cpp pointprocess/optimizers/InverseGaussianOptimizer.H ext/csv.h pointprocess/optimizers/GaussianOptimizer.cpp pointprocess/optimizers/GaussianOptimizer.h pointprocess/optimizers/LogNormalOptimizer.cpp pointprocess/optimizers/LogNormalOptimizer.h pointprocess/serialize.h)

# ------------------------------------------- Setup dependencies---------------------------------------------
find_package(Boost REQUIRED)
message(STATUS "Boost Found: ${Boost_FOUND}")
target_link_libraries(pointprocess Boost::boost)

find_package (Eigen3 3.3 REQUIRED)
message(STATUS "Eigen3 Found: ${Eigen3_FOUND}")
target_link_libraries (pointprocess Eigen3::Eigen)

# --------------------------------------------- Setup testing ---------------------------------------------
# I followed the instructions at:
# https://www.jetbrains.com/help/clion/creating-google-test-run-debug-configuration-for-test.html#add-google-tests

add_subdirectory(ext/Google_tests)
