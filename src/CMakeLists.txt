


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_NO_DEBUG")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-math-errno")
set(CMAKE_CXX_LAGS "${CMAKE_CXX_FLAGS} -fno-math-errno")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -openmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")

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
        pointprocess SHARED library.cpp
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
        external/csv.h
        pointprocess/optimizers/GaussianOptimizer.cpp
        pointprocess/optimizers/GaussianOptimizer.h
        pointprocess/optimizers/LogNormalOptimizer.cpp
        pointprocess/optimizers/LogNormalOptimizer.h
        pointprocess/serialize.h
)

add_executable(singleRegression examples/singleRegression.cpp
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
        external/csv.h
        pointprocess/optimizers/GaussianOptimizer.cpp
        pointprocess/optimizers/GaussianOptimizer.h
        pointprocess/optimizers/LogNormalOptimizer.cpp
        pointprocess/optimizers/LogNormalOptimizer.h
        pointprocess/serialize.h)

add_executable(fullRegression examples/fullRegression.cpp
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
        external/csv.h
        pointprocess/optimizers/GaussianOptimizer.cpp
        pointprocess/optimizers/GaussianOptimizer.h
        pointprocess/optimizers/LogNormalOptimizer.cpp
        pointprocess/optimizers/LogNormalOptimizer.h
        pointprocess/serialize.h)

# ------------------------------------------- Setup dependencies---------------------------------------------
find_package(Boost REQUIRED)
target_link_libraries(pointprocess Boost::boost)
target_link_libraries(singleRegression Boost::boost)
target_link_libraries(fullRegression Boost::boost)


# set(Eigen3_DIR "C:/Eigen3/build")   # On Windows...
find_package (Eigen3 3.3 REQUIRED)
target_link_libraries (pointprocess Eigen3::Eigen)
target_link_libraries(singleRegression Eigen3::Eigen)
target_link_libraries(fullRegression Eigen3::Eigen)