



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
        pointprocess/optimizers/InverseGaussianOptimizer.h
        external/csv.h
        pointprocess/optimizers/GaussianOptimizer.cpp
        pointprocess/optimizers/GaussianOptimizer.h
        pointprocess/optimizers/LogNormalOptimizer.cpp
        pointprocess/optimizers/LogNormalOptimizer.h
        pointprocess/OptimizersFactory.h
        pointprocess/PointProcessUtils.h pointprocess/PointProcessUtils.cpp
        pointprocess/DatasetBuffer.cpp pointprocess/DatasetBuffer.h
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
        pointprocess/optimizers/InverseGaussianOptimizer.h
        external/csv.h
        pointprocess/optimizers/GaussianOptimizer.cpp
        pointprocess/optimizers/GaussianOptimizer.h
        pointprocess/optimizers/LogNormalOptimizer.cpp
        pointprocess/optimizers/LogNormalOptimizer.h
        pointprocess/OptimizersFactory.h
        pointprocess/PointProcessUtils.h pointprocess/PointProcessUtils.cpp
        pointprocess/DatasetBuffer.cpp pointprocess/DatasetBuffer.h)

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
        pointprocess/optimizers/InverseGaussianOptimizer.h
        external/csv.h
        pointprocess/optimizers/GaussianOptimizer.cpp
        pointprocess/optimizers/GaussianOptimizer.h
        pointprocess/optimizers/LogNormalOptimizer.cpp
        pointprocess/optimizers/LogNormalOptimizer.h
        pointprocess/OptimizersFactory.h
        pointprocess/PointProcessUtils.h
        pointprocess/PointProcessUtils.cpp
        pointprocess/DatasetBuffer.cpp
        pointprocess/DatasetBuffer.h)

# ------------------------------------------- Setup dependencies---------------------------------------------
find_package(Boost REQUIRED)
if(Boost_FOUND)
    message(STATUS "Boost found!")
endif()
target_link_libraries(pointprocess Boost::boost)
target_link_libraries(singleRegression Boost::boost)
target_link_libraries(fullRegression Boost::boost)


find_package (Eigen3 3.3 REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "Eigen found!")
endif()

IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

if (WIN32)
    target_link_libraries (pointprocess)
    target_link_libraries(singleRegression)
    target_link_libraries(fullRegression)
else()
    target_link_libraries (pointprocess Eigen3::Eigen)
    target_link_libraries(singleRegression Eigen3::Eigen)
    target_link_libraries(fullRegression Eigen3::Eigen)
endif()

