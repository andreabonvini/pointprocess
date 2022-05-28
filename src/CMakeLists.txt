
# TODO: Should the following be added to CMakePresets.json? (Probably yes)
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


# ------------------------------------------- Setup dependencies---------------------------------------------
find_package(Boost REQUIRED)
find_package (Eigen3 3.3 REQUIRED)

IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

if (WIN32)
    target_link_libraries(pointprocess Boost::boost)
else()
    target_link_libraries (pointprocess Boost::boost Eigen3::Eigen)
endif()

