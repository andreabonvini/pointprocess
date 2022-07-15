# FetchContent is one of the modules which are part of the CMake distribution, we use it mainly to load an external
#  project such as, in this case, GoogleTest.
include(FetchContent)

FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG        main
)
FetchContent_MakeAvailable(googlebenchmark)

find_package(Boost REQUIRED)
find_package (Eigen3 3.3 REQUIRED)

IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

# We define a single executable (RunBenchmarks.cpp) that will run all the available benchmarks.
add_executable(
        RunBenchmarks RunBenchmarks.cpp
        ../src/pointprocess/optimizers/BaseOptimizer.cpp
        ../src/pointprocess/optimizers/BaseOptimizer.h
        ../src/pointprocess/optimizers/GaussianOptimizer.cpp
        ../src/pointprocess/optimizers/GaussianOptimizer.h
        ../src/pointprocess/optimizers/InverseGaussianOptimizer.cpp
        ../src/pointprocess/optimizers/InverseGaussianOptimizer.h
        ../src/pointprocess/optimizers/LogNormalOptimizer.cpp
        ../src/pointprocess/optimizers/LogNormalOptimizer.h
        ../src/pointprocess/PointProcessUtils.cpp
        ../src/pointprocess/PointProcessUtils.h
        ../src/pointprocess/WeightsProducer.cpp
        ../src/pointprocess/WeightsProducer.h
        ../src/pointprocess/PointProcessDataset.cpp
        ../src/pointprocess/PointProcessDataset.h
        ../src/pointprocess/OptimizersFactory.cpp
        ../src/pointprocess/OptimizersFactory.h
        ../src/pointprocess/DatasetBuffer.cpp
        ../src/pointprocess/DatasetBuffer.h
        ../src/pointprocess/spectral/spectral.cpp
        ../src/pointprocess/spectral/spectral.h
        ../src/external/indicators.h
)

if (WIN32)
    target_link_libraries(RunBenchmarks Boost::boost benchmark::benchmark)
else()
    target_link_libraries(RunBenchmarks Eigen3::Eigen Boost::boost benchmark::benchmark)
endif()



