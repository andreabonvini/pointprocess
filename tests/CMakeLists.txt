
# FetchContent is one of the modules which are part of the CMake distribution, we use it mainly to load an external
#  project such as, in this case, GoogleTest.
include(FetchContent)

# Automatically fetch googletest from github...
# Note that the following declaration is valid today (December 4th, 2021), if CMake gives you trouble it may be 'cause
#  the GIT_TAG changed (e.g. once it was "master" instead of "main"), so just check it on the official repo...
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        main
)
FetchContent_MakeAvailable(googletest)

find_package(Boost REQUIRED)
find_package (Eigen3 3.3 REQUIRED)

IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()

include(GoogleTest)
enable_testing()

# We define a single executable (RunTests.cpp) that will run all the available tests.


add_executable(
        RunTests RunTests.cpp
        ../src/pointprocess/InterEventDistributions.h
        ../src/pointprocess/WeightsProducer.cpp
        ../src/pointprocess/WeightsProducer.h
        ../src/pointprocess/PointProcessDataset.cpp
        ../src/pointprocess/PointProcessDataset.h
        ../src/pointprocess/RegressionPipeline.cpp
        ../src/pointprocess/RegressionPipeline.h
        ../src/pointprocess/optimizers/BaseOptimizer.cpp
        ../src/pointprocess/optimizers/BaseOptimizer.h
        ../src/pointprocess/optimizers/InverseGaussianOptimizer.cpp
        ../src/pointprocess/optimizers/InverseGaussianOptimizer.h
        ../src/pointprocess/optimizers/GaussianOptimizer.cpp
        ../src/pointprocess/optimizers/GaussianOptimizer.h
        ../src/pointprocess/optimizers/LogNormalOptimizer.cpp
        ../src/pointprocess/optimizers/LogNormalOptimizer.h
        ../src/pointprocess/OptimizersFactory.cpp
        ../src/pointprocess/OptimizersFactory.h
        ../src/pointprocess/PointProcessUtils.h
        ../src/pointprocess/PointProcessUtils.cpp
        ../src/pointprocess/DatasetBuffer.cpp
        ../src/pointprocess/DatasetBuffer.h
        ../src/pointprocess/spectral/spectral.cpp
        ../src/pointprocess/spectral/spectral.h
        ../src/external/indicators.h
)


if (WIN32)
    target_link_libraries(RunTests PRIVATE gtest gtest_main Boost::boost)
    target_include_directories(RunTests PRIVATE ${EIGEN3_INCLUDE_DIR})
else()
    target_link_libraries(RunTests PRIVATE gtest gtest_main Eigen3::Eigen Boost::boost)
endif()

gtest_discover_tests(RunTests)