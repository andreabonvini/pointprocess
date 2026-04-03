
# --- Fetch GoogleTest (pinned to a stable release) ---
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.17.0
)
FetchContent_MakeAvailable(googletest)

# --- Dependencies ---
# Boost is used header-only (Boost.Math distributions).
# Try CMake config-mode first (Homebrew, vcpkg, libboost-dev >= 1.70);
# fall back to the legacy FindBoost module for setups that lack BoostConfig.cmake
# (e.g. Chocolatey on Windows, which sets BOOST_ROOT instead).
find_package(Boost CONFIG QUIET)
if(NOT Boost_FOUND)
    if(POLICY CMP0167)
        cmake_policy(SET CMP0167 OLD)  # keep FindBoost module available
    endif()
    find_package(Boost REQUIRED)
endif()
find_package(Eigen3 CONFIG REQUIRED)

include(GoogleTest)
enable_testing()

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

target_link_libraries(RunTests PRIVATE gtest gtest_main Eigen3::Eigen Boost::boost)

gtest_discover_tests(RunTests)