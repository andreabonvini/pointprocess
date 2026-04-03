
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

# --- Fetch pybind11 (pinned to a stable release) ---
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v3.0.3
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(
        pointprocess SHARED pointprocess.cpp
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
        pointprocess/optimizers/GaussianOptimizer.cpp
        pointprocess/optimizers/GaussianOptimizer.h
        pointprocess/optimizers/LogNormalOptimizer.cpp
        pointprocess/optimizers/LogNormalOptimizer.h
        pointprocess/OptimizersFactory.cpp
        pointprocess/OptimizersFactory.h
        pointprocess/PointProcessUtils.h pointprocess/PointProcessUtils.cpp
        pointprocess/DatasetBuffer.cpp pointprocess/DatasetBuffer.h
        pointprocess/spectral/spectral.cpp pointprocess/spectral/spectral.h
        external/indicators.h
)

# ------------------------------------------- Setup dependencies---------------------------------------------

# Boost is used header-only (Boost.Math distributions).
# Detection priority:
#   1. CMake config-mode  – works on Homebrew, vcpkg, libboost-dev >= 1.70
#   2. Legacy FindBoost   – works when BOOST_ROOT env var is set
#                           (Chocolatey on Windows sets it as a machine env var;
#                            GitHub Actions doesn't propagate that automatically,
#                            so we mirror it into the CMake cache first.)
if(DEFINED ENV{BOOST_ROOT} AND NOT DEFINED CACHE{BOOST_ROOT})
    set(BOOST_ROOT "$ENV{BOOST_ROOT}" CACHE PATH "Boost root directory (from BOOST_ROOT env var)")
endif()

find_package(Boost CONFIG QUIET)
if(NOT Boost_FOUND)
    if(POLICY CMP0167)
        cmake_policy(SET CMP0167 OLD)  # keep legacy FindBoost module available
    endif()
    find_package(Boost REQUIRED)
endif()

# Find Eigen3 (config-mode, works with both Eigen 3.x and 5.x)
find_package(Eigen3 CONFIG REQUIRED)

# Link libraries
target_include_directories(pointprocess PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(pointprocess PRIVATE Eigen3::Eigen)

