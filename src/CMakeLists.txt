
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

# FetchContent is one of the modules which are part of the CMake distribution, we use it mainly to load an external
#  project such as, in this case, PyBind11 (needed to create Python bindings)
include(FetchContent)

# Automatically fetch pybind11 from github...
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        master
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

# Find Boost (header-only libraries for math distributions)
find_package(Boost QUIET)
if(NOT Boost_FOUND)
    message(FATAL_ERROR 
        "Boost not found. Please install Boost:\n"
        "  macOS (Homebrew):     brew install boost\n"
        "  Ubuntu/Debian:        sudo apt-get install libboost-all-dev\n"
        "  Fedora/RHEL:          sudo dnf install boost-devel\n"
        "  Windows (Chocolatey): choco install boost-msvc-14.3\n"
        "  Or download from:     https://www.boost.org/doc/libs/release/more/getting_started/\n"
    )
endif()

# Find Eigen3
find_package(Eigen3 3.3 QUIET)
if(NOT Eigen3_FOUND)
    message(FATAL_ERROR 
        "Eigen3 not found. Please install Eigen3:\n"
        "  macOS (Homebrew):     brew install eigen\n"
        "  Ubuntu/Debian:        sudo apt-get install libeigen3-dev\n"
        "  Fedora/RHEL:          sudo dnf install eigen3-devel\n"
        "  Windows (Chocolatey): choco install eigen\n"
        "  Or download from:     https://eigen.tuxfamily.org/\n"
    )
endif()

# Link libraries
target_include_directories(pointprocess PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(pointprocess PRIVATE Eigen3::Eigen)

