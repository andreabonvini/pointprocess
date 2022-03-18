# FetchContent is one of the modules which are part of the CMake distribution, we use it mainly to load an external
#  project such as, in this case, GoogleTest.
include(FetchContent)

FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG        main
)
FetchContent_MakeAvailable(googlebenchmark)

find_package (Eigen3 3.3 REQUIRED)
find_package(Boost REQUIRED)

# We define a single executable (RunBenchmarks.cpp) that will run all the available benchmarks.
add_executable(
        RunBenchmarks RunBenchmarks.cpp
        ../src/pointprocess/optimizers/BaseOptimizer.cpp
        ../src/pointprocess/optimizers/BaseOptimizer.h
        ../src/pointprocess/PointProcessUtils.cpp
        ../src/pointprocess/PointProcessUtils.h
)

target_link_libraries(RunBenchmarks  Eigen3::Eigen Boost::boost benchmark::benchmark)


