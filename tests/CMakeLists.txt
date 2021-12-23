
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

find_package (Eigen3 3.3 REQUIRED)
find_package(Boost REQUIRED)


# We define a single executable (RunTests.cpp) that will run all the available tests.
add_executable(
        RunTests RunTests.cpp
        ../src/pointprocess/optimizers/BaseOptimizer.cpp
        ../src/pointprocess/optimizers/BaseOptimizer.h
        ../src/pointprocess/PointProcessUtils.cpp
        ../src/pointprocess/PointProcessUtils.h
)

target_link_libraries(RunTests gtest gtest_main Eigen3::Eigen Boost::boost)
