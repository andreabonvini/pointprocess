 # 'Google_test' is the subproject name
 project(Google_tests)

 # 'lib' is the folder with Google Test sources
 add_subdirectory(lib)
 include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

 # 'Google_Tests_run' is the target name
 add_executable(Google_Tests_run ../../../pointprocess/WeightsProducer.cpp ../../../pointprocess/WeightsProducer.h ../../../pointprocess/PointProcessDataset.cpp ../../../pointprocess/PointProcessDataset.h ../../../tests/TestWeightsProducer.h ../../../tests/TestPointProcessDataset.h ../../../tests/TestPointProcessDataset.cpp ../../../tests/TestWeightsProducer.cpp ../../../pointprocess/RegressionPipeline.cpp ../../../pointprocess/RegressionPipeline.h ../../../pointprocess/optimizers/BaseOptimizer.cpp ../../../pointprocess/optimizers/BaseOptimizer.h ../../../pointprocess/optimizers/InverseGaussianOptimizer.cpp ../../../pointprocess/optimizers/InverseGaussianOptimizer.h ../../../tests/testData.h ../../../tests/TestOptimizers.cpp ../../../tests/TestOptimizers.h ../../pointprocess/serialize.h)
 target_link_libraries(Google_Tests_run gtest gtest_main)


 find_package(Boost REQUIRED)
 message(STATUS "Boost Found: ${Boost_FOUND}")
 target_link_libraries(Google_Tests_run Boost::boost)

 find_package (Eigen3 3.3 REQUIRED)
 message(STATUS "Eigen3 Found: ${Eigen3_FOUND}")
 target_link_libraries (Google_Tests_run Eigen3::Eigen)