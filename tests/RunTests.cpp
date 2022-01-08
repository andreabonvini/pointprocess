//
// Created by Andrea Bonvini on 04/12/21.
//

// A good resource to understand how you should structure your C++ tests with GoogleTest is
//   https://github.com/google/googletest/blob/main/docs/primer.md

#include <gtest/gtest.h>

#include "TestWeightsProducer.h"
#include "TestPointProcessDataset.h"
#include "TestOptimizers.h"
#include "TestPipeline.h"
#include "TestDatasetBuffer.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}