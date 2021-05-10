//
// Created by Andrea Bonvini on 21/04/21.
//

#include "gtest/gtest.h"
#include "TestWeightsProducer.h"


TEST(tests,testAlphaOne)
{
    ASSERT_EQ(testAlphaOne(),1);
}
TEST(tests,testAlphaDot5)
{
    ASSERT_EQ(testAlphaDot9(),1);
}