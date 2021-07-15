//
// Created by Andrea Bonvini on 21/04/21.
//

#include "gtest/gtest.h"
#include "TestWeightsProducer.h"


TEST(tests,testAlphaZero)
{
    ASSERT_EQ(testAlphaZero(),1);
}
TEST(tests,testAlphaDot9)
{
    ASSERT_EQ(testAlphaDot9(),1);
}