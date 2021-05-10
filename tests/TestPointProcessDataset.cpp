//
// Created by Andrea Bonvini on 21/04/21.
//

#include "gtest/gtest.h"
#include "TestPointProcessDataset.h"

TEST(tests,testToeplitz)
{
    ASSERT_EQ(testToeplitz1(),1);
    ASSERT_EQ(testToeplitz2(),1);
    ASSERT_EQ(testToeplitz3(),1);
}

TEST(tests,testPointProcessDataset)
{
    ASSERT_EQ(testPointProcessDataset_hasTheta0(),1);
    ASSERT_EQ(testPointProcessDataset(),1);
}