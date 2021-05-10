//
// Created by Andrea Bonvini on 07/05/21.
//

#include "gtest/gtest.h"
#include "TestOptimizers.h"


TEST(tests,testGradients)
{
    ASSERT_EQ(testInverseGaussianGradient(),1);
    ASSERT_EQ(testGaussianGradient(),1);
    ASSERT_EQ(testLogNormalGradient(),1);
}

TEST(tests,testGradientsRc)
{
    ASSERT_EQ(testInverseGaussianGradientRc(),1);
    ASSERT_EQ(testGaussianGradientRc(),1);
    ASSERT_EQ(testLogNormalGradientRc(),1);
}

TEST(tests,testHessians)
{
    ASSERT_EQ(testInverseGaussianHessian(),1);
    ASSERT_EQ(testGaussianHessian(),1);
    ASSERT_EQ(testLogNormalHessian(),1);
}