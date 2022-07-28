//
// Created by Andrea Bonvini on 15/07/22.
//

#ifndef POINTPROCESS_TESTOPTIMIZERSFACTORY_H
#define POINTPROCESS_TESTOPTIMIZERSFACTORY_H

#include "../src/pointprocess/InterEventDistributions.h"
#include "../src/pointprocess/OptimizersFactory.h"
#include <Eigen/Core>
#include <gtest/gtest.h>


namespace pointprocess {
    namespace tests {
        TEST(TestOptimizersFactory, TestGaussian) {
            auto factory = OptimizersFactory();
            auto opt = factory.create(pointprocess::Distributions::Gaussian);
            EXPECT_TRUE(opt->distribution == pointprocess::Distributions::Gaussian);
        }
        TEST(TestOptimizersFactory, TestInverseGaussian) {
            auto factory = OptimizersFactory();
            auto opt = factory.create(pointprocess::Distributions::InverseGaussian);
            EXPECT_TRUE(opt->distribution == pointprocess::Distributions::InverseGaussian);
        }
        TEST(TestOptimizersFactory, TestLogNormal) {
            auto factory = OptimizersFactory();
            auto opt = factory.create(pointprocess::Distributions::LogNormal);
            EXPECT_TRUE(opt->distribution == pointprocess::Distributions::LogNormal);
        }
    }
}

#endif //POINTPROCESS_TESTOPTIMIZERSFACTORY_H
