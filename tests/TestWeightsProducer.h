//
// Created by Andrea Bonvini on 20/04/21.
//

#ifndef POINTPROCESS_TESTWEIGHTSPRODUCER_H
#define POINTPROCESS_TESTWEIGHTSPRODUCER_H
#include "../src/pointprocess/WeightsProducer.h"
#include <Eigen/Core>

#include <gtest/gtest.h>


namespace pointprocess::tests {
    namespace {
        // When multiple tests in a test suite need to share common objects and subroutines,
        //  you can put them into a test fixture class.
        class TestWeightsProducerFixture : public ::testing::Test {
        protected:
            std::vector<double> target_distances = {2.0, 1.5, 1.0, 0.5, 0.0};
        };

        TEST_F(TestWeightsProducerFixture, testAlphaZero) {
            auto wp = WeightsProducer(0.0);
            Eigen::VectorXd res = wp.produce(target_distances);
            Eigen::VectorXd expected(5);
            expected << 1.0, 1.0, 1.0, 1.0, 1.0;
            ASSERT_TRUE(res.isApprox(expected));
        }

        TEST_F(TestWeightsProducerFixture, testAlphaDot9) {
            auto wp = WeightsProducer(0.9);
            Eigen::VectorXd res = wp.produce(target_distances);
            Eigen::VectorXd expected(5);
            expected << 0.165299, 0.25924, 0.40657, 0.637628, 1.0;
            ASSERT_TRUE(res.isApprox(expected, 5));
        }
    }
}




#endif //POINTPROCESS_TESTWEIGHTSPRODUCER_H
