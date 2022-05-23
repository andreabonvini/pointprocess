//
// Created by Andrea Bonvini on 30/05/21.
//

#ifndef POINTPROCESS_TESTPIPELINE_H
#define POINTPROCESS_TESTPIPELINE_H

#include "../src/pointprocess/InterEventDistributions.h"
#include "../src/pointprocess/RegressionPipeline.h"
#include "../src/pointprocess/WeightsProducer.h"
#include <vector>
#include <gtest/gtest.h>


namespace pointprocess::tests {
    namespace {

        TEST(TestPipeline,TestPipelineSetup){
            std::vector<double> testEvents = {0.0,0.997,2.0,3.0,4.0,4.997,6.0};
            double windowLength = 5.0;
            double delta = 0.005;
            unsigned char AR_ORDER = 2;
            bool hasTheta0 = true;
            auto setup = pp::utils::getPipelineSetup(testEvents, hasTheta0, AR_ORDER, windowLength, delta, WeightsProducer(1.00));

            EXPECT_EQ(setup.last_event_index, 5);
            EXPECT_EQ(setup.bins, 1200);
            EXPECT_EQ(setup.bins_in_window, 1000);
        }
        TEST(TestPipeline, TestComputeTaus){
            double delta = 0.005;
            std::vector<double> testEvents = {0.0,0.997,2.0,3.0,4.0,4.997,6.0,7.0,7.883,9.0,10.222};
            bool rc = false;
            bool hasTheta0 = true;
            unsigned char ar_order = 2;
            unsigned long last_event_index = 5;
            unsigned long bins = 2045;
            unsigned long bins_in_window = 1000;
            unsigned long maxIter = 1000;
            auto wp = WeightsProducer();
            auto setup = pp::PipelineSetup(delta, testEvents, hasTheta0, ar_order, last_event_index, bins, bins_in_window, wp);

            // Here I mock by hand all the hazard-rates (lambdas).
            std::vector<double> lambdas(setup.bins - setup.bins_in_window + 1,1.0);

            /*
             * In this case we have the following fully observed intervals:
             * 6.000 to  7.000 --> Tau = 1.0 * dt = 1.0 * 1.000 = 1.000
             * 7.000 to  7.883 --> Tau = 1.0 * dt = 1.0 * 0.883 = 0.883
             * 7.885 to  9.000 --> Tau = 1.0 * dt = 1.0 * 1.115 = 1.115
             * 9.000 to 10.222 --> Tau = 1.0 * dt = 1.0 * 1.222 = 1.222
             * Note that we do not include the interval between 4.997 to 6.0 since our moving window starts at 5.000
             * (the interval is not fully observed)
             * Moreover we ignore the leftmost segment of each observed interval since the hazard-rate for that segment will
             * always be 0.0 (e.g we ignore the segment between 7.885 and 7.885)
             */
            std::vector<double> expected_taus = {1.000, 0.883, 1.115, 1.222};
            std::vector<double> taus;
            pp::utils::computeTaus(taus, lambdas, setup);
            EXPECT_EQ(taus.size(), expected_taus.size());
            for(int i = 0 ; i < taus.size(); i++) {
                EXPECT_FLOAT_EQ(taus[i], expected_taus[i]);
            }
        }
    }
}

#endif //POINTPROCESS_TESTPIPELINE_H
