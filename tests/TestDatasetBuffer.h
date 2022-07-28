//
// Created by Andrea Bonvini on 02/01/22.
//

#ifndef POINTPROCESS_TESTDATASETBUFFER_H
#define POINTPROCESS_TESTDATASETBUFFER_H

#include "../src/pointprocess/DatasetBuffer.h"
#include <Eigen/Core>


#include <gtest/gtest.h>


namespace pointprocess::tests {
    namespace {
        // When multiple tests in a test suite need to share common objects and subroutines,
        //  you can put them into a test fixture class.
        class TestDatasetBufferFixture : public ::testing::Test {
        protected:
            std::vector<double> events = {0.9,1.8,2.7,3.6,4.5,5.6,6.6,7.2,8.1,9.0,9.9,10.8};
            double delta = 0.5;
            unsigned long bins_in_window = 10; // i.e. windowLength = 5 s (given that we have delta=0.5s)
            unsigned long last_event_index=4; // events[4] = 4.5, which is the last event observed in the first window.
            unsigned long bins = 22; // ceil(10.8/5)
            /* At each time step we should observe the following:
             *
             *     TIME_STEP OBSERVED_EVENTS            dataset.wt                      resetParameters                    eventHappened
             *     0         0.9,1.8,2.7,3.6,4.5         0.5   (currentTime: 5.0)       true  (it's the first time step)   false
             *     1         0.9,1.8,2.7,3.6,4.5         1.0   (currentTime: 5.5)       false (obs. events are the same)   false
             *     2         1.8,2.7,3.6,4.5,5.6         0.4   (currentTime: 6.0)       true  (obs. events have changed)   true
             *     3         1.8,2.7,3.6,4.5,5.6         0.9   (currentTime: 6.5)       false (obs. events are the same)   false
             *     4         2.7,3.6,4.5,5.6,6.6         0.4   (currentTime: 7.0)       true  (obs. events have changed)   true
             *     5         2.7,3.6,4.5,5.6,6.6,7.2     0.3   (currentTime: 7.5)       true  (obs. events have changed)   true
             *     6         3.6,4.5,5.6,6.6.7.2         0.8   (currentTime: 8.0)       true  (obs. events have changed)   false
             *     7         3.6,4.5,5.6,6.6.7.2,8.1     0.4   (currentTime: 8.5)       true  (obs. events have changed)   true
             *     8         4.5,5.6,6.6.7.2,8.1,9.0     0.0   (currentTime: 9.0)       true  (obs. events have changed)   true
             *     9         4.5,5.6,6.6.7.2,8.1,9.0     0.5   (currentTime: 9.5)       false (obs. events are the same)   false
             *     10        5.6,6.6.7.2,8.1,9.0,9.9     0.1   (currentTime: 10.0)      true  (obs. events have changed)   true
             *     11        5.6,6.6.7.2,8.1,9.0,9.9     0.6   (currentTime: 10.5)      false (obs. events are the same)   false
             *     12        6.6.7.2,8.1,9.0,9.9,10.8    0.2   (currentTime: 11.0)      true  (obs. events have changed)   true
             */
            unsigned long expectedSize = 13;
            std::vector<double> expectedCts = {5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0};
            std::vector<bool> expectedEhs = {false,false,true,false,true,true,false,true,true,false,true,false,true};
            std::vector<double> expectedWts = {0.5,1.0,0.4,0.9,0.4,0.3,0.8,0.4,0.0,0.5,0.1,0.6,0.2};
            std::vector<bool> expectedRps = {true,false,true,false,true,true,true,true,true,false,true,false,true};

            // These variables are needed just to create the setup object but aren't specific to this test.
            bool rightCensoring = true;
            bool hasTheta0 = true;
            unsigned char AR_ORDER = 3;
            unsigned int maxIter = 1000;
            WeightsProducer weightsProducer= WeightsProducer();

            // Create setup object to initialize DatasetBuffer
            pointprocess::PipelineSetup setup = pointprocess::PipelineSetup(delta,events,hasTheta0,AR_ORDER,last_event_index,bins,bins_in_window,weightsProducer);

        };

        TEST_F(TestDatasetBufferFixture, testDatasetGenerator1) {

            auto dataBuffer = DatasetBuffer(setup);
            EXPECT_EQ(dataBuffer.size(), expectedSize);

            std::vector<double> wts;
            std::vector<double> cts;
            std::vector<bool> ehs;
            std::vector<bool> rps;

            // FIXME: dataset here is a copy?
            for(auto [currentTime, eventHappened, resetParameters, dataset]: dataBuffer){
                cts.push_back(currentTime);
                ehs.push_back(eventHappened);
                rps.push_back(resetParameters);
                wts.push_back(dataset.wt);
            }

            for (int i = 0; i < dataBuffer.size(); ++i) {
                EXPECT_FLOAT_EQ(expectedCts[i], cts[i]) << "Vectors expectedCts and cts differ at index " << i;
            }
            for (int i = 0; i < dataBuffer.size(); ++i) {
                EXPECT_FLOAT_EQ(expectedEhs[i], ehs[i]) << "Vectors expectedEhs and ehs differ at index " << i;
            }
            for (int i = 0; i < dataBuffer.size(); ++i) {
                EXPECT_FLOAT_EQ(expectedWts[i], wts[i]) << "Vectors expectedWts and wts differ at index " << i;
            }
            for (int i = 0; i < dataBuffer.size(); ++i) {
                EXPECT_FLOAT_EQ(expectedRps[i], rps[i]) << "Vectors expectedRps and rps differ at index " << i;
            }

        }

    }
}

#endif //POINTPROCESS_TESTDATASETBUFFER_H
