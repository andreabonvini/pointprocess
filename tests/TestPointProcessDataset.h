//
// Created by Andrea Bonvini on 21/04/21.
//

#ifndef POINTPROCESS_TESTPOINTPROCESSDATASET_H
#define POINTPROCESS_TESTPOINTPROCESSDATASET_H
#include "../src/pointprocess/PointProcessDataset.h"
#include <Eigen/Core>


#include <gtest/gtest.h>


namespace pointprocess::tests {
    namespace {
        // When multiple tests in a test suite need to share common objects and subroutines,
        //  you can put them into a test fixture class.
        class TestPointProcessDatasetFixture : public ::testing::Test {
        protected:
            std::deque<double> events= {
                    727.391,
                    728.297,
                    729.188,
                    730.062,
                    730.984,
                    731.93,
                    732.875,
                    733.828,
                    734.781,
                    735.711,
            };
            // note that diff(events): 0.906, 0.891, 0.874, 0.922, 0.946, 0.945, 0.953, 0.953, 0.93
            unsigned char ar_order = 3;
            unsigned char n_samples = 6;
            Eigen::VectorXd expected_wn;
            Eigen::VectorXd expected_eta;

            WeightsProducer wp = WeightsProducer(0.0);
            double current_time = 736.0;
            double expected_wt = current_time - events[events.size() - 1];

            void SetUp() override {
                expected_eta.resize(n_samples);
                expected_eta << 1.000, 1.000, 1.000, 1.000, 1.000, 1.000;
                expected_wn.resize(n_samples);
                expected_wn << 0.922, 0.946, 0.945, 0.953, 0.953, 0.930;
            }
        };

        TEST_F(TestPointProcessDatasetFixture, testPointProcessDataset_hasTheta0False) {

            bool hasTheta0 = false;

            Eigen::VectorXd expected_xt(ar_order + hasTheta0);
            Eigen::MatrixXd expected_xn(n_samples, ar_order + hasTheta0);

            expected_xt << 0.930, 0.953, 0.953;
            expected_xn << 0.874, 0.891, 0.906,
                    0.922, 0.874, 0.891,
                    0.946, 0.922, 0.874,
                    0.945, 0.946, 0.922,
                    0.953, 0.945, 0.946,
                    0.953, 0.953, 0.945;


            auto dataset = PointProcessDataset::load(events,ar_order,hasTheta0,wp,current_time);

            EXPECT_TRUE(dataset.wn.isApprox(expected_wn));
            EXPECT_TRUE(dataset.xn.isApprox(expected_xn));
            EXPECT_TRUE(dataset.xt.isApprox(expected_xt));
            EXPECT_TRUE(dataset.eta.isApprox(expected_eta));
            EXPECT_TRUE(dataset.wt == expected_wt);

        }

        TEST_F(TestPointProcessDatasetFixture, testPointProcessDataset_hasTheta0True) {

            bool hasTheta0 = true;

            Eigen::VectorXd expected_xt(ar_order + hasTheta0);
            Eigen::MatrixXd expected_xn(n_samples, ar_order + hasTheta0);

            expected_xt << 1.000, 0.930, 0.953, 0.953;
            expected_xn << 1.000, 0.874, 0.891, 0.906,
                    1.000, 0.922, 0.874, 0.891,
                    1.000, 0.946, 0.922, 0.874,
                    1.000, 0.945, 0.946, 0.922,
                    1.000, 0.953, 0.945, 0.946,
                    1.000, 0.953, 0.953, 0.945;


            auto dataset = PointProcessDataset::load(events,ar_order,hasTheta0,wp,current_time);

            EXPECT_TRUE(dataset.wn.isApprox(expected_wn));
            EXPECT_TRUE(dataset.xn.isApprox(expected_xn));
            EXPECT_TRUE(dataset.xt.isApprox(expected_xt));
            EXPECT_TRUE(dataset.eta.isApprox(expected_eta));
            EXPECT_TRUE(dataset.wt == expected_wt);


        }

        TEST(TestToeplitz, testToeplitz1){

            std::vector<double> col {1.0,2.0,3.0};
            std::vector<double> row {1.0,5.0,6.0};
            Eigen::MatrixXd expected(col.size(),row.size());
            expected << 1.0,5.0,6.0,
                    2.0,1.0,5.0,
                    3.0,2.0,1.0;

            Eigen::MatrixXd res(col.size(),row.size());
            PointProcessDataset::toeplitz(col,row,res);
            EXPECT_TRUE(res.isApprox(expected));
        }

        TEST(TestToeplitz, testToeplitz2){
            std::vector<double> col {1.0,2.0,3.0,4.0};
            std::vector<double> row {1.0,6.0,7.0};
            Eigen::MatrixXd expected(col.size(),row.size());
            expected << 1.0,6.0,7.0,
                    2.0,1.0,6.0,
                    3.0,2.0,1.0,
                    4.0,3.0,2.0;
            Eigen::MatrixXd res(col.size(),row.size());
            PointProcessDataset::toeplitz(col,row,res);

            EXPECT_TRUE(res.isApprox(expected));
        }
        TEST(TestToeplitz, testToeplitz3){
            std::vector<double> col {5.0,6.0,7.0};
            std::vector<double> row {5.0,2.0,3.0,4.0};
            Eigen::MatrixXd expected(col.size(),row.size());
            expected << 5.0,2.0,3.0,4.0,
                    6.0,5.0,2.0,3.0,
                    7.0,6.0,5.0,2.0;
            Eigen::MatrixXd res(col.size(),row.size());
            PointProcessDataset::toeplitz(col,row,res);
            EXPECT_TRUE(res.isApprox(expected));
        }
    }
}

#endif //POINTPROCESS_TESTPOINTPROCESSDATASET_H
