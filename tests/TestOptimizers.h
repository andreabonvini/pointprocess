//
// Created by Andrea Bonvini on 07/05/21.
//

#ifndef POINTPROCESS_TESTOPTIMIZERS_H
#define POINTPROCESS_TESTOPTIMIZERS_H

#include "data/testData.h"
#include "../src/pointprocess/InterEventDistributions.h"
#include "../src/pointprocess/optimizers/GaussianOptimizer.h"

#include <unordered_map>
#include <gtest/gtest.h>


namespace pointprocess::tests {
    namespace {
        // This is a Value-Parameterized fixture, check https://github.com/google/googletest/blob/main/docs/advanced.md
        //  for more details.
        class TestOptimizersFixture : public testing::TestWithParam<PointProcessDistributions> {
        public:
            PointProcessDataset dataset = getTestDataset();
            int a = 3;

            std::vector<double> events = {};
            WeightsProducer wp = WeightsProducer(0.0);
            PipelineSetup set = PipelineSetup(0.005,events,true,true,8,10,100,10,10000,wp);
            std::unordered_map<PointProcessDistributions, std::shared_ptr<BaseOptimizer>> dist2opt = {
                    {PointProcessDistributions::Gaussian, std::make_shared<GaussianOptimizer>(GaussianOptimizer())},
                    {PointProcessDistributions::InverseGaussian, std::make_shared<InverseGaussianOptimizer>(InverseGaussianOptimizer())},
                    {PointProcessDistributions::LogNormal, std::make_shared<LogNormalOptimizer>(LogNormalOptimizer())}
            };
        };


        TEST_P(TestOptimizersFixture, TestGradientRc){

            PointProcessDistributions dist = GetParam();

            Eigen::VectorXd gradientRc(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::VectorXd approxGradientRc(dataset.AR_ORDER + dataset.hasTheta0 + 1);

            Eigen::VectorXd xStart(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::VectorXd xRight(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::VectorXd xLeft(dataset.AR_ORDER + dataset.hasTheta0 + 1);

            auto dummyOpt = dist2opt[dist];
            dummyOpt->populateStartingPoint(xStart);
            if(dist == PointProcessDistributions::Gaussian){
                xStart[0] = 0.1;
            }

            // Compute exact gradient
            dummyOpt->updateGradientRc(xStart,dataset,gradientRc);

            // Approximate gradient
            double machinePrecision = std::cbrt(1e-15);
            double step;
            double negloglikelRcRight;
            double negloglikelRcLeft;

            xRight = xStart;
            xLeft = xStart;
            for(long i = 0; i < xStart.size(); i++){
                step = machinePrecision * xStart[i];
                // Move right and left for parameter i
                xRight[i] += step;
                xLeft[i] -= step;
                // Compute tmp negloglikel and update gradient through finite difference approximation.
                negloglikelRcRight = dummyOpt->computeLikelRc(xRight, dataset);
                negloglikelRcLeft = dummyOpt->computeLikelRc(xLeft, dataset);
                approxGradientRc[i] = (negloglikelRcRight - negloglikelRcLeft) / (2.0 * step);
                // Reset...
                xRight[i] = xStart[i];
                xLeft[i] = xStart[i];
            }
            // Compare
            EXPECT_TRUE(approxGradientRc.isApprox(gradientRc,1e-5));


        }

        TEST_P(TestOptimizersFixture, TestGradient){

            PointProcessDistributions dist = GetParam();

            Eigen::VectorXd gradient(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::VectorXd approxGradient(dataset.AR_ORDER + dataset.hasTheta0 + 1);

            Eigen::VectorXd xStart(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::VectorXd xRight(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::VectorXd xLeft(dataset.AR_ORDER + dataset.hasTheta0 + 1);

            auto dummyOpt = dist2opt[dist];
            dummyOpt->populateStartingPoint(xStart);

            // Compute exact gradient
            dummyOpt->updateGradient(xStart,dataset,gradient);

            // Approximate gradient
            double machinePrecision = std::cbrt(1e-15);
            double step;
            double negloglikelRcRight;
            double negloglikelRcLeft;

            xRight = xStart;
            xLeft = xStart;
            for(long i = 0; i < xStart.size(); i++){
                step = machinePrecision * xStart[i];
                // Move right and left for parameter i
                xRight[i] += step;
                xLeft[i] -= step;
                // Compute tmp negloglikel and update gradient through finite difference approximation.
                negloglikelRcRight = dummyOpt->computeLikel(xRight, dataset);
                negloglikelRcLeft = dummyOpt->computeLikel(xLeft, dataset);
                approxGradient[i] = (negloglikelRcRight - negloglikelRcLeft) / (2.0 * step);
                // Reset...
                xRight[i] = xStart[i];
                xLeft[i] = xStart[i];
            }
            // Compare
            EXPECT_TRUE(approxGradient.isApprox(gradient,1e-5));
        }

        TEST_P(TestOptimizersFixture, TestHessians){



            Eigen::MatrixXd hessian(dataset.AR_ORDER + dataset.hasTheta0 + 1, dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::MatrixXd approxHessian(dataset.AR_ORDER + dataset.hasTheta0 + 1, dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::VectorXd xStart(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            auto gradRight = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            auto gradLeft = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::VectorXd xRight(dataset.AR_ORDER + dataset.hasTheta0 + 1);
            Eigen::VectorXd xLeft(dataset.AR_ORDER + dataset.hasTheta0 + 1);

            std::vector<double> events = {};
            auto wp = WeightsProducer(0.0);
            auto set = PipelineSetup(0.005,events,true,true,8,10,100,10,10000,wp);

            PointProcessDistributions dist = GetParam();
            auto dummyOpt = dist2opt[dist];
            dummyOpt->populateStartingPoint(xStart);

            // Compute exact gradient
            dummyOpt->updateHessian(xStart,dataset,hessian);

            // Approximate gradient
            double machinePrecision = std::cbrt(1e-15);
            double step;
            xRight = xStart;
            xLeft = xStart;
            for (long i = 0; i < xStart.size(); i++){
                step = machinePrecision * xStart[i];
                // Move right and left for parameter i
                xRight[i] += step;
                xLeft[i] -= step;
                // Compute tmp gradient and update hessianRc through finite difference approximation.
                dummyOpt->updateGradient(xRight, dataset,gradRight);
                dummyOpt->updateGradient(xLeft, dataset,gradLeft);
                for (long j = 0; j < xStart.size() ; j++){
                    approxHessian(i,j) = (gradRight[j] - gradLeft[j]) / (2.0 * step);
                }
                // Reset...
                xRight[i] = xStart[i];
                xLeft[i] = xStart[i];
            }
            // Compare
            EXPECT_TRUE(approxHessian.isApprox(hessian,1e-5));
        }

        std::string getParamName(const ::testing::TestParamInfo<TestOptimizersFixture::ParamType>& info){
            // We need this function just to display the distribution information when checking the tests output.
            // https://google.github.io/googletest/reference/testing.html Check INSTANTIATE_TEST_SUITE_P section
            //  for more details.
            std::unordered_map<PointProcessDistributions, std::string> map = {
                {PointProcessDistributions::Gaussian, std::string("Gaussian")},
                {PointProcessDistributions::InverseGaussian, std::string("InverseGaussian")},
                {PointProcessDistributions::LogNormal, std::string("LogNormal")}
            };
            return map[info.param];
        }

        INSTANTIATE_TEST_SUITE_P(
                TestGradientsAndHessians,
                TestOptimizersFixture,
                testing::Values(
                        PointProcessDistributions::Gaussian,
                        PointProcessDistributions::InverseGaussian,
                        PointProcessDistributions::LogNormal
                ),
                getParamName
        );

    }
}


//bool testTrain(){
//
//
//    std::vector<double> testEvents = {0.0,0.997,2.0,3.0,4.0,4.997,6.0};
//    double windowLength = 5.0;
//    double delta = 0.005;
//    unsigned char AR_ORDER = 2;
//    bool hasTheta0 = true;
//    auto setup = getPipelineSetup(testEvents, false, hasTheta0, AR_ORDER, windowLength, delta, 1, WeightsProducer(1.00));
//    /*
//     * TODO: In this test I'm just testing the currentTime, it should be possible to add debugging
//     *  information to the RegressionResults e.g. by adding the current observed_events and dataset used in each
//     *  precise time bin.
//     * Since we will use a windowLength of 5.0 seconds,
//     * In the first window we expect:
//     *    currentTime     : 5.0
//     *    observed_events : 0.0, 0.997, 2.0, 3.0, 4.0, 4.997
//     *    dataset.wt      : 5.0 - 4.997  = 0.003
//     * In the last window we expect:
//     *    currentTime     : 6.0
//     *    observed_events : 2.0, 3.0, 4.0, 4.997, 6.0
//     *    dataset.wt      : 6.0 - 6.0  = 0.0
//     */
//
//    auto opt = InverseGaussianOptimizer();
//    auto results = opt.train(setup);
//    return (results[0]->time == 5.0 && results[results.size() - 1]->time == 6.0);
//}


#endif //POINTPROCESS_TESTOPTIMIZERS_H
