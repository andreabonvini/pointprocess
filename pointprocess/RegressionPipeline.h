//
// Created by Andrea Bonvini on 21/04/21.
//

#ifndef POINTPROCESS_REGRESSIONPIPELINE_H
#define POINTPROCESS_REGRESSIONPIPELINE_H

#include "InterEventDistributions.h"
#include "optimizers/BaseOptimizer.h"
#include "optimizers/InverseGaussianOptimizer.h"
#include "optimizers/GaussianOptimizer.h"
#include "optimizers/LogNormalOptimizer.h"
#include "WeightsProducer.h"
#include <utility>
#include <vector>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <iostream>


static void computeTaus(std::vector<double>& taus ,const std::vector<std::shared_ptr<RegressionResult>>& results, const OptimizerSetup& setup){

    double currentTime;
    bool eventHappened;
    unsigned long last_event_index = setup.last_event_index;
    double dt;
    double intL = 0.0;
    double pr = 0.0;
    bool wait = true;

    unsigned long offset = setup.bins_in_window;
    for (unsigned long bin_index = offset; bin_index <= setup.bins + 1; bin_index++) {
        currentTime = (double) (bin_index - 1) * setup.delta;
        eventHappened = setup.events[last_event_index + 1] <= currentTime;
        if (eventHappened){
            last_event_index++;
            dt = setup.events[last_event_index] - (currentTime - setup.delta);
            intL = intL + dt * pr;
            if (!wait) {
                taus.push_back(intL);
            }
            // we set wait to false when we observe the FIRST event after the starting timeWindow, this way we'll
            // be able to compute the integral of lambda (i.e. tau) only for fully observed inter-event intervals.
            /* e.g.
             *              1  2    3   4    5    6  7
             * events:    --|--|----|---|----|----|--|
             *   bins:    ----------------------------
             * timeWindow:  \_____________/XX       <- WRONG
             *                  \_____________/     <- RIGHT
             *  We have to wait the first observable event (5) in order to compute a reliable estimate
             *  of the first tau (integral of lambda)
             */
            wait = false;
            intL = 0.0;
        }
        else{
            intL = intL + setup.delta * results[bin_index - offset]->lambda;
        }
        if (bin_index <= results.size()){
            pr =  results[bin_index - offset]->lambda;
        }
    }
}



class RegressionPipeline{
public:
    PointProcessDistributions interEventDistribution;
    unsigned char AR_ORDER;
    bool hasTheta0;
    RegressionPipeline(
            PointProcessDistributions interEventDistribution,
            unsigned char AR_ORDER,
            bool hasTheta0
            ){
        /******************************************************************
         * Parameters:
         *     interEventDistribution_: One of the enumeration values of PointProcessDistributions, either
         *         1) InverseGaussian
         *         2) LogNormal
         *         3) Gaussian
         *     AR_ORDER_: AR order to use for the estimation of the first moment of the given distribution.
         *     hasTheta0_: if the AR model takes account for a mean/theta0 parameter.
         *
         *****************************************************************/
        this->interEventDistribution = interEventDistribution;
        this->AR_ORDER = AR_ORDER;
        this->hasTheta0 = hasTheta0;

    }

    [[nodiscard]] PointProcessResult fullRegression(
            const std::vector<double>& events_times,
            double windowLength = 60.0,
            double delta = 0.005,
            bool rightCensoring = true,
            unsigned int maxIter = 1000,
            WeightsProducer weightsProducer = WeightsProducer()
            ) const{
        /**************************************************************************************************************
         * This function implements part of the pipeline suggested by Riccardo Barbieri, Eric C. Matten,
         * Abdul Rasheed A. Alabi and Emery N. Brown in the paper:
         * "A point-process model of human heartbeat intervals:
         *  new definitions of heart rate and heart rate variability".
         *
         * Check the various Optimizer objects to see how the optimization process is carried out for the different
         * distributions.
         *
         * Parameters:
         *     events_times: event times expressed in seconds.
         *     windowLength: time window used for the local likelihood maximization.
         *     delta: how much the local likelihood time interval is shifted to compute the next parameter update,
         *            be careful: time_resolution must be little enough s.t. at most ONE event can happen in each
         *            time bin. Moreover the smaller it is the better since we use it to approximate the integral
         *            of the lambda function.
         *     rightCensoring: if the regression should take into account right-censoring or not, if true we should have
         *                     more accurate estimates for the first and second moment of the selected distribution.
         *     maxIter: maximum number of iterations allowed for each optimization procedure.
         ************************************************************************************************************/

        // We want the first event to be at time 0 (events = events_times - events_times[0])
        std::vector<double> events = events_times;
        double t0 = events_times[0];
        std::transform(events_times.begin(), events_times.end(), events.begin(),[&](auto& value){ return value - t0;});
        auto pipelineSetup = getOptimizerSetup(t0, events, rightCensoring, hasTheta0, AR_ORDER, windowLength, delta,
                                               maxIter, weightsProducer);

        // TODO: factorize...
        switch (this->interEventDistribution) {
            case PointProcessDistributions::InverseGaussian: {
                auto optimizer = InverseGaussianOptimizer(pipelineSetup);
                auto results = optimizer.train();
                std::vector<double> taus;
                computeTaus(taus, results, pipelineSetup);
                // TODO: compute ksDistance & Co.
                double percOut = 0.0;
                double ksDistance = 0.0;
                double autoCorr = 0.0;
                return PointProcessResult(results,taus,percOut,ksDistance, pipelineSetup.t0, autoCorr);
            }
            case PointProcessDistributions::LogNormal: {
                auto optimizer = GaussianOptimizer(pipelineSetup);
                auto results = optimizer.train();
                std::vector<double> taus;
                computeTaus(taus, results, pipelineSetup);
                // TODO: compute ksDistance & Co.
                double percOut = 0.0;
                double ksDistance = 0.0;
                double autoCorr = 0.0;
                return PointProcessResult(results,taus,percOut,ksDistance, pipelineSetup.t0, autoCorr);
            }
            case PointProcessDistributions::Gaussian: {
                auto optimizer = GaussianOptimizer(pipelineSetup);
                auto results = optimizer.train();
                std::vector<double> taus;
                computeTaus(taus, results, pipelineSetup);
                // TODO: compute ksDistance & Co.
                double percOut = 0.0;
                double ksDistance = 0.0;
                double autoCorr = 0.0;
                return PointProcessResult(results,taus,percOut,ksDistance, pipelineSetup.t0, autoCorr);
            }
            default:
                throw std::logic_error("Please, insert a valid InterEvent distribution.");

        }


    }

    [[nodiscard]] std::shared_ptr<IGRegressionResult> singleRegression(const std::vector<double>& events_times, unsigned long maxIter = 1000, WeightsProducer weightsProducer = WeightsProducer()) const{

        switch (this->interEventDistribution) {
            case PointProcessDistributions::InverseGaussian: {
//                std::deque<double> et;
//                for (const auto& e: events_times)
//                    et.push_back(e);
//
//                auto dataset = PointProcessDataset::load(et, AR_ORDER, hasTheta0, weightsProducer, 0.0);
//                double kappaStart = 1500.0;
//                auto thetaStart = Eigen::VectorXd(AR_ORDER + hasTheta0);
//                thetaStart.setConstant(1.0 / double (AR_ORDER + hasTheta0));
//                auto gradient = Eigen::VectorXd(AR_ORDER + hasTheta0 + 1);
//                auto hessian =  Eigen::MatrixXd(AR_ORDER + hasTheta0 + 1,AR_ORDER + hasTheta0 + 1);
//                auto rcGradient = Eigen::VectorXd(AR_ORDER + hasTheta0 + 1);
//                auto rcHessian = Eigen::MatrixXd(AR_ORDER + hasTheta0 + 1,AR_ORDER + hasTheta0 + 1);
//                auto vars = igTmpVars(
//                        gradient,
//                        hessian,
//                        rcGradient,
//                        rcHessian,
//                        &rcData,
//                        &rcTmpData
//                );
//                return optimizeNewton(dataset, false, maxIter, kappaStart, thetaStart, vars);
            }
            case LogNormal:
                throw std::logic_error("LogNormal singleRegression hasn't been implemented yet!");
            case Gaussian:
                throw std::logic_error("Gaussian singleRegression hasn't been implemented yet!");
        }
    }
};


#endif //POINTPROCESS_REGRESSIONPIPELINE_H

