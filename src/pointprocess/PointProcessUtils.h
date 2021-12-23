//
// Created by Andrea Bonvini on 23/12/21.
//

#ifndef POINTPROCESS_POINTPROCESSUTILS_H
#define POINTPROCESS_POINTPROCESSUTILS_H

#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "WeightsProducer.h"
#include "InterEventDistributions.h"


namespace pp {

    struct Stats {
        double ksDistance = 0.0;
        double percOut = 0.0;
        double autoCorr = 0.0;

        Stats(double ksDistance, double percOut, double autoCorr) : ksDistance(ksDistance), percOut(percOut),
                                                                    autoCorr(autoCorr) {}
    };

    struct TmpVars {
        Eigen::VectorXd gradient;
        Eigen::MatrixXd hessian;
        Eigen::VectorXd rcGradient;
        Eigen::MatrixXd rcHessian;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;
        Eigen::VectorXd xold;
        Eigen::VectorXd alpha;

        TmpVars(
                Eigen::VectorXd &gradient,
                Eigen::MatrixXd &hessian,
                Eigen::VectorXd &rcGradient,
                Eigen::MatrixXd &rcHessian,
                Eigen::VectorXd &xold,
                Eigen::VectorXd &alpha) {
            this->gradient = gradient;
            this->hessian = hessian;
            this->rcGradient = rcGradient;
            this->rcHessian = rcHessian;
            this->xold = xold;
            this->alpha = alpha;
        }
    };

    struct PipelineSetup {
        double delta;
        std::vector<double> events;
        bool rightCensoring;
        bool hasTheta0;
        unsigned char AR_ORDER;
        unsigned long last_event_index;
        unsigned long bins;
        unsigned long bins_in_window;
        unsigned int maxIter;
        WeightsProducer weightsProducer;

        PipelineSetup(
                double delta,
                std::vector<double> events,
                bool rc,
                bool hasTheta0,
                unsigned char AR_ORDER,
                unsigned long last_event_index,
                unsigned long b,
                unsigned long biw,
                unsigned long maxIter,
                WeightsProducer weightsProducer
        ) :
                delta(delta),
                events(std::move(events)),
                rightCensoring(rc),
                hasTheta0(hasTheta0),
                AR_ORDER(AR_ORDER),
                last_event_index(last_event_index),
                bins(b),
                bins_in_window(biw),
                maxIter(maxIter),
                weightsProducer(weightsProducer) {}
    };

    struct RegressionResult {
        double theta0;
        Eigen::VectorXd thetaP;
        double mu;
        double sigma;
        double lambda;
        double meanInterval;
        double time = 0.0;
        unsigned long nIter;
        double likelihood;
        bool eventHappened = false;
        double maxGrad;

        RegressionResult(
                double theta0,
                Eigen::VectorXd thetaP,
                double mu,
                double sigma,
                double lambda,
                double meanInterval,
                long nIter,
                double likelihood,
                double maxGrad
        ) :
                theta0(theta0),
                thetaP(std::move(thetaP)),
                mu(mu),
                sigma(sigma),
                lambda(lambda),
                meanInterval(meanInterval),
                nIter(nIter),
                likelihood(likelihood),
                maxGrad(maxGrad) {}

        // I declare a virtual destructor just to have run-time type information (RTTI), which is needed
        // to guarantee polymorphic behaviour.
        virtual ~RegressionResult() = default;
    };

    struct IGRegressionResult : RegressionResult {
        double kappa;

        IGRegressionResult(double theta0_,
                           const Eigen::VectorXd &thetaP_,
                           double mu,
                           double sigma,
                           double lambda,
                           double meanInterval,
                           long nIter,
                           double likelihood,
                           double maxGrad,
                           double kappa)
                :
                kappa(kappa),
                RegressionResult(theta0_, thetaP_, mu, sigma, lambda, meanInterval, nIter, likelihood, maxGrad) {}
    };

    struct Result { // TODO: add Documentation
        std::vector<std::shared_ptr<RegressionResult>> results;
        std::vector<double> taus;
        PointProcessDistributions distribution;
        unsigned char AR_ORDER;
        bool hasTheta0;
        double windowLength;
        double delta;
        double t0;
        Stats stats;

        Result(
                std::vector<std::shared_ptr<RegressionResult>> results,
                std::vector<double> taus,
                PointProcessDistributions distribution,
                unsigned char AR_ORDER,
                bool hasTheta0,
                double windowLength,
                double delta,
                double t0,
                Stats stats
        ) :
                results(std::move(results)),
                taus(std::move(taus)),
                distribution(distribution),
                AR_ORDER(AR_ORDER),
                hasTheta0(hasTheta0),
                windowLength(windowLength),
                delta(delta),
                t0(t0),
                stats(stats) {}
    };

    namespace utils {

        /*
         * This function computes a series of useful metrics to assess goodness of fit of a point process model
         * starting from the computed taus (check scientific docs for details):
         *  ksDistance: TODO: add documentation
         *  percOut: TODO: add documentation
         *  autoCorr: TODO: Implement.
         */
        Stats computeStats(std::vector<double> &taus);

        /*
         * This function returns a PipelineSetup object containing a series of useful parameters for a fullRegression(),
         * such as:
         * last_event_index:
         *     index of the last event within the first time window
         *     e.g. if events = [0.0, 1.3, 2.1, 3.2, 3.9, 4.5] and window_length = 3.5 then last_event_index = 3
         *     (since events[3] = 3.2 and events[4] =3.9)
         * bins:
         *     total number of bins we can discretize our events in (given our time_resolution)
         * bins_in_window:
         *     number of bins in a single time window.
         */
        PipelineSetup getPipelineSetup(const std::vector<double> &events, bool rc, bool hasTheta0_, unsigned char AR_ORDER_,
                         double windowLength, double delta, unsigned long maxIter, WeightsProducer weightsProducer);

        // TODO: add documentation
        void computeTaus(std::vector<double> &taus, const std::vector<double> &lambdas, const PipelineSetup &setup);

        namespace serialize {

            void ppResData2csv(Result &ppRes, const std::string &outputResultsName);

            void ppResTaus2csv(Result &ppRes, const std::string &outputTausName);
        }
    }
}

#endif //POINTPROCESS_POINTPROCESSUTILS_H