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

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


namespace pp {

    struct Stats {
        double ksDistance = 0.0;
        double percOut = 0.0;
        double autoCorr = 0.0;

        Stats(double ksDistance, double percOut, double autoCorr);
    };

    struct TmpVars {
        Eigen::VectorXd gradient;
        Eigen::MatrixXd hessian;
        Eigen::VectorXd rcGradient;
        Eigen::MatrixXd rcHessian;
        Eigen::VectorXd xold;
        Eigen::VectorXd alpha;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;

        TmpVars(
                Eigen::VectorXd gradient,
                Eigen::MatrixXd hessian,
                Eigen::VectorXd rcGradient,
                Eigen::MatrixXd rcHessian,
                Eigen::VectorXd xold,
                Eigen::VectorXd alpha
        );
    };

    struct PipelineSetup {
        double delta;
        std::vector<double> events;
        bool hasTheta0;
        unsigned char AR_ORDER;
        unsigned long last_event_index;
        unsigned long bins;
        unsigned long bins_in_window;
        WeightsProducer weightsProducer;

        PipelineSetup(
                double delta,
                std::vector<double> events,
                bool hasTheta0,
                unsigned char AR_ORDER,
                unsigned long last_event_index,
                unsigned long b,
                unsigned long biw,
                WeightsProducer weightsProducer
        );
    };
    // LCOV_EXCL_START
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
        bool eventHappened;
        double maxGrad;
        bool converged;
        bool cdfIsOne;

        RegressionResult(
                double theta0,
                Eigen::VectorXd thetaP,
                double mu,
                double sigma,
                double lambda,
                double meanInterval,
                long nIter,
                double likelihood,
                double maxGrad,
                bool converged,
                bool cdfIsOne
        );

        // I declare a virtual destructor just to have run-time type information (RTTI), which is needed
        // to guarantee polymorphic behaviour.
        virtual ~RegressionResult();
    };
    // LCOV_EXCL_STOP

    // LCOV_EXCL_START
    struct IGRegressionResult : RegressionResult {
        double kappa;

        IGRegressionResult(
                double theta0_,
                const Eigen::VectorXd &thetaP_,
                double mu,
                double sigma,
                double lambda,
                double meanInterval,
                long nIter,
                double likelihood,
                double maxGrad,
                double kappa,
                bool converged,
                bool cdfIsOne
        );
    };
    // LCOV_EXCL_STOP

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
        );
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
        PipelineSetup getPipelineSetup(const std::vector<double> &events, bool hasTheta0_, unsigned char AR_ORDER_,
                         double windowLength, double delta, WeightsProducer weightsProducer);

        // TODO: add documentation
        void computeTaus(std::vector<double> &taus, const std::vector<double> &lambdas, const PipelineSetup &setup);

        namespace serialize {

            void ppResData2csv(Result &ppRes, const std::string &outputResultsName);

            void ppResTaus2csv(Result &ppRes, const std::string &outputTausName);
        }

        namespace logging {
            void printProgress(double currentTime, double percentage);
        }
    }
}

#endif //POINTPROCESS_POINTPROCESSUTILS_H