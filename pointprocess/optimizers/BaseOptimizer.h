//
// Created by Andrea Bonvini on 22/04/21.
//

#ifndef POINTPROCESS_BASEOPTIMIZER_H
#define POINTPROCESS_BASEOPTIMIZER_H


#include <utility>
#include <vector>
#include <numeric>
#include <iostream>

#include "../WeightsProducer.h"
#include "../PointProcessDataset.h"
#include "../InterEventDistributions.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <deque>
#include <utility>


struct TmpVars{
    Eigen::VectorXd gradient;
    Eigen::MatrixXd hessian;
    Eigen::VectorXd rcGradient;
    Eigen::MatrixXd rcHessian;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;
    Eigen::VectorXd xold;
    Eigen::VectorXd alpha;
    TmpVars(
            Eigen::VectorXd& gradient,
            Eigen::MatrixXd& hessian,
            Eigen::VectorXd& rcGradient,
            Eigen::MatrixXd& rcHessian,
            Eigen::VectorXd& xold,
            Eigen::VectorXd& alpha){
        this->gradient = gradient;
        this->hessian = hessian;
        this->rcGradient = rcGradient;
        this->rcHessian = rcHessian;
        this->xold = xold;
        this->alpha = alpha;
    }
};

struct PipelineSetup{
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
            double delta_,
            std::vector<double> e,
            bool rc,
            bool hasTheta0_,
            unsigned char AR_ORDER_,
            unsigned long lei,
            unsigned long b,
            unsigned long biw,
            unsigned long maxIter_,
            WeightsProducer weightsProducer_
            ){
        delta = delta_;
        events = std::move(e);
        rightCensoring = rc;
        hasTheta0 = hasTheta0_;
        AR_ORDER = AR_ORDER_;
        last_event_index = lei;
        bins = b;
        bins_in_window = biw;
        maxIter = maxIter_;
        weightsProducer = weightsProducer_;
    }

};


struct RegressionResult{
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
            double theta0_,
            Eigen::VectorXd thetaP_,
            double mu_,
            double sigma_,
            double lambda_,
            double meanInterval_,
            long nIter_,
            double likelihood_,
            double maxGrad_
            ){
        theta0 = theta0_;
        thetaP = std::move(thetaP_);
        mu = mu_;
        sigma = sigma_;
        lambda = lambda_;
        meanInterval = meanInterval_;
        nIter = nIter_;
        likelihood = likelihood_;
        maxGrad = maxGrad_;
    }
    virtual ~RegressionResult() = default;
};


class BaseOptimizer{
public:
    PointProcessDistributions distribution;

    explicit BaseOptimizer(PointProcessDistributions distribution){
        this->distribution = distribution;
    }

    virtual std::shared_ptr<RegressionResult> optimizeNewton(
            const PointProcessDataset& dataset,
            bool rightCensoring,
            unsigned long maxIter,
            Eigen::VectorXd& x,
            TmpVars& vars
    ){
        double oldNegloglikel = INFINITY;
        double negloglikel;
        double negloglikelRc;
        double tmpNegloglikel;
        double tmpNegloglikelRc;

        double gradTol = 1e-4;
        double maxGrad = INFINITY;
        unsigned long iter = 0;
        unsigned long newtonIter = 0;
        unsigned long gradientDescentIter = 0;
        bool NEWTON_WORKED;
        bool GRADIENT_DESCENT_WORKED;
        bool cdfIsOne = false;

        auto oldold = vars.xold;

        while (maxGrad > gradTol && iter < maxIter){
            vars.xold = x;
            // -------- compute gradient, hessian and negloglikelihood with the current parameters theta and kappa ---------
            updateGradient(x, dataset, vars.gradient);
            updateHessian(x, dataset, vars.hessian);
            negloglikel = computeLikel(x,dataset);

            if (rightCensoring && dataset.wt > 0.0){
                updateGradientRc(x, dataset, vars.rcGradient);
                updateHessianRc(x, dataset, vars.rcHessian);
                vars.gradient += vars.rcGradient;
                vars.hessian += vars.rcHessian;
                negloglikelRc = computeLikelRc(x, dataset);
                negloglikel = negloglikel + negloglikelRc;
            }

            if (isinf(negloglikelRc)){
//                if (!cdfIsOne){ TODO: UNCOMMENT. // put a break?
//                    std::cout << "Detected cdf == 1.0 during right censoring. Maybe you lost an event while annotating the data... " << std::endl;
//                }
                cdfIsOne = true;
            }
            else{
                cdfIsOne = false;
            }

            // ---------- update theta and kappa with Raphson-Newton method or gradient descent ---------------------------
            vars.eigenSolver.compute(vars.hessian);
            NEWTON_WORKED = false;
            GRADIENT_DESCENT_WORKED = false;
            if (vars.eigenSolver.eigenvalues().real().minCoeff() > 0.0 && !cdfIsOne ){
                // Newton-Raphson with line-search
                for (char i = 0; i < 10; i++){
                    vars.alpha.setConstant(1.0 / pow(2.0, (double) i));
                    x = vars.xold - (vars.alpha.array() * (vars.hessian.inverse() * vars.gradient).array()).matrix();
                    // Compute new likelihood
                    tmpNegloglikel = computeLikel(x, dataset);
                    if (tmpNegloglikel != INFINITY && rightCensoring){
                        tmpNegloglikelRc = computeLikelRc(x, dataset);
                        tmpNegloglikel += tmpNegloglikelRc;
                    }
                    else if (!rightCensoring){
                        // Even if we don't apply rightCensoring we want to assure that the estimate for the current
                        // time bin is greater than 0.0
                        if (x.segment(1,x.size() - 1).dot(dataset.xt) < 0.0) tmpNegloglikel = INFINITY;
                    }
                    // If the new point is better than the old one we go on, otherwise we decrease the learning rate.
                    if (tmpNegloglikel != INFINITY && tmpNegloglikel < negloglikel){
                        negloglikel = tmpNegloglikel;
                        NEWTON_WORKED = true;
                        break;
                    }
                    else{
                        continue;
                    }
                }
            }
            if (!NEWTON_WORKED && !cdfIsOne){
                for (char i = 0; i < 10; i++) {
                    vars.alpha.setConstant(0.0005 / pow(2.0, (double) i));
                    x = vars.xold - (vars.alpha.array() * vars.gradient.array()).matrix();
                    // Compute new likelihood
                    tmpNegloglikel = computeLikel(x, dataset);
                    if (tmpNegloglikel != INFINITY && rightCensoring) {
                        tmpNegloglikelRc = computeLikelRc(x, dataset);
                        tmpNegloglikel = tmpNegloglikel + tmpNegloglikelRc;
                    }
                    else if (!rightCensoring){
                        // Even if we don't apply rightCensoring we want to assure that the estimate for the current
                        // time bin is greater than 0.0
                        if (x.segment(1,x.size() - 1).dot(dataset.xt) < 0.0) tmpNegloglikel = INFINITY;
                    }
                    // If the new point is better than the old one we go on.
                    if (tmpNegloglikel != INFINITY && tmpNegloglikel < negloglikel) {
                        negloglikel = tmpNegloglikel;
                        GRADIENT_DESCENT_WORKED = true;
                        break;
                    }
                    else {
                        continue;
                    }
                }
            }

            if (NEWTON_WORKED){
                newtonIter++;
            }

            else if (GRADIENT_DESCENT_WORKED){
                gradientDescentIter++;
            }
            else{
                x = vars.xold;
                break;
            }
            iter++;
            oldNegloglikel = negloglikel;
            maxGrad = vars.gradient.array().abs().maxCoeff();
        }

        if (x.segment(1,x.size() - 1).dot(dataset.xt) < 0.0){
            std::cout << "dataset.wt: " << dataset.wt << std::endl;
            std::cout << "rcMu: " << x.segment(1,x.size() - 1).dot(dataset.xt) << std::endl;
            std::cout << "oldold x: \n" << oldold << std::endl;
            std::cout << "current x: \n" << x << std::endl;
            std::cout << "maxGrad: " << maxGrad << std::endl;
            std::cout << "gradient: \n" << vars.gradient << std::endl;
            std::cout << "rcGradient: \n" << vars.rcGradient << std::endl;
            std::cout << "NEWTON WORKED: " << NEWTON_WORKED << std::endl;
            std::cout << "GRADIENT DESCENT WORKED: " << GRADIENT_DESCENT_WORKED << std::endl;
            std::cout << "cdfIsOne: " << cdfIsOne << std::endl;
            std::cout << "iter: " << iter << std::endl;
            std::cout << "negloglikel: " << negloglikel << std::endl;
        }


        return packResult(x,dataset,rightCensoring, iter, maxGrad);
    }

    virtual std::shared_ptr<RegressionResult> packResult(const Eigen::VectorXd& x, const PointProcessDataset& dataset, bool rightCensoring, unsigned long nIter, double maxGrad) {

        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)

        // Check constraints
        assert (x[0] > 0.0);
        assert ((dataset.xn * x.segment(1,x.size() - 1)).minCoeff() > 0.0 );
        assert (x.segment(1,x.size() - 1).dot(dataset.xt) > 0.0);

        double meanInterval = dataset.eta.dot(dataset.wn) / dataset.eta.array().sum();
        double mu = dataset.xt.dot(x.segment(1,x.size() - 1));
        double sigma = x[0];


        return std::make_shared<RegressionResult>(
                dataset.hasTheta0 ? x[1] : 0.0,
                dataset.hasTheta0 ? x.segment(2,x.size() - 2) : x.segment(1,x.size() - 1),
                mu,
                sigma,
                (dataset.wt > 0.0) ? computeLambda(x,dataset) : 0.0,
                meanInterval,
                nIter,
                computeLikel(x, dataset) + (rightCensoring? computeLikelRc(x, dataset) : 0.0),
                maxGrad);
    };

    // The following functions are distribution-specific.
    virtual void populateStartingPoint(Eigen::VectorXd& startingPoint) = 0;
    virtual void updateGradient(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradient) = 0;
    virtual void updateGradientRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradientRc) = 0;
    virtual void updateHessian(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessian) = 0;
    virtual double computePDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) = 0;
    virtual double computeCDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) = 0;
    virtual double computeLikel(const Eigen::VectorXd& x, const PointProcessDataset& dataset) = 0;

    // Once the computePDF and computeCDF functions are implemented in a derived class, the computation of the Hazard Rate (Lambda)
    // and the censored negative log likelihood can be carried out here in the BaseOptimizer class. Moreover since we always approximate the
    // hessian for the right censoring term, we can generalize its computation here.
    double computeLambda(const Eigen::VectorXd& x, const PointProcessDataset& dataset){
        return computePDF(x,dataset) / (1.0 - computeCDF(x,dataset));
    }

    double computeLikelRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset){
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
        if (x.segment(1,x.size() - 1).dot(dataset.xt) < 0.0) return INFINITY;
        return - dataset.eta[dataset.eta.size() - 1] * log(1.0 - computeCDF(x,dataset));
    };

    void updateHessianRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessianRc) {
        auto gradRight = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        auto gradLeft = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        auto xRight = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        auto xLeft = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        double machinePrecision = std::cbrt(1e-15);
        double step;
        xRight = x;
        xLeft = x;
        for (long i = 0; i < x.size(); i++){
            step = machinePrecision * x[i];
            // Move right and left for parameter i
            xRight[i] += step;
            xLeft[i] -= step;
            // Compute tmp gradient and update hessianRc through finite difference approximation.
            updateGradientRc(xRight, dataset,gradRight);
            updateGradientRc(xLeft, dataset,gradLeft);
            for (long j = 0; j <x.size() ;j++){
                hessianRc(i,j) = (gradRight[j] - gradLeft[j]) / (2.0 * step);
            }
            // Reset...
            xRight[i] = x[i];
            xLeft[i] = x[i];
        }
    };

    std::shared_ptr<RegressionResult> singleRegression(PointProcessDataset& dataset, bool rightCensoring = false, unsigned int maxIter = 1000){

        auto startingPoint = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        populateStartingPoint(startingPoint);
        auto x = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        x = startingPoint;

        // Declare some useful variables that will be used in the updateNewton() routine.
        auto gradient = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        auto hessian =  Eigen::MatrixXd(dataset.AR_ORDER + dataset.hasTheta0 + 1, dataset.AR_ORDER + dataset.hasTheta0 + 1);
        auto rcGradient = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        auto rcHessian = Eigen::MatrixXd(dataset.AR_ORDER + dataset.hasTheta0 + 1,dataset.AR_ORDER + dataset.hasTheta0 + 1);
        auto xold = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        auto alpha =Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
        auto vars =  TmpVars(
                gradient,
                hessian,
                rcGradient,
                rcHessian,
                xold,
                alpha
        );


        return optimizeNewton(dataset, rightCensoring && dataset.wt > 0.0, maxIter, x, vars);

    }

    std::vector<std::shared_ptr<RegressionResult>> train(PipelineSetup& setup) {

        unsigned long last_event_index = setup.last_event_index;
        auto observed_events = std::deque<double>(setup.events.begin(), setup.events.begin() + (long) last_event_index + 1);

        /* observed_events here is the subset of events observed during the first window, this std::deque will keep track
         * of the events used for local regression at each time bin, discarding old events and adding new ones.
         * It works as a buffer for our regression pipeline.
         */

        // Initialize results vector and some useful variables
        std::vector<std::shared_ptr<RegressionResult>> results;
        results.reserve(setup.bins - setup.bins_in_window);
        double currentTime;
        double windowLength = setup.delta * (double) (setup.bins_in_window);
        bool resetParameters = true;
        bool eventHappened;
        unsigned int tmpIter;
        assert (fabs(std::remainder(windowLength,setup.delta)) < 1e-10); // no better way to check if windowLength is a multiple of delta...

        auto startingPoint = Eigen::VectorXd(setup.AR_ORDER + setup.hasTheta0 + 1);
        populateStartingPoint(startingPoint);
        auto x = Eigen::VectorXd(setup.AR_ORDER + setup.hasTheta0 + 1);
        x = startingPoint;
        // Note that x will be updated inside updateNewton() and will represent the optimal distribution
        // parameters at each time bin.

        // Declare some useful variables that will be used in the updateNewton() routine.
        auto gradient = Eigen::VectorXd(setup.AR_ORDER + setup.hasTheta0 + 1);
        auto hessian =  Eigen::MatrixXd(setup.AR_ORDER + setup.hasTheta0 + 1,setup.AR_ORDER + setup.hasTheta0 + 1);
        auto rcGradient = Eigen::VectorXd(setup.AR_ORDER + setup.hasTheta0 + 1);
        auto rcHessian = Eigen::MatrixXd(setup.AR_ORDER + setup.hasTheta0 + 1,setup.AR_ORDER + setup.hasTheta0 + 1);
        auto xold = Eigen::VectorXd(setup.AR_ORDER + setup.hasTheta0 + 1);
        auto alpha = Eigen::VectorXd(setup.AR_ORDER + setup.hasTheta0 + 1);
        TmpVars vars = TmpVars(
                gradient,
                hessian,
                rcGradient,
                rcHessian,
                xold,
                alpha
                );

        // Main loop
        for (unsigned long bin_index = setup.bins_in_window; bin_index <= setup.bins; bin_index ++){
            // TODO: Add some kind of progress bar.
            // TODO: Parallelize?.
            if (bin_index % 10000 == 0){ // TODO: Remove.
                std::cout << bin_index << " / " << setup.bins << "\n";
            }
            currentTime = (double) bin_index * setup.delta;
            /* If the first element of observed_events happened before the
             * time window between (current_time - window_length) and (current_time)
             * we can discard it since it will not be part of the current optimization process.
             */
            if (!observed_events.empty() && observed_events[0] < currentTime - windowLength){
                observed_events.pop_front();
                // Force re-evaluation of starting point for theta and kappa.
                resetParameters = true;
            }
            // We check whether an event happened in ((bin_index - 1) * delta, bin_index * delta]
            eventHappened = setup.events[last_event_index + 1] <= currentTime;
            if (eventHappened){
                last_event_index++;
                observed_events.push_back(setup.events[last_event_index]);
                resetParameters = true;
            }
            // We create a PointProcessDataset for the current time bin
            auto dataset = PointProcessDataset::load(observed_events,setup.AR_ORDER,setup.hasTheta0,setup.weightsProducer, currentTime);
            if (resetParameters){
                // The uncensored solution is a good starting point.
                std::shared_ptr<RegressionResult> tmpRes = optimizeNewton(dataset, false,( (bin_index == setup.bins_in_window) ? 50000 : setup.maxIter), x, vars);
                tmpIter = tmpRes -> nIter;
                resetParameters = false;
            }

            // Compute the actual parameters by applying right-censoring (if specified)
            std::shared_ptr<RegressionResult> result = optimizeNewton(dataset, setup.rightCensoring && dataset.wt > 0.0, setup.maxIter, x, vars);

            // Append metadata to result
            result->time = currentTime;
            result->eventHappened = eventHappened;
            if(eventHappened){
                result->nIter += tmpIter;
            }
            results.push_back(result);
        }

        return results;
    }
};

#endif //POINTPROCESS_BASEOPTIMIZER_H
