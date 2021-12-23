//
// Created by Andrea Bonvini on 22/04/21.
//

#include "BaseOptimizer.h"



BaseOptimizer::BaseOptimizer(PointProcessDistributions distribution){
    this->distribution = distribution;
}


std::shared_ptr<pp::RegressionResult> BaseOptimizer::optimizeNewton(
        const PointProcessDataset& dataset,
        bool rightCensoring,
        unsigned long maxIter,
        Eigen::VectorXd& x,
        pp::TmpVars& vars
){
    double oldNegloglikel = INFINITY;
    double negloglikel;
    double negloglikelRc;
    double tmpNegloglikel;
    double tmpNegloglikelRc;

    double gradTol = 1e-5;
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
            if (!cdfIsOne){ // TODO: put a break? Log this information somewhere...
                // FIXME: this shouldn't happen
                // Detected cdf == 1.0 during right censoring. Maybe you lost an event while annotating the data...
                // Maybe you should use a distribution with a longer tail (e.g. LogNormal)
                // std::cout << "Detected cdf == 1.0 during right censoring. Maybe you lost an event while annotating the data... " << std::endl;
            }
            cdfIsOne = true;
        }
        else{
            cdfIsOne = false;
        }

        // ---------- update parameters with Raphson-Newton method or gradient descent ---------------------------
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
                    // time bin the sigma/kappa parameters are greater than 0.0
                    if (x.segment(1,x.size() - 1).dot(dataset.xt) < 0.0 || x[0] < 0.0) tmpNegloglikel = INFINITY;
                }
                // If the new point is better than the old one we can stop and go on with the next time bin.
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

    return packResult(x,dataset,rightCensoring, iter, maxGrad);
}

std::shared_ptr<pp::RegressionResult> BaseOptimizer::packResult(const Eigen::VectorXd& x, const PointProcessDataset& dataset, bool rightCensoring, unsigned long nIter, double maxGrad) {

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

    return std::make_shared<pp::RegressionResult>(
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




// Once the computePDF and computeCDF functions are implemented in a derived class, the computation of the Hazard Rate (Lambda)
// and the censored negative log likelihood can be carried out here in the BaseOptimizer class. Moreover since we always approximate the
// hessian for the right censoring term, we can generalize its computation here.
double BaseOptimizer::computeLambda(const Eigen::VectorXd& x, const PointProcessDataset& dataset){
    return computePDF(x,dataset) / (1.0 - computeCDF(x,dataset));
}

double BaseOptimizer::computeLikelRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset){
    // Sigma = x[0]
    // Theta = x.segment(1,x.size() - 1)
    // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
    double rcEta = 1.0; // dataset.eta[dataset.eta.size() - 1];
    if (
            x.segment(1,x.size() - 1).dot(dataset.xt) < 0.0
            ||
            x.segment(1,x.size() - 1).dot(dataset.xt) > dataset.wn.array().maxCoeff()
            ) return INFINITY;
    return - rcEta * log(1.0 - computeCDF(x,dataset));
};

void BaseOptimizer::updateHessianRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessianRc) {
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

double BaseOptimizer::estimate_x0(const PointProcessDataset& dataset) {
    assert(dataset.wn.size() > 2); // the number of target events (dataset.wn) is < 2, can't estimate x0!
    double mu_hat = dataset.eta.dot(dataset.wn) / dataset.eta.array().sum();
    double var = 1.0 / ((double)dataset.wn.size() - 1.0) * (dataset.wn.array() - mu_hat).pow(2.0).sum();
    // Variance = 1 / (n - 1) * sum([w - mu_hat for w in wn])
    // sigma = sqrt(Variance)
    return sqrt(var);
}

std::shared_ptr<pp::RegressionResult> BaseOptimizer::singleRegression(PointProcessDataset& dataset, bool rightCensoring = false, unsigned int maxIter = 1000){

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
    auto alpha = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    auto vars =  pp::TmpVars(
            gradient,
            hessian,
            rcGradient,
            rcHessian,
            xold,
            alpha
    );


    return optimizeNewton(dataset, rightCensoring && dataset.wt > 0.0, maxIter, x, vars);

}

std::vector<std::shared_ptr<pp::RegressionResult>> BaseOptimizer::train(pp::PipelineSetup& setup) {

    unsigned long last_event_index = setup.last_event_index;
    auto observed_events = std::deque<double>(setup.events.begin(), setup.events.begin() + (long) last_event_index + 1);

    /* observed_events here is the subset of events observed during the first window, this std::deque will keep track
     * of the events used for local regression at each time bin, discarding old events and adding new ones.
     * It works as a buffer for our regression pipeline.
     */

    // Initialize results vector and some useful variables
    std::vector<std::shared_ptr<pp::RegressionResult>> results;
    results.reserve(setup.bins - setup.bins_in_window + 1);
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
    auto vars = pp::TmpVars(
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
            if (setup.AR_ORDER == 0 && setup.hasTheta0){
                // If we are directly estimating the mu parameter, we can trivially set as a starting point the
                // mean inter-event interval
                x[1] = dataset.eta.dot(dataset.wn) / dataset.eta.array().sum();
                // TODO: do something like that for the AutoRegressive estimate too.
            }
            startingPoint = x;
            // We can estimate sigma (or kappa) given our target intervals.
            double tmp_x0 = dataset.wn.size() <= 2 ? x[0]: estimate_x0(dataset);
            startingPoint[0] = tmp_x0;

            if (computeLikel(startingPoint,dataset) < computeLikel(x,dataset)){
                // If our estimate for the first parameter is better w.r.t. the old value, we set our estimate as x[0].
                // Remember that this is a negative log likelihood, the lower the better.
                x[0] = tmp_x0;
            }

            std::shared_ptr<pp::RegressionResult> tmpRes = optimizeNewton(dataset, false,( (bin_index == setup.bins_in_window) ? 50000 : setup.maxIter), x, vars);
            tmpIter = tmpRes -> nIter;
            resetParameters = false;
        }

        // Compute the actual parameters by applying right-censoring (if specified)
        std::shared_ptr<pp::RegressionResult> result = optimizeNewton(dataset, setup.rightCensoring && dataset.wt > 0.0, setup.maxIter, x, vars);

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
