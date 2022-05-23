//
// Created by Andrea Bonvini on 22/04/21.
//

#include "BaseOptimizer.h"


BaseOptimizer::BaseOptimizer(PointProcessDistributions distribution){
    this->distribution = distribution;
}

// LCOV_EXCL_START
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
    bool NOTHING_WORKED;
    bool cdfIsOne = false;

    auto oldold = vars.xold;
    bool rightCensoring_ = rightCensoring && dataset.wt > 1e-5;

    while (iter < maxIter){
        vars.xold = x;
        // -------- compute gradient, hessian and negloglikelihood with the current parameters theta and kappa ---------

        // For efficiency reasons we first update the gradient and the negloglikelihood alone,
        //  and see if there's any need for further optimization.
        updateGradient(x, dataset, vars.gradient);
        negloglikel = computeLikel(x,dataset);
        if (rightCensoring_){
            updateGradientRc(x, dataset, vars.rcGradient);
            vars.gradient += vars.rcGradient;
            negloglikelRc = computeLikelRc(x, dataset);
            negloglikel = negloglikel + negloglikelRc;
        }
        maxGrad = vars.gradient.array().abs().maxCoeff();

        if (std::isinf(negloglikelRc)){
            // FIXME: If we enter this branch it may be cause of numerical problems (CDF becomes 1.0 since dataset.wt
            //  is too high and the current distribution (defined by the parameters x) isn't long-tailed enough.
       /*     std::cout << "\nSomething is wrong, the negative log-likelihood is INFINITY." << std::endl;
            std::cout << "This is probably due to the right censoring term." << std::endl;
            std::cout << "negloglikelRc: " << negloglikelRc << std::endl;
            std::cout << "dataset.wt: " << dataset.wt << std::endl;
            std::cout << "maxGrad: " << maxGrad << std::endl;*/
            negloglikel = 0.0; // FIXME:  Completely arbitrary!
            cdfIsOne = true;
            // FIXME: Are we sure this is what's happening?
            //  In case, how should we behave? We could force the distribution to have a longer tail.
            //  (or just ignore the problem as we're doing currently).
        }


    if (maxGrad <= gradTol || cdfIsOne){
            // There's no need to optimize further, the starting point x is already good enough.
            // TODO: It would be useful to log something in case cdfIsOne is true
            //  (maybe directly serialize it in the RegressionResult)
            break;
        }

        // We can now compute the hessian matrix and proceed with the optimization procedure
        updateHessian(x, dataset, vars.hessian);

        if (rightCensoring_){
            updateHessianRc(x, dataset, vars.rcHessian);
            vars.hessian += vars.rcHessian;
        }

        // ---------- update parameters with Raphson-Newton method or gradient descent ---------------------------
        vars.eigenSolver.compute(vars.hessian);
        NEWTON_WORKED = false;
        GRADIENT_DESCENT_WORKED = false;
        NOTHING_WORKED = false;
        if (vars.eigenSolver.eigenvalues().real().minCoeff() > 0.0){
            // Newton-Raphson with line-search
            for (char i = 0; i < 15; i++){
                vars.alpha.setConstant(1.0 / pow(2.0, (double) i));
                x = vars.xold - (vars.alpha.array() * (vars.hessian.inverse() * vars.gradient).array()).matrix();
                // Compute new likelihood
                tmpNegloglikel = computeLikel(x, dataset);
                if (tmpNegloglikel != INFINITY && rightCensoring_){
                    tmpNegloglikelRc = computeLikelRc(x, dataset);
                    tmpNegloglikel += tmpNegloglikelRc;
                }
                else if (!(rightCensoring_)){
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
        if (!NEWTON_WORKED){
            for (char i = 0; i < 15; i++) {
                vars.alpha.setConstant(0.0005 / pow(2.0, (double) i));
                x = vars.xold - (vars.alpha.array() * vars.gradient.array()).matrix();
                // Compute new likelihood
                tmpNegloglikel = computeLikel(x, dataset);
                if (tmpNegloglikel != INFINITY && rightCensoring_) {
                    tmpNegloglikelRc = computeLikelRc(x, dataset);
                    tmpNegloglikel = tmpNegloglikel + tmpNegloglikelRc;
                }
                else if (!(rightCensoring_)){
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
            /*std::cout << "\nBAD -> negloglikel:" << negloglikel << std::endl;
            std::cout <<   "BAD ->     maxGrad:" << maxGrad << std::endl;*/
            NOTHING_WORKED = true;
            break;
        }
        iter++;
    }

    return packResult(x,dataset,negloglikel, rightCensoring, iter, maxGrad, maxGrad < gradTol, cdfIsOne);
}
// LCOV_EXCL_STOP

std::shared_ptr<pp::RegressionResult> BaseOptimizer::packResult(const Eigen::VectorXd& x, const PointProcessDataset& dataset, double negloglikelihood, bool rightCensoring, unsigned long nIter, double maxGrad, bool converged, bool cdfIsOne) {

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
            negloglikelihood,
            maxGrad,
            converged,
            cdfIsOne);
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


    return optimizeNewton(dataset, rightCensoring, maxIter, x, vars);

}

std::vector<std::shared_ptr<pp::RegressionResult>> BaseOptimizer::train(DatasetBuffer& datasetBuffer, bool rightCensoring, unsigned long maxIter) {

    unsigned int numberOfParams = datasetBuffer.getNumberOfRegressionParameters() + getNumberOfAdditionalParams();
    // Initialize results vector and some useful variables
    std::vector<std::shared_ptr<pp::RegressionResult>> results;
    results.reserve(datasetBuffer.size());

    unsigned int tmpIter;

    auto startingPoint = Eigen::VectorXd(numberOfParams);
    populateStartingPoint(startingPoint);

    auto x = Eigen::VectorXd(numberOfParams);
    x = startingPoint;
    // Note that x will be updated inside updateNewton() and will represent the optimal distribution
    // parameters at each time bin.

    // Declare some useful variables that will be used in the updateNewton() routine.
    auto gradient = Eigen::VectorXd(numberOfParams);
    auto hessian =  Eigen::MatrixXd(numberOfParams,numberOfParams);
    auto rcGradient = Eigen::VectorXd(numberOfParams);
    auto rcHessian = Eigen::MatrixXd(numberOfParams,numberOfParams);
    auto xold = Eigen::VectorXd(numberOfParams);
    auto alpha = Eigen::VectorXd(numberOfParams);
    auto vars = pp::TmpVars(
            gradient,
            hessian,
            rcGradient,
            rcHessian,
            xold,
            alpha
    );

    bool firstRegression = true;
    unsigned long time_step = 1;
    // Main loop
    for (auto [currentTime, eventHappened, resetParameters, dataset] : datasetBuffer){
        pp::utils::logging::printProgress(currentTime, (double) time_step/ (double) datasetBuffer.size());
        time_step++;
        if (resetParameters){

            // The uncensored solution is a good starting point.
            if (dataset.AR_ORDER==0 && dataset.hasTheta0){
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

            std::shared_ptr<pp::RegressionResult> tmpRes = optimizeNewton(dataset, false,( firstRegression ? 50000 : maxIter), x, vars);
            tmpIter = tmpRes -> nIter;
            resetParameters = false;
            firstRegression = false;

        }

        // Compute the actual parameters by applying right-censoring (if specified)
        std::shared_ptr<pp::RegressionResult> result = optimizeNewton(dataset, rightCensoring, maxIter, x, vars);

        // Append metadata to result
        result->time = currentTime; // TODO: Not really good practice, should be in the constructor.
        result->eventHappened = eventHappened;
        if(resetParameters){
            result->nIter += tmpIter;
        }
        results.push_back(result);
    }
    return results;
}
