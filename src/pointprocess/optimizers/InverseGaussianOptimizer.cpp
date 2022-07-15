//
// Created by Andrea Bonvini on 22/04/21.
//

#include "InverseGaussianOptimizer.h"


InverseGaussianOptimizer::InverseGaussianOptimizer(): BaseOptimizer(pointprocess::Distributions::InverseGaussian){};

void InverseGaussianOptimizer::populateStartingPoint(Eigen::VectorXd& startingPoint){
    // Kappa = x[0]
    // Theta = x.segment(1,x.size() - 1)
    startingPoint[0] = 1000.0; // k0
    startingPoint.segment(1,startingPoint.size() - 1).setConstant(1.0 / double (startingPoint.size()));
}

void InverseGaussianOptimizer::updateGradient(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradient){
    // Kappa = x[0]
    // Theta = x.segment(1,x.size() - 1)
    gradient <<
    0.5 * dataset.eta.dot((- 1.0 / x[0] +(dataset.wn - (dataset.xn * x.segment(1,x.size() - 1))).array().pow(2.0)/((dataset.xn * x.segment(1,x.size() - 1)).array().pow(2.0) * dataset.wn.array())).matrix()),  // derivative w.r.t. kappa
    dataset.xn.transpose() * ( - 1 * x[0] * dataset.eta.array() * ((dataset.wn.array() - (dataset.xn * x.segment(1,x.size() - 1)).array())/(dataset.xn * x.segment(1,x.size() - 1)).array().pow(3.0))).matrix() ;   // derivatives w.r.t. theta
}

void InverseGaussianOptimizer::updateGradientRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradientRc){
    // Kappa = x[0]
    // Theta = x.segment(1,x.size() - 1)
    // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)

    boost::math::normal norm;

    double rcEta = 1.0; // dataset.eta[dataset.eta.size() - 1];

    gradientRc <<
    rcEta / (1.0 - computeCDF(x,dataset)) * (
    (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) - 1.0)/(2.0 * dataset.wt * sqrt(x[0]/dataset.wt))
    *
    pdf(norm, sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) - 1.0))
    +
    2.0 / x.segment(1,x.size() - 1).dot(dataset.xt)
    *
    exp(2.0 * x[0] / x.segment(1,x.size() - 1).dot(dataset.xt) + log(cdf(norm,- sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) + 1.0))))
    -
    (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) + 1.0)/(2.0 * dataset.wt * sqrt(x[0]/dataset.wt))
    *
    exp(2.0 * x[0] / x.segment(1,x.size() - 1).dot(dataset.xt) + log(pdf(norm,- sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) + 1.0))))
    )
    ,
    rcEta / (1.0 - computeCDF(x,dataset)) * (
    (- sqrt(x[0] / dataset.wt) * dataset.wt / pow(x.segment(1,x.size() - 1).dot(dataset.xt),2.0))
    *
    pdf(norm,sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) - 1.0))
    -
    2.0 * x[0] / pow(x.segment(1,x.size() - 1).dot(dataset.xt),2.0)
    *
    exp(2.0 * x[0] / x.segment(1,x.size() - 1).dot(dataset.xt) + log(cdf(norm,- sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) + 1.0))))
    +
    sqrt(x[0] / dataset.wt) * dataset.wt / pow(x.segment(1,x.size() - 1).dot(dataset.xt),2.0)
    *
    exp(2.0 * x[0] / x.segment(1,x.size() - 1).dot(dataset.xt) + log(pdf(norm,- sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) + 1.0))))
    )
    * dataset.xt;

}

void InverseGaussianOptimizer::updateHessian(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessian){
    // Kappa = x[0]
    // Theta = x.segment(1,x.size() - 1)

    // 2nd order derivative w.r.t. k
    hessian(0,0) = 0.5 * dataset.eta.sum() * (1 / pow(x[0],2.0));

    // partial derivative w.r.t k and theta (or theta and k)
    hessian.row(0).segment(1,dataset.AR_ORDER + dataset.hasTheta0) <<  - (dataset.eta.array() * (dataset.wn - (dataset.xn * x.segment(1,x.size() - 1))).array() / (dataset.xn * x.segment(1,x.size() - 1)).array().pow(3.0)).matrix().transpose() * dataset.xn;
    hessian.col(0).segment(1,dataset.AR_ORDER + dataset.hasTheta0 ) << hessian.row(0).segment(1,dataset.AR_ORDER + dataset.hasTheta0).transpose();

    // 2nd order derivative w.r.t. theta
    hessian.bottomRightCorner(dataset.AR_ORDER + dataset.hasTheta0,dataset.AR_ORDER + dataset.hasTheta0)
    <<
    dataset.xn.transpose() *
    (x[0] * dataset.eta.array() * ((3.0 * dataset.wn - 2.0 * (dataset.xn * x.segment(1,x.size() - 1))).array() / (dataset.xn * x.segment(1,x.size() - 1)).array().pow(4.0))).matrix().asDiagonal()
    * dataset.xn;
}

double InverseGaussianOptimizer::computePDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset){
    // Kappa = x[0]
    // Theta = x.segment(1,x.size() - 1)
    // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
    return exp(0.5 * log(x[0] / (2.0 * M_PI * pow(dataset.wt, 3.0))) - 0.5 * (x[0] * pow(dataset.wt - x.segment(1,x.size() - 1).dot(dataset.xt), 2.0)) / (pow(x.segment(1,x.size() - 1).dot(dataset.xt), 2.0) * dataset.wt));
}

double InverseGaussianOptimizer::computeCDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset){
    // Kappa = x[0]
    // Theta = x.segment(1,x.size() - 1)
    // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
    boost::math::normal norm;
    double arg1 = sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) - 1.0);
    double arg2 = - sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) + 1.0);
    return cdf(norm, arg1) + exp( 2.0 * x[0] / x.segment(1,x.size() - 1).dot(dataset.xt) + log(cdf(norm, arg2)));
}

double InverseGaussianOptimizer::computeLikel(const Eigen::VectorXd& x, const PointProcessDataset& dataset){
    // Kappa = x[0]
    // Theta = x.segment(1,x.size() - 1)
    // Mus = (dataset.xn * x.segment(1,x.size() - 1))
    if ((x[0] < 0.0)||( (dataset.xn * x.segment(1,x.size() - 1)).array().minCoeff() < 0 ))
        return INFINITY; // Check constraints.
    else {
        return -dataset.eta.dot(
                (((x[0] / (2.0 * M_PI * dataset.wn.array().pow(3.0))).sqrt().log() -
                ((x[0] * (dataset.wn - (dataset.xn * x.segment(1, x.size() - 1))).array().pow(2.0)) /
                (2.0 * (dataset.xn * x.segment(1, x.size() - 1)).array().pow(2.0) *
                dataset.wn.array()))).matrix()));
    }
}

std::shared_ptr<pointprocess::RegressionResult> InverseGaussianOptimizer::packResult(const Eigen::VectorXd& x, const PointProcessDataset& dataset, double negloglikelihood, bool rightCensoring, unsigned long nIter, double maxGrad, bool converged, bool cdfIsOne, double time, bool eventHappened){

    // Kappa = x[0]
    // Theta = x.segment(1,x.size() - 1)
    // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)

    // Check constraints
    assert (x[0] > 0.0);
    assert ((dataset.xn * x.segment(1,x.size() - 1)).minCoeff() > 0.0 );
    assert (x.segment(1,x.size() - 1).dot(dataset.xt) > 0.0);

    double meanInterval = dataset.eta.dot(dataset.wn) / dataset.eta.array().sum();
    double mu = dataset.xt.dot(x.segment(1,x.size() - 1));
    double sigma = sqrt(pow(mu,3.0) / x[0]); // Variance = Mu^3 /Kappa
    double kappa = x[0];

    return std::make_shared<pointprocess::IGRegressionResult>(
            dataset.hasTheta0 ? x[1] : 0.0,
            dataset.hasTheta0 ? x.segment(2,x.size() - 2) : x.segment(1,x.size() - 1),
            mu,
            sigma,
            (dataset.wt > 0.0) ? computeLambda(x,dataset) : 0.0,
            meanInterval,
            nIter,
            negloglikelihood,
            maxGrad,
            x[0],
            converged,
            cdfIsOne,
            eventHappened,
            time
            );
}

double InverseGaussianOptimizer::estimate_x0(const PointProcessDataset& dataset) {
    // We can use the sample variance in order to have a coarse estimate for kappa (x0 for the InverseGaussian distribution)
    assert(dataset.wn.size() > 2); // the number of target events (dataset.wn) is < 2, can't estimate x0!
    double mu_hat = dataset.eta.dot(dataset.wn) / dataset.eta.array().sum();
    double var = 1.0 / ((double)dataset.wn.size() - 1.0) * (dataset.wn.array() - mu_hat).pow(2.0).sum();
    // Variance = 1 / (n - 1) * sum([w - mu_hat for w in wn])
    // Kappa = mu_hat^3 / Variance
    return pow(mu_hat,3.0) / var;
    }

unsigned int InverseGaussianOptimizer::getNumberOfAdditionalParams() {
    // kappa (scale parameter)
    return 1;
}