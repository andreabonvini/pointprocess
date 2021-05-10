//
// Created by Andrea Bonvini on 22/04/21.
//

#ifndef POINTPROCESS_INVERSEGAUSSIANOPTIMIZER_H
#define POINTPROCESS_INVERSEGAUSSIANOPTIMIZER_H

#include "BaseOptimizer.h"
#include "../InterEventDistributions.h"
#include "../PointProcessDataset.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <deque>
#include <utility>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <fstream>
#include <boost/math/distributions/normal.hpp>

struct IGRegressionResult : RegressionResult {

    double kappa;
    IGRegressionResult(double theta0_,
                       const Eigen::VectorXd& thetaP_,
                       double mu,
                       double sigma,
                       double lambda,
                       double meanInterval,
                       long nIter,
                       double kappa_)
                       : RegressionResult(theta0_, thetaP_, mu, sigma, lambda, meanInterval, nIter) {
        kappa = kappa_;
    }

};


class InverseGaussianOptimizer : public BaseOptimizer{
public:

    explicit InverseGaussianOptimizer(OptimizerSetup setup): BaseOptimizer(setup){
        distribution = PointProcessDistributions::InverseGaussian;
    };

    void populateStartingPoint(Eigen::VectorXd& startingPoint) override{
        // Kappa = x[0]
        // Theta = x.segment(1,x.size() - 1)
        startingPoint[0] = 1500.0; // k0
        startingPoint.segment(1,startingPoint.size() - 1).setConstant(1.0 / double (startingPoint.size()));
    };

    void updateGradient(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradient) override{
        // Kappa = x[0]
        // Theta = x.segment(1,x.size() - 1)
        gradient <<
        0.5 * dataset.eta.dot((- 1.0 / x[0] +(dataset.wn - (dataset.xn * x.segment(1,x.size() - 1))).array().pow(2.0)/((dataset.xn * x.segment(1,x.size() - 1)).array().pow(2.0) * dataset.wn.array())).matrix()),  // derivative w.r.t. kappa
        dataset.xn.transpose() * ( - 1 * x[0] * dataset.eta.array() * ((dataset.wn.array() - (dataset.xn * x.segment(1,x.size() - 1)).array())/(dataset.xn * x.segment(1,x.size() - 1)).array().pow(3.0))).matrix() ;   // derivatives w.r.t. theta

    };

    void updateGradientRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradientRc) override{
        // Kappa = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)

        boost::math::normal norm;

        gradientRc <<
        dataset.eta[dataset.eta.size() - 1] / (1.0 - computeCDF(x,dataset)) * (
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
        dataset.eta[dataset.eta.size() - 1] / (1.0 - computeCDF(x,dataset)) * (
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

    };

    void updateHessian(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessian) override{
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
    };

    double computePDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override {
        // Kappa = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
        return exp(0.5 * log(x[0] / (2.0 * M_PI * pow(dataset.wt, 3.0))) - 0.5 * (x[0] * pow(dataset.wt - x.segment(1,x.size() - 1).dot(dataset.xt), 2.0)) / (pow(x.segment(1,x.size() - 1).dot(dataset.xt), 2.0) * dataset.wt));
    };

    double computeCDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override {
        // Kappa = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
        boost::math::normal norm;
        return cdf(norm, sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) - 1.0)) + exp( 2.0 * x[0] / x.segment(1,x.size() - 1).dot(dataset.xt) + log(cdf(norm, - sqrt(x[0] / dataset.wt) * (dataset.wt / x.segment(1,x.size() - 1).dot(dataset.xt) + 1.0))));
    };

    double computeLikel(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override{
        // Kappa = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // Mus = (dataset.xn * x.segment(1,x.size() - 1))
        if ( x[0] < 0.0 || (dataset.xn * x.segment(1,x.size() - 1)).minCoeff() < 0.0 ) return INFINITY; // Check constraints.
        return - dataset.eta.dot(
                (((x[0] / (2.0 * M_PI * dataset.wn.array().pow(3.0))).sqrt().log() - (( x[0] * (dataset.wn - (dataset.xn * x.segment(1,x.size() - 1))).array().pow(2.0))/(2.0 * (dataset.xn * x.segment(1,x.size() - 1)).array().pow(2.0) * dataset.wn.array()))).matrix()));
    };

    std::shared_ptr<RegressionResult> packResult(const Eigen::VectorXd& x, const PointProcessDataset& dataset, unsigned long nIter) override{

        // Kappa = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)

        // Check constraints
        assert (x[0] > 0.0);
        assert ((dataset.xn * x.segment(1,x.size() - 1)).minCoeff() > 0.0 );

        double meanInterval = dataset.eta.dot(dataset.wn) / dataset.eta.array().sum();
        double mu = dataset.xt.dot(x.segment(1,x.size() - 1));
        double sigma = sqrt(pow(mu,3.0) / x[0]); // Variance = Mu^3 /Kappa

        return std::make_shared<IGRegressionResult>(
                dataset.hasTheta0 ? x[1] : 0.0,
                dataset.hasTheta0 ? x.segment(2,x.size() - 2) : x.segment(1,x.size() - 1),
                mu,
                sigma,
                (dataset.wt > 0.0) ? computeLambda(x,dataset) : 0.0,
                meanInterval,
                nIter,
                x[0]);
    };

};


#endif //POINTPROCESS_INVERSEGAUSSIANOPTIMIZER_H
