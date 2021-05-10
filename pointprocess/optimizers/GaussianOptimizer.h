//
// Created by Andrea Bonvini on 05/05/21.
//

#ifndef POINTPROCESS_GAUSSIANOPTIMIZER_H
#define POINTPROCESS_GAUSSIANOPTIMIZER_H

#include "BaseOptimizer.h"
#include "../PointProcessDataset.h"
#include <Eigen/Core>
#include <boost/math/distributions/inverse_gaussian.hpp>


class GaussianOptimizer : public BaseOptimizer{
public:
    explicit GaussianOptimizer(OptimizerSetup setup): BaseOptimizer(setup){
        distribution = PointProcessDistributions::Gaussian;
    };

    void populateStartingPoint(Eigen::VectorXd& startingPoint) override{
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        startingPoint[0] = 0.0003; // sigma0
        startingPoint.segment(1,startingPoint.size() - 1).setConstant(1.0 / double (startingPoint.size()));
    };

    void updateGradient(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradient) override{

        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // Mus = dataset.xn * x.segment(1,x.size() - 1)
        gradient <<
        - dataset.eta.dot((- 1.0 / x[0] + 0.5 * (dataset.wn - dataset.xn * x.segment(1,x.size() - 1)).array().pow(2.0)/pow(x[0],2.0)).matrix()),
        - dataset.xn.transpose() * (dataset.eta.array() * ( 1.0 / x[0] * (dataset.wn.array() - (dataset.xn * x.segment(1,x.size() - 1)).array()))).matrix();

    };

    void updateGradientRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradientRc) override{
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
        // rcEta = dataset.eta[dataset.eta.size() - 1]

        boost::math::normal norm;

        gradientRc <<
        dataset.eta[dataset.eta.size() - 1] / (1.0 - computeCDF(x,dataset)) * (
                - (dataset.wt - x.segment(1,x.size() - 1).dot(dataset.xt)) / pow(x[0],2.0) * pdf(norm,(dataset.wt - x.segment(1,x.size() - 1).dot(dataset.xt)) / x[0])
        )
        ,
        dataset.eta[dataset.eta.size() - 1] / (1.0 - computeCDF(x,dataset)) * (
                - 1.0 / x[0] * pdf(norm,(dataset.wt - x.segment(1,x.size() - 1).dot(dataset.xt)) / x[0])
        )
        * dataset.xt;
    };

    void updateHessian(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessian) override{
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // Mus = dataset.xn * x.segment(1,x.size() - 1)

        // 2nd order derivative w.r.t. Sigma
        hessian(0,0) = -dataset.eta.dot( (1.0 / pow(x[0],2.0) - (dataset.wn - dataset.xn * x.segment(1,x.size() - 1)).array().pow(2.0) / pow(x[0],3.0)).matrix() );

        // partial derivative w.r.t Sigma and Theta (or Theta and Sigma)
        hessian.row(0).segment(1,dataset.AR_ORDER + dataset.hasTheta0) <<  (dataset.eta.array() * ( 1.0 / pow(x[0],2.0) * (dataset.wn.array() - (dataset.xn * x.segment(1,x.size() - 1)).array()))).matrix().transpose() * dataset.xn;
        hessian.col(0).segment(1,dataset.AR_ORDER + dataset.hasTheta0 ) << hessian.row(0).segment(1,dataset.AR_ORDER + dataset.hasTheta0).transpose();

        // 2nd order derivative w.r.t. Theta
        hessian.bottomRightCorner(dataset.AR_ORDER + dataset.hasTheta0,dataset.AR_ORDER + dataset.hasTheta0) <<
        dataset.xn.transpose()
        *
        (dataset.eta.array() / x[0]).matrix().asDiagonal()
        * dataset.xn;
    };

    double computePDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override {
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
//        std::cout.precision(30);
//        std::cout << exp(log( 1.0 / (x[0] * sqrt(2.0 * M_PI))) - 0.5 * (x[0] * pow((dataset.wt - x.segment(1,x.size() - 1).dot(dataset.xt)) / x[0], 2.0)) ) << std::endl ;
        return exp(log( 1.0 / (x[0] * sqrt(2.0 * M_PI))) - 0.5 * (x[0] * pow((dataset.wt - x.segment(1,x.size() - 1).dot(dataset.xt)) / x[0], 2.0)) );
    };

    double computeCDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override {
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
//        std::cout << "CDF" << std::endl;
//        std::cout << "rcMu: " << x.segment(1,x.size() - 1).dot(dataset.xt) << std::endl;
//        std::cout << "Sigma: " << x[0] << std::endl;
//        boost::math::normal gauss(x.segment(1,x.size() - 1).dot(dataset.xt),x[0]);

        double res = (1.0 + std::erf( (dataset.wt - x.segment(1,x.size() - 1).dot(dataset.xt)) / (x[0] * sqrt(2.0)) ))/ 2.0;

        if (res == 1.0){
//            std::cout << "AAAAAAA" << std::endl;
//            std::cout << "Sigma: " << x[0] << std::endl;
//            std::cout << "rcMu: " << x.segment(1,x.size() - 1).dot(dataset.xt) << std::endl;
//            std::cout << "wt: " << dataset.wt << std::endl;
//            std::cout << "res: " << res << std::endl;
        }
        return std::min(res, 1.0 - 1e-15);
    };


    double computeLikel(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override{
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // Mus = dataset.xn * x.segment(1,x.size() - 1)
        if ( x[0] < 0.0 || (dataset.xn * x.segment(1,x.size() - 1)).minCoeff() < 0.0 ) return INFINITY; // Check constraints.
        return - dataset.eta.dot(((
                log(1.0 / (x[0]*sqrt(2.0 * M_PI)))
                -
                0.5 * (dataset.wn.array() - (dataset.xn * x.segment(1,x.size() - 1)).array()).pow(2.0) / x[0]
                ).matrix()));
    };


};


#endif //POINTPROCESS_GAUSSIANOPTIMIZER_H
