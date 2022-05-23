//
// Created by Andrea Bonvini on 05/05/21.
//

#ifndef POINTPROCESS_LOGNORMALOPTIMIZER_H
#define POINTPROCESS_LOGNORMALOPTIMIZER_H

#include "BaseOptimizer.h"
#include "../PointProcessDataset.h"
#include <Eigen/Core>
#include <boost/math/distributions/inverse_gaussian.hpp>


class LogNormalOptimizer : public BaseOptimizer{
public:
    explicit LogNormalOptimizer(): BaseOptimizer(PointProcessDistributions::LogNormal){}

    void populateStartingPoint(Eigen::VectorXd& startingPoint) override{
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
         startingPoint[0] = 0.02; // sigma0
         startingPoint.segment(1,startingPoint.size() - 1).setConstant(1.0 / double (startingPoint.size()));
    }

    void updateGradient(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradient) override{
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // Mus = (dataset.xn * x.segment(1,x.size() - 1))
        gradient <<
        - dataset.eta.dot((- 1.0 / x[0] + (dataset.wn.array().log() - (dataset.xn * x.segment(1,x.size() - 1)).array().log()).pow(2.0) / pow(x[0],3.0) ).matrix()),
        dataset.xn.transpose() * (- dataset.eta.array() * 1.0 / (pow(x[0], 2.0) * (dataset.xn * x.segment(1,x.size() - 1)).array()) * (dataset.wn.array().log() - (dataset.xn * x.segment(1,x.size() - 1)).array().log()) ).matrix();
    }

    void updateGradientRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradientRc) override{
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
        boost::math::normal norm;

        double rcEta = 1.0; // dataset.eta[dataset.eta.size() - 1];
        gradientRc <<
               - rcEta / (1.0 - computeCDF(x,dataset)) * (
                       pdf(norm, (log(dataset.wt) - log(x.segment(1,x.size() - 1).dot(dataset.xt))) / x[0] ) * (log(dataset.wt) - log(x.segment(1,x.size() - 1).dot(dataset.xt))) / pow(x[0],2.0)
               )
               ,
               - rcEta / (1.0 - computeCDF(x,dataset)) * (
                       pdf(norm, (log(dataset.wt) - log(x.segment(1,x.size() - 1).dot(dataset.xt))) / x[0] ) * 1.0 / (x[0] * x.segment(1,x.size() - 1).dot(dataset.xt))
               )
               * dataset.xt;
    }

    void updateHessian(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessian) override{
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // Mus = (dataset.xn * x.segment(1,x.size() - 1))

        // 2nd order derivative w.r.t. Sigma
        hessian(0,0) = -dataset.eta.dot( (1.0 / pow(x[0],2.0) - 3.0 * (dataset.wn.array().log() - (dataset.xn * x.segment(1,x.size() - 1)).array().log()).pow(2.0) / pow(x[0], 4.0) ).matrix() );

        // partial derivative w.r.t Sigma and Theta (or Theta and Sigma)
        hessian.row(0).segment(1,dataset.AR_ORDER + dataset.hasTheta0) <<
        (
                dataset.eta.array()
                *
                ( 2.0 * (dataset.wn.array().log() - (dataset.xn * x.segment(1,x.size() - 1)).array().log()) / (pow(x[0], 3.0) * (dataset.xn * x.segment(1,x.size() - 1)).array()))
        ).matrix().transpose()
        * dataset.xn;

        hessian.col(0).segment(1,dataset.AR_ORDER + dataset.hasTheta0 ) << hessian.row(0).segment(1,dataset.AR_ORDER + dataset.hasTheta0).transpose();

        // 2nd order derivative w.r.t. Theta
        hessian.bottomRightCorner(dataset.AR_ORDER + dataset.hasTheta0,dataset.AR_ORDER + dataset.hasTheta0) <<
        dataset.xn.transpose()
        *
        (
                dataset.eta.array() * (1.0 + dataset.wn.array().log() - (dataset.xn * x.segment(1,x.size() - 1)).array().log()) / (pow(x[0],2.0) * (dataset.xn * x.segment(1,x.size() - 1)).array().pow(2.0) )
        ).matrix().asDiagonal()
        * dataset.xn;
    }

    double computePDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override {
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
        return exp(
                log(1.0 / (dataset.wt * x[0] * sqrt(2.0 * M_PI)))
                -
                pow((log(dataset.wt) - log(x.segment(1,x.size() - 1).dot(dataset.xt))), 2.0) / (2.0 * pow(x[0],2.0))
                )  ;
    }

    double computeCDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override {
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // rcMu = x.segment(1,x.size() - 1).dot(dataset.xt)
        boost::math::normal norm;
        return cdf(norm, (log(dataset.wt) - log(x.segment(1,x.size() - 1).dot(dataset.xt))) /  x[0] );
    }

    double computeLikel(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override{
        // Sigma = x[0]
        // Theta = x.segment(1,x.size() - 1)
        // Mus = (dataset.xn * x.segment(1,x.size() - 1))
        if (
                (x[0] < 0.0)
                ||
                ( (dataset.xn * x.segment(1,x.size() - 1)).array().minCoeff() < 0 )
            ) return INFINITY; // Check constraints.
        else{
            return - dataset.eta.dot( ((1.0 / (dataset.wn.array() * x[0] * sqrt(2.0 * M_PI))).log() - 0.5 * (dataset.wn.array().log() - (dataset.xn * x.segment(1,x.size() - 1)).array().log()).pow(2.0) / pow(x[0],2.0)).matrix() );
        }
    }

    unsigned int getNumberOfAdditionalParams() override {
        // sigma (std of natural logarithm)
        return 1;
    }
};


#endif //POINTPROCESS_LOGNORMALOPTIMIZER_H
