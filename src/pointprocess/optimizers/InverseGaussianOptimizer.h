//
// Created by Andrea Bonvini on 22/04/21.
//

#ifndef POINTPROCESS_INVERSEGAUSSIANOPTIMIZER_H
#define POINTPROCESS_INVERSEGAUSSIANOPTIMIZER_H

#include "BaseOptimizer.h"
#include "../PointProcessUtils.h"
#include "../InterEventDistributions.h"
#include "../PointProcessDataset.h"
#include "../PointProcessUtils.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <deque>
#include <utility>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <fstream>
#include <boost/math/distributions/normal.hpp>


class InverseGaussianOptimizer : public BaseOptimizer {
public:

    InverseGaussianOptimizer();

    void populateStartingPoint(Eigen::VectorXd& startingPoint) override;

    void updateGradient(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradient) override;

    void updateGradientRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradientRc) override;

    void updateHessian(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessian) override;

    double computePDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override;

    double computeCDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override;

    double computeLikel(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override;

    std::shared_ptr<pointprocess::RegressionResult> packResult(const Eigen::VectorXd& x, const PointProcessDataset& dataset, double negloglikelihood, bool rightCensoring, unsigned long nIter, double maxGrad, bool converged, bool cdfIsOne, double time, bool eventHappened) override;

    double estimate_x0(const PointProcessDataset& dataset) override;

    unsigned int getNumberOfAdditionalParams() override ;

};


#endif //POINTPROCESS_INVERSEGAUSSIANOPTIMIZER_H
