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
    GaussianOptimizer();

    void populateStartingPoint(Eigen::VectorXd& startingPoint) override;

    void updateGradient(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradient) override;

    void updateGradientRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradientRc) override;

    void updateHessian(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessian) override;

    double computePDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override;

    double computeCDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override;

    double computeLikel(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override;

    unsigned int getNumberOfAdditionalParams() override;
};


#endif //POINTPROCESS_GAUSSIANOPTIMIZER_H
