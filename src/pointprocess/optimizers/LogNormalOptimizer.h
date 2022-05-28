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

    LogNormalOptimizer();

    void populateStartingPoint(Eigen::VectorXd& startingPoint) override;

    void updateGradient(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradient) override;

    void updateGradientRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradientRc) override;

    void updateHessian(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessian) override;

    double computePDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override;

    double computeCDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override;

    double computeLikel(const Eigen::VectorXd& x, const PointProcessDataset& dataset) override;

    unsigned int getNumberOfAdditionalParams() override;
};


#endif //POINTPROCESS_LOGNORMALOPTIMIZER_H
