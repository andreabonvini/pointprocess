//
// Created by Andrea Bonvini on 22/04/21.
//

#ifndef POINTPROCESS_BASEOPTIMIZER_H
#define POINTPROCESS_BASEOPTIMIZER_H


#include <utility>
#include <vector>
#include <numeric>
#include <iostream>
#include <cmath>
#include <memory>

#include "../WeightsProducer.h"
#include "../PointProcessDataset.h"
#include "../PointProcessUtils.h"
#include "../InterEventDistributions.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <deque>
#include <utility>



class BaseOptimizer{
    // TODO: The only public methods should be singleRegression() and train(), for now everything is public
    //  for testing purposes.
public:
    /* BaseOptimizer is the base class which handles the optimization procedure to retrieve
     * the optimal distribution parameters for a given PointProcessDataset.
     * You should refer to either the library scientific documentation, or to the original paper
     * "A point-process model of human heartbeat intervals: new definitions of heart rate and heart rate variability"
     * (Riccardo Barbieri, Eric C. Matten, Abdul Rasheed A. Alabi and Emery N. Brown) for further details.
     */
    PointProcessDistributions distribution;
    explicit BaseOptimizer(PointProcessDistributions distribution);
    virtual ~BaseOptimizer() = default;

    /*
     * This is the most important function of the BaseOptimizer class. It returns a RegressionResult object which
     * contains, among other things, the optimal distribution parameters for a given PointProcessDataset. Currently, the
     * implemented optimization algorithm consists of a combination of both the Raphson-Newton method or gradient
     * descent. TODO: Add Blog-Post to the scientific docs which explains the algorithm in detail.
     */
    std::shared_ptr<pp::RegressionResult> optimizeNewton(
            const PointProcessDataset& dataset,
            bool rightCensoring,
            unsigned long maxIter,
            Eigen::VectorXd& x,
            pp::TmpVars& vars
    );

    /*
     * This function constructs a RegressionResult object starting from the optimal parameters x and some additional
     * variables. It's called at the end ot the optimizeNewton() procedure. It's virtual since the number of parameters
     * can change depending on the distribution (e.g. for the Inverse Gaussian we have to save the scale parameter kappa
     * too)
     */
    virtual std::shared_ptr<pp::RegressionResult> packResult(const Eigen::VectorXd& x, const PointProcessDataset& dataset, bool rightCensoring, unsigned long nIter, double maxGrad);

    /*
     * Compute a trivial estimate for the standard deviation (or kappa parameter for the Inverse Gaussian) based on the
     * sample variance computed on the observed inter-event intervals. It's virtual since the first parameter x[0]
     * doesn't always correspond to the standard deviation (e.g. for the Inverse Gaussian it corresponds to the scale
     * parameter kappa).
     */
    virtual double estimate_x0(const PointProcessDataset& dataset);

    /*
     * The following pure virtual methods are distribution-specific. In order to add a new Optimizer (i.e. a new distribution)
     * it is sufficient to override the following methods (and, if needed, packResult and estimate_x0).
     */
    // First initialization of distribution parameters used in the train() method.
    // TODO: The initialization should be dependent on the current dataset.
    virtual void populateStartingPoint(Eigen::VectorXd& startingPoint) = 0;
    // Compute gradient vector based on the current parameters x and dataset.
    virtual void updateGradient(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradient) = 0;
    // Compute gradient vector for the right-censoring term based on the current parameters x and dataset.
    virtual void updateGradientRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::VectorXd& gradientRc) = 0;
    // Compute hessian matrix based on the current parameters x and dataset.
    virtual void updateHessian(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessian) = 0;
    // Compute the probability density function for the right censoring term based on the current parameters x.
    virtual double computePDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) = 0;
    // Compute the cumulative density function for the right censoring term based on the current parameters x.
    virtual double computeCDF(const Eigen::VectorXd& x, const PointProcessDataset& dataset) = 0;
    // Compute the negative log likelihood function based on the current parameters x.
    virtual double computeLikel(const Eigen::VectorXd& x, const PointProcessDataset& dataset) = 0;

    /*
     * Once the computePDF and computeCDF functions are implemented in a derived class,
     * the computation of the Hazard Rate (Lambda) can be easily computed (see definition in BaseOptimizer.cpp)
     */
    double computeLambda(const Eigen::VectorXd& x, const PointProcessDataset& dataset);

    /*
     * Once the computePDF and computeCDF functions are implemented in a derived class,
     * the computation of the negative log-likelihood for the right-censoring term can be easily computed
     * (see definition in BaseOptimizer.cpp)
     */
    double computeLikelRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset);

    /*
     * The computation of the hessian matrix for the right-censoring term is always approximated through the
     * finite-difference method (see definition in BaseOptimizer.cpp)
     */
    void updateHessianRc(const Eigen::VectorXd& x, const PointProcessDataset& dataset, Eigen::MatrixXd& hessianRc);

    /*
     * Call the optimizeNewton() method taking into account rightCensoring (or not)
     */
    std::shared_ptr<pp::RegressionResult> singleRegression(PointProcessDataset& dataset, bool rightCensoring, unsigned int maxIter);

    /*
     * Call the optimizeNewton() method multiple times based on the options defined in the PipelineSetup.
     */
    std::vector<std::shared_ptr<pp::RegressionResult>> train(pp::PipelineSetup& setup);
};

#endif //POINTPROCESS_BASEOPTIMIZER_H
