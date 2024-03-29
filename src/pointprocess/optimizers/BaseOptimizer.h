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
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <deque>
#include <utility>

#include "../WeightsProducer.h"
#include "../PointProcessDataset.h"
#include "../InterEventDistributions.h"
#include "../DatasetBuffer.h"
#include "../../external/indicators.h"



#ifndef M_PI  // For Visual Studio
    #define M_PI 3.14159265358979323846
#endif



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
    pointprocess::Distributions distribution;
    explicit BaseOptimizer(pointprocess::Distributions distribution);
    virtual ~BaseOptimizer();

    /*
     * This is the most important function of the BaseOptimizer class. It returns a RegressionResult object which
     * contains, among other things, the optimal distribution parameters for a given PointProcessDataset. Currently, the
     * implemented optimization algorithm consists of a combination of both the Raphson-Newton method or gradient
     * descent. TODO: Add Blog-Post to the scientific docs which explains the algorithm in detail.
     */
    std::shared_ptr<pointprocess::RegressionResult> optimizeNewton(
            const PointProcessDataset& dataset,
            bool rightCensoring,
            unsigned long maxIter,
            Eigen::VectorXd& x,
            pointprocess::TmpVars& vars,
            double time = 0.0,
            bool eventHappened = false
    );

    /*
     * This function constructs a RegressionResult object starting from the optimal parameters x and some additional
     * variables. It's called at the end ot the optimizeNewton() procedure. It's virtual since the number of parameters
     * can change depending on the distribution (e.g. for the Inverse Gaussian we have to save the scale parameter kappa
     * too)
     */
    virtual std::shared_ptr<pointprocess::RegressionResult> packResult(const Eigen::VectorXd& x, const PointProcessDataset& dataset, double negloglikelihood,  bool rightCensoring, unsigned long nIter, double maxGrad, bool converged, bool cdfIsOne, double time, bool eventHappened);

    /*
     * Compute a trivial estimate for the standard deviation (or kappa parameter for the Inverse Gaussian) based on the
     * sample variance computed on the observed inter-event intervals. It's virtual since the first parameter x[0]
     * doesn't always correspond to the standard deviation (e.g. for the Inverse Gaussian it corresponds to the scale
     * parameter kappa).
     * TODO: In the future the number of additional parameters mey be higher than 1 (e.g. see the Beta distribution),
     *  so it may be reasonable to change the following function in estimateAdditionalParameters() and return a
     *  std::vector<double>
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
    // Every distribution, have a variable number of parameters other than the mean (e.g. the sigma for the Normal distribution
    //  or the scale parameter for the Inverse Gaussian), this function simply represents the number of parameters other than
    virtual unsigned int getNumberOfAdditionalParams() = 0;

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
    std::shared_ptr<pointprocess::RegressionResult> singleRegression(PointProcessDataset& dataset, bool rightCensoring, unsigned long maxIter);

    /*
     * Call the optimizeNewton() method multiple times based on the options defined in the PipelineSetup.
     */
    std::vector<std::shared_ptr<pointprocess::RegressionResult>> train(DatasetBuffer& datasetBuffer, bool rightCensoring, unsigned long maxIter);
};

#endif //POINTPROCESS_BASEOPTIMIZER_H
