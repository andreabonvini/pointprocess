//
// Created by Andrea Bonvini on 21/04/21.
//

#ifndef POINTPROCESS_REGRESSIONPIPELINE_H
#define POINTPROCESS_REGRESSIONPIPELINE_H


#include "PointProcessUtils.h"
#include "OptimizersFactory.h"


#include <vector>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <iostream>


class RegressionPipeline{
private:
  // One of the enumeration values of PointProcessDistributions, either,  InverseGaussian, LogNormal, Gaussian
  PointProcessDistributions distribution;
  
  // AR order to use for the estimation of the first moment of the given distribution.
  unsigned char AR_ORDER;
  
  // if the AR model takes account for a mean/theta0 parameter.
  bool hasTheta0;
  
public:
  RegressionPipeline(
    PointProcessDistributions distribution,
    unsigned char AR_ORDER,
    bool hasTheta0
    );
  
    /**************************************************************************************************************
         * This function implements the pipeline suggested by Riccardo Barbieri, Eric C. Matten,
         * Abdul Rasheed A. Alabi and Emery N. Brown in the paper:
         * "A point-process model of human heartbeat intervals:
         *  new definitions of heart rate and heart rate variability".
         *
         * Check the various Optimizer objects to see how the optimization process is carried out for the different
         * distributions.
         *
         * Parameters:
         *     events_times: event times expressed in seconds.
         *     windowLength: time window used for the local likelihood maximization.
         *     delta: how much the local likelihood time interval is shifted to compute the next parameter update,
         *            be careful: time_resolution must be little enough s.t. at most ONE event can happen in each
         *            time bin. Moreover the smaller it is the better since we use it to approximate the integral
         *            of the lambda function.
         *     rightCensoring: if the regression should take into account right-censoring or not, if true we should have
         *                     more accurate estimates for the first and second moment of the selected distribution.
         *     maxIter: maximum number of iterations allowed for each optimization procedure.
    ************************************************************************************************************/
    [[nodiscard]] pp::Result fullRegression(
            const std::vector<double>& events_times,
            double windowLength,
            double delta,
            bool rightCensoring,
            unsigned int maxIter,
            WeightsProducer weightsProducer
      ) const;
};


#endif //POINTPROCESS_REGRESSIONPIPELINE_H

