//
// Created by Andrea Bonvini on 20/04/21.
//

#ifndef POINTPROCESS_POINTPROCESSDATASET_H
#define POINTPROCESS_POINTPROCESSDATASET_H

#include <utility>
#include <vector>
#include <deque>
#include <numeric>
#include <Eigen/Core>
#include "WeightsProducer.h"

#include <iostream>


class PointProcessDataset{
public:
    unsigned char N_SAMPLES;
    unsigned char AR_ORDER;
    bool hasTheta0;
    Eigen::MatrixXd xn;
    Eigen::VectorXd wn;
    Eigen::VectorXd eta;
    Eigen::VectorXd xt;
    double wt;
    PointProcessDataset(
            unsigned char N_SAMPLES_,
            unsigned char AR_ORDER_,
            bool hasTheta0_,
            Eigen::MatrixXd xn_,
            Eigen::VectorXd wn_,
            Eigen::VectorXd eta_,
            Eigen::VectorXd xt_,
            double wt_
    );

    static PointProcessDataset load(
            std::deque<double> events_times,
            unsigned char AR_ORDER_,
            bool hasTheta0_,
            WeightsProducer& weightsProducer,
            double current_time = 0.0
    );

    static void toeplitz(const std::vector<double>& col, const std::vector<double>& row, Eigen::MatrixXd& toep);
};


#endif //POINTPROCESS_POINTPROCESSDATASET_H
