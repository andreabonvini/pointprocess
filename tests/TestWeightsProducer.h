//
// Created by Andrea Bonvini on 20/04/21.
//

#ifndef POINTPROCESS_TESTWEIGHTSPRODUCER_H
#define POINTPROCESS_TESTWEIGHTSPRODUCER_H
#include "../pointprocess/WeightsProducer.h"
#include <Eigen/Core>

bool testAlphaOne(){
    std::vector<double> target_distances = {2.0, 1.5, 1.0, 0.5, 0.0};
    auto wp = WeightsProducer(1.0);
    Eigen::VectorXd res = wp.produce(target_distances);
    Eigen::VectorXd expected(5);
    expected << 1.0,1.0,1.0,1.0,1.0;
    return (expected - res).norm() < 1e-10;
}

bool testAlphaDot9(){
    std::vector<double> target_distances = {2.0, 1.5, 1.0, 0.5, 0.0};
    auto wp = WeightsProducer(0.9);
    Eigen::VectorXd res = wp.produce(target_distances);
    Eigen::VectorXd expected(5);
    expected << 0.897507,0.946056,0.997223,1.051117,1.108003;
    return (expected - res).norm() < 1e-4;
}


#endif //POINTPROCESS_TESTWEIGHTSPRODUCER_H
