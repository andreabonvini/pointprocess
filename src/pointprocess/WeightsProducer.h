//
// Created by Andrea Bonvini on 20/04/21.
//

#ifndef POINTPROCESS_WEIGHTSPRODUCER_H
#define POINTPROCESS_WEIGHTSPRODUCER_H

#include <Eigen/Core>
#include <vector>

#include <iostream>

class WeightsProducer {
    public:
        double alpha;
        explicit WeightsProducer(double alpha_ = 0.02);
        [[nodiscard]] Eigen::VectorXd produce(std::vector<double> target_distances) const;
};

#endif //POINTPROCESS_WEIGHTSPRODUCER_H
