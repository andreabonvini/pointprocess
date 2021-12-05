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
        explicit WeightsProducer(double alpha_ = 0.02){
            alpha = alpha_;
        }
        //
        [[nodiscard]] Eigen::VectorXd produce(std::vector<double> target_distances) const{
            Eigen::VectorXd eta = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(target_distances.data(), (long) target_distances.size());
            // Paper's implementation below: w(t-u) = exp(-alpha * (t-u))
            // IMPORTANT: Differently from the paper in this case we will use a weight equal to 1.0 for the right-censoring part
            // (the paper used the same weight of the last observed event).
            eta = Eigen::exp(- alpha * eta.array());
            // MATLAB implementation below:
            // eta = Eigen::exp(eta.array() * log(alpha));
            // eta = eta / eta.sum() * target_distances.size();
            return eta;
        }
};

#endif //POINTPROCESS_WEIGHTSPRODUCER_H
