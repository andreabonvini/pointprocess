//
// Created by Andrea Bonvini on 20/04/21.
//

#ifndef POINTPROCESS_WEIGHTSPRODUCER_H
#define POINTPROCESS_WEIGHTSPRODUCER_H

#include <Eigen/Core>
#include <vector>

class WeightsProducer {
    public:
        double alpha;
        explicit WeightsProducer(double alpha_ = 0.98){
            alpha = alpha_;
        }
        [[nodiscard]] Eigen::VectorXd produce(std::vector<double> target_distances) const{
            Eigen::VectorXd eta = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(target_distances.data(), (long) target_distances.size());
            eta = Eigen::exp(eta.array() * log(alpha));
            eta = eta / eta.sum() * target_distances.size();
            return eta;
        }
};

#endif //POINTPROCESS_WEIGHTSPRODUCER_H
