//
// Created by Andrea Bonvini on 15/07/22.
//

#include "OptimizersFactory.h"

OptimizersFactory::OptimizersFactory() {
    optimizer_factories_map[pointprocess::Distributions::Gaussian] = [] { return std::make_unique<GaussianOptimizer>(); };
    optimizer_factories_map[pointprocess::Distributions::InverseGaussian] = [] { return std::make_unique<InverseGaussianOptimizer>(); };
    optimizer_factories_map[pointprocess::Distributions::LogNormal] = [] { return std::make_unique<LogNormalOptimizer>(); };
}
std::unique_ptr<BaseOptimizer> OptimizersFactory::create(pointprocess::Distributions dist) { return optimizer_factories_map[dist](); }


