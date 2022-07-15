#ifndef POINTPROCESS_OPTIMIZERSFACTORY_H
#define POINTPROCESS_OPTIMIZERSFACTORY_H


#include <utility>
#include <unordered_map>

#include "optimizers/BaseOptimizer.h"
#include "optimizers/InverseGaussianOptimizer.h"
#include "optimizers/GaussianOptimizer.h"
#include "optimizers/LogNormalOptimizer.h"

//Factory design pattern (https://hackernoon.com/desing-patterns-exploring-factory-method-in-modern-c-hi1h3uvw)
class OptimizersFactory {
    std::unordered_map<pointprocess::Distributions, std::function<std::unique_ptr<BaseOptimizer>() >>  optimizer_factories_map;

public:
    OptimizersFactory();
    std::unique_ptr<BaseOptimizer> create(pointprocess::Distributions dist);
};

#endif //POINTPROCESS_OPTIMIZERSFACTORY_H