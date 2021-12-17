#include <utility>
#include <unordered_map>

//Factory design pattern (https://hackernoon.com/desing-patterns-exploring-factory-method-in-modern-c-hi1h3uvw)
class OptimizersFactory {
    std::unordered_map<PointProcessDistributions, std::function<std::unique_ptr<BaseOptimizer>() >>  optimizer_factories_map;

public:
    OptimizersFactory() {
        optimizer_factories_map[PointProcessDistributions::Gaussian] = [] { return std::make_unique<GaussianOptimizer>(); };
        optimizer_factories_map[PointProcessDistributions::InverseGaussian] = [] { return std::make_unique<InverseGaussianOptimizer>(); };
        optimizer_factories_map[PointProcessDistributions::LogNormal] = [] { return std::make_unique<LogNormalOptimizer>(); };
    }
    std::unique_ptr<BaseOptimizer> create(PointProcessDistributions dist) { return optimizer_factories_map[dist](); }
};
