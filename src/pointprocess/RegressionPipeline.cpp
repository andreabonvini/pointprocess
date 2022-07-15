//
// Created by Andrea Bonvini on 21/04/21.
//

#include "RegressionPipeline.h"

// LCOV_EXCL_START
RegressionPipeline::RegressionPipeline(
        pointprocess::Distributions distribution,
        unsigned char AR_ORDER,
        bool hasTheta0
  ) : distribution(distribution), AR_ORDER(AR_ORDER), hasTheta0(hasTheta0) {
}
// LCOV_EXCL_STOP

// TODO: test!
// LCOV_EXCL_START
pointprocess::Result RegressionPipeline::fullRegression(
            const std::vector<double>& events_times,
            double windowLength = 60.0,
            double delta = 0.005,
            bool rightCensoring = true,
            unsigned int maxIter = 1000,
            double alpha = 0.0
      ){
        auto weightsProducer = WeightsProducer(alpha);

        // We want the first event to be at time 0 (events = events_times - events_times[0])
        std::vector<double> events = events_times;
        double t0 = events_times[0];
        std::transform(events_times.begin(), events_times.end(), events.begin(),[&](auto& value){ return value - t0;});

        auto pipelineSetup = pointprocess::utils::getPipelineSetup(events, hasTheta0, AR_ORDER, windowLength, delta, weightsProducer);
        auto datasetBuffer = DatasetBuffer(pipelineSetup);


        auto factory = OptimizersFactory();
        auto optimizer = factory.create(distribution);

        auto results = optimizer->train(datasetBuffer, rightCensoring, maxIter);

        std::vector<double> taus;
        // Creating a new vector containing just the instantaneous lambda values to compute the taus.
        std::vector<double> lambdas;
        lambdas.reserve(results.size());
        std::transform(results.begin(), results.end(), std::back_inserter(lambdas),
                       [](auto& res) { return res->lambda; });

        pointprocess::utils::computeTaus(taus, lambdas, pipelineSetup);
        auto stats = pointprocess::utils::computeStats(taus);
        return {
                results,
                taus,
                this->distribution,
                this->AR_ORDER,
                this-> hasTheta0,
                windowLength,
                delta,
                t0,
                stats
                };
    }
// LCOV_EXCL_STOP
