//
// Created by Andrea Bonvini on 18/05/21.
//

#include "pointprocess/OptimizersFactory.h"
#include "pointprocess/PointProcessUtils.h"
#include "pointprocess/RegressionPipeline.h"


extern "C" void regrlikel(
        unsigned int n_events,
        double c_events_pointer[],
        double windowLength,
        double delta,
        unsigned char AR_ORDER,
        bool hasTheta0,
        bool rightCensoring,
        double alpha,
        PointProcessDistributions distribution,
        unsigned int maxIter,
        bool serializeData,
        char * outputDataName,
        char * outputTausName
){

/*    auto pip = RegressionPipeline(
            PointProcessDistributions::InverseGaussian, // distribution
            9, // AR_ORDER
            true // hasTheta0
    );
    auto ppRes = pip.fullRegression(
            td.testEvents, // event times
            60.0,  // windowLength
            0.005, // delta
            true,  // rightCensoring
            1000,  // maxIter
            WeightsProducer(0.02)
    );*/

    std::vector<double> events;
    events.reserve(n_events);
    for (unsigned int i = 0; i < n_events; i++){
        events.push_back(c_events_pointer[i]);
    }
    auto pip = RegressionPipeline(distribution, AR_ORDER, hasTheta0);
    auto ppRes = pip.fullRegression(
            events,
            windowLength,
            delta,
            rightCensoring,
            maxIter,
            WeightsProducer(alpha)
            );

    if (serializeData) {
        pp::utils::serialize::ppResData2csv(ppRes, std::string(outputDataName));
    }

    pp::utils::serialize::ppResTaus2csv(ppRes, std::string(outputTausName));
}



