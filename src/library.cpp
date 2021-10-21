//
// Created by Andrea Bonvini on 18/05/21.
//

#include "../pointprocess/optimizers/BaseOptimizer.h"
#include "../pointprocess/optimizers/InverseGaussianOptimizer.h"
#include "../pointprocess/optimizers/LogNormalOptimizer.h"
#include "../pointprocess/optimizers/GaussianOptimizer.h"
#include "../pointprocess/serialize.h"


//#include "../markedpointprocess/optimizers/BaseOptimizerMarked.h"
//#include "../markedpointprocess/optimizers/InverseGaussianOptimizerMarked.h"
//#include "../markedpointprocess/optimizers/LogNormalOptimizerMarked.h"
//#include "../markedpointprocess/optimizers/GaussianOptimizerMarked.h"
//#include "../markedpointprocess/RegressionPipelineMarked.h"
//#include "../markedpointprocess/serializeMarked.h"



extern "C" void regrlikel(
        unsigned int n_events,
        double c_events_pointer[n_events],
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

    std::vector<double> events;
    for (unsigned int i = 0; i < n_events; i++){
        events.push_back(c_events_pointer[i]);
    }
    auto pip = RegressionPipeline(distribution, AR_ORDER, hasTheta0);
    auto ppRes = pip.fullRegression(events,windowLength,delta,rightCensoring,maxIter, WeightsProducer(alpha));
    if (serializeData) {
        ppResData2csv(ppRes, std::string(outputDataName));
    }

    ppResTaus2csv(ppRes, std::string(outputTausName));
}



//extern "C" void regrlikel_marked(
//        unsigned int n_events,
//        double c_events_pointer[n_events],
//        double c_amplitudes_pointer[n_events],
//        double windowLength,
//        double delta,
//        unsigned char AR_ORDER,
//        bool hasTheta0,
//        double alpha,
//        double beta,
//        PointProcessDistributions distribution,
//        unsigned int maxIter,
//        bool serializeData,
//        char * outputDataName,
//        char * outputTausName
//){
//
//    std::vector<double> events;
//    for (unsigned int i = 0; i < n_events; i++){
//        events.push_back(c_events_pointer[i]);
//    }
//    std::vector<double> amplitudes;
//    for (unsigned int i = 0; i < n_events; i++){
//        amplitudes.push_back(c_amplitudes_pointer[i]);
//    }
//    auto pipMark = RegressionPipelineMarked(distribution, AR_ORDER, hasTheta0);
//    auto ppRes = pipMark.fullRegressionMarked(events,amplitudes,windowLength,delta,maxIter, WeightsProducer(alpha), WeightsProducer(beta));
//
//    if (serializeData) {
//        ppResData2csvMarked(ppRes, std::string(outputDataName));
//    }
//    // ppResTaus2csv(ppRes, std::string(outputTausName));
//}


