//
// Created by Andrea Bonvini on 30/05/21.
//

#ifndef POINTPROCESS_TESTPIPELINE_H
#define POINTPROCESS_TESTPIPELINE_H

#include "../pointprocess/InterEventDistributions.h"
#include "../pointprocess/RegressionPipeline.h"
#include "../pointprocess/WeightsProducer.h"
#include <vector>


/*
 * TODO:
 *  - Factorize...
 *  - Mock
 */


bool testPipelineSetup(){

    std::vector<double> testEvents = {0.0,0.997,2.0,3.0,4.0,4.997,6.0};
    double windowLength = 5.0;
    double delta = 0.005;
    unsigned char AR_ORDER = 2;
    bool hasTheta0 = true;
    auto setup = getPipelineSetup(testEvents, false, hasTheta0, AR_ORDER, windowLength, delta, 1, WeightsProducer(1.00));
    return (setup.last_event_index == 5 && setup.bins == 1200 && setup.bins_in_window == 1000);

}

bool testTrain(){

    // TODO: Move to TestOptimizers.h

    std::vector<double> testEvents = {0.0,0.997,2.0,3.0,4.0,4.997,6.0};
    double windowLength = 5.0;
    double delta = 0.005;
    unsigned char AR_ORDER = 2;
    bool hasTheta0 = true;
    auto setup = getPipelineSetup(testEvents, false, hasTheta0, AR_ORDER, windowLength, delta, 1, WeightsProducer(1.00));
    /*
     * TODO: In this test I'm just testing the currentTime, it should be possible to add debugging
     *  information to the RegressionResults e.g. by adding the current observed_events and dataset used in each
     *  precise time bin.
     * Since we will use a windowLength of 5.0 seconds,
     * In the first window we expect:
     *    currentTime     : 5.0
     *    observed_events : 0.0, 0.997, 2.0, 3.0, 4.0, 4.997
     *    dataset.wt      : 5.0 - 4.997  = 0.003
     * In the last window we expect:
     *    currentTime     : 6.0
     *    observed_events : 2.0, 3.0, 4.0, 4.997, 6.0
     *    dataset.wt      : 6.0 - 6.0  = 0.0
     */

    auto opt = InverseGaussianOptimizer();
    auto results = opt.train(setup);
    return (results[0]->time == 5.0 && results[results.size() - 1]->time == 6.0);
}


bool testTaus(){
    std::vector<double> testEvents = {0.0,0.997,2.0,3.0,4.0,4.997,6.0,7.0,7.883,9.0,10.222};
    double windowLength = 5.0;
    double delta = 0.005;
    unsigned char AR_ORDER = 2;
    bool hasTheta0 = true;
    auto setup = getPipelineSetup(testEvents, false, hasTheta0, AR_ORDER, windowLength, delta, 1, WeightsProducer(1.00));
    auto opt = InverseGaussianOptimizer();
    auto results = opt.train(setup);

    // Here I mock by hand all the hazard-rates (lambda) contained in results.
    for (auto & result : results){
        result->lambda = 1.0;
    }

    /*
     * In this case we have the following fully observed intervals:
     * 6.000 to  7.000 --> Tau = 1.0 * dt = 1.0 * 1.000 = 1.000
     * 7.000 to  7.883 --> Tau = 1.0 * dt = 1.0 * 0.883 = 0.883
     * 7.885 to  9.000 --> Tau = 1.0 * dt = 1.0 * 1.115 = 1.115
     * 9.000 to 10.222 --> Tau = 1.0 * dt = 1.0 * 1.222 = 1.222
     * Note that we do not include the interval between 4.997 to 6.0 since our moving window starts at 5.000
     * (the interval is not fully observed)
     * Moreover we ignore the leftmost segment of each observed interval since the hazard-rate for that segment will
     * always be 0.0 (e.g we ignore the segment between 7.885 and 7.885)
     */
    std::vector<double> expected_taus = {1.000, 0.883, 1.115, 1.222};
    std::vector<double> taus;
    computeTaus(taus, results, setup);
    return ( fabs(taus[0] - 1.000) < 1e-5 && fabs(taus[1] - 0.883) < 1e-5 && fabs(taus[2] - 1.115) < 1e-5 && fabs(taus[3] - 1.222) < 1e-5 );

}


#endif //POINTPROCESS_TESTPIPELINE_H
