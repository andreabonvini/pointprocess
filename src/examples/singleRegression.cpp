//
// Created by Andrea Bonvini on 19/05/21.
//

#include "../pointprocess/InterEventDistributions.h"
#include "../pointprocess/optimizers/InverseGaussianOptimizer.h"
#include "../pointprocess/WeightsProducer.h"
#include "../pointprocess/serialize.h"

// TODO: include it directly as a .csv...
#include "../../tests/data/testData.h"



int main() {
    // Retrieve some data...
    auto td = getTestData();
    // Since td.testEvents is a std::vector<double>, in this case we need to copy it inside a std::deque<double>
    // in order to build a PointProcessDataset.
    std::deque<double> events;
    for (int i = 75; i < 300; i++){
        events.push_back(td.testEvents[i]);
    }

    auto wp = WeightsProducer(0.02);

    auto dataset = PointProcessDataset::load(
            events, // event_times
            0, // AR_ORDER
            true, // hasTheta0
            wp
    );

    auto optimizer = InverseGaussianOptimizer();
    auto result = optimizer.singleRegression(
            dataset,
            false, // rightCensoring
            10000 // maxIter
    );

    std::cout << "Mu: " << result->mu << std::endl;
    // Do stuff...
}