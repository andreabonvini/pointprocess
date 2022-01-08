#include <vector>
#include "../pointprocess/PointProcessUtils.h"
#include "../pointprocess/InterEventDistributions.h"
#include "../pointprocess/RegressionPipeline.h"
#include "../pointprocess/WeightsProducer.h"
#include "../external/csv.h"

// TODO: include it directly as a .csv...
#include "../../tests/data/testData.h"



int main() {
    auto td = getTestData();
    auto pip = RegressionPipeline(
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
            );
    pp::utils::serialize::ppResData2csv(ppRes, std::string("myData.csv"));
    pp::utils::serialize::ppResTaus2csv(ppRes, std::string("myTaus.csv"));
}

