#include "../../pointprocess/InterEventDistributions.h"
#include "../../pointprocess/RegressionPipeline.h"
#include "../../pointprocess/WeightsProducer.h"
#include "../../tests/testData.h"

#include <vector>
#include "../../pointprocess/serialize.h"



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
            WeightsProducer(0.98)
            );
    ppResData2csv(ppRes, std::string("myData.csv"));
    ppResTaus2csv(ppRes, std::string("myTaus.csv"));
}
