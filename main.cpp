#include "pointprocess/InterEventDistributions.h"
#include "pointprocess/RegressionPipeline.h"
#include "pointprocess/WeightsProducer.h"
#include "tests/testData.h"

#include <vector>
#include "pointprocess/serialize.h"


int main()
{
    auto td = getTestData();
    auto pip = RegressionPipeline(PointProcessDistributions::LogNormal, 9, true);
    auto fullRes = pip.fullRegression(td.testEvents,60.0,0.005,true,20000, WeightsProducer(0.98));

    ppRes2csv(fullRes, std::string("../../notebooks/TEST/LogNormalData.csv"), std::string("../../notebooks/TEST/LogNormalTaus.csv"));
}
