#include "../../pointprocess/InterEventDistributions.h"
#include "../../pointprocess/RegressionPipeline.h"
#include "../../pointprocess/WeightsProducer.h"
#include "../../tests/testData.h"

#include <vector>
#include "../../pointprocess/serialize.h"

#include "../../ext/csv.h"


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
            WeightsProducer(0.05)
            );
    ppResData2csv(ppRes, std::string("myData.csv"));
    ppResTaus2csv(ppRes, std::string("myTaus.csv"));
}
//
//int main(){
//    // The file only contains column "a"
//    io::CSVReader<2> in(std::string("../DATA/propofol/EDA/EDA_eeganes02new4.csv"));
//    in.read_header(io::ignore_missing_column, "Times", "Amplitudes");
//    double time,amp;
//    std::vector<double> events_times;
//    while(in.read_row(time,amp)){
//        events_times.push_back(time);
//    }
//    auto pip = RegressionPipeline(
//            PointProcessDistributions::InverseGaussian, // distribution
//            1, // AR_ORDER
//            true // hasTheta0
//            );
//    auto ppRes = pip.fullRegression(
//            events_times, // event times
//            660.0,  // windowLength
//            0.005, // delta
//            true,  // rightCensoring
//            20,  // maxIter
//            WeightsProducer(0.98)
//            );
//}
