#include <iostream>
#include "pointprocess/InterEventDistributions.h"
#include "pointprocess/RegressionPipeline.h"
#include "pointprocess/WeightsProducer.h"
#include "tests/testData.h"
#include <Eigen/Core>

#include <memory>
#include <vector>
#include <fstream>
#include "csv.h"

void produceCSV_Ig(PointProcessResult& ppRes, bool ig = false){
//    std::ofstream csv("../../notebooks/INVERSE_GAUSSIAN_pphrv_data.csv");
//    std::ofstream taus("../../notebooks/INVERSE_GAUSSIAN_pphrv_taus.csv");
    std::ofstream csv("../../notebooks/tilt-rest/12734_Python.csv");
    std::ofstream taus("../../notebooks/tilt-rest/12734_taus.csv");

    // Define header
    csv << "Time,Mu,Sigma,Kappa,Lambda,eventHappened,nIter,Theta0";
    auto tmp = dynamic_cast<IGRegressionResult*>(ppRes.results[0].get());
    for(long i = 0 ; i < tmp->thetaP.size(); i++) {
        csv << "," << "Theta" << std::to_string(i + 1);
    }
    csv << "\n";
    for (auto& res: ppRes.results){
        auto tmp = dynamic_cast<IGRegressionResult*>(res.get());
        csv << tmp->time << "," << tmp->mu << "," << tmp->sigma << "," << tmp->kappa << "," << tmp->lambda << "," << tmp->eventHappened << "," << tmp->nIter << "," << tmp->theta0 ;
        for (long i = 0 ; i < tmp->thetaP.size(); i++){
            csv << "," << tmp->thetaP(i);
        }
        csv << "\n";
    }
    taus << "Taus\n";
    for (auto& tau: ppRes.taus){
        taus<< tau << "\n";
    }
    csv.close();
    taus.close();

    std::cout << "Done!" << std::endl;
}


void produceCSV(PointProcessResult& ppRes, bool ig = false){
//    std::ofstream csv("../../notebooks/LOGNORMAL_pphrv_data.csv");
//    std::ofstream taus("../../notebooks/LOGNORMAL_pphrv_taus.csv");

//    std::ofstream csv("../../notebooks/GAUSSIAN_pphrv_data.csv");
//    std::ofstream taus("../../notebooks/GAUSSIAN_pphrv_taus.csv");

    std::ofstream csv("../../notebooks/GAUSSIAN_SYNTHETIC_data.csv");
    std::ofstream taus("../../notebooks/GAUSSIAN_SYNTHETIC_taus.csv");

//      std::ofstream csv("../../notebooks/LOGNORMAL_SYNTHETIC_data.csv");
//      std::ofstream taus("../../notebooks/LOGNORMAL_SYNTHETIC_taus.csv");

    // Define header
    csv << "Time,Mu,Sigma,Lambda,eventHappened,nIter,Theta0";
    for(long i = 0 ; i < ppRes.results[0].get()->thetaP.size(); i++) {
        csv << "," << "Theta" << std::to_string(i + 1);
    }
    csv << "\n";
    for (auto& res: ppRes.results){
        csv << res->time << "," << res->mu << "," << res->sigma << "," << res->lambda << "," << res->eventHappened << "," << res->nIter << "," << res->theta0 ;
        for (long i = 0 ; i < res->thetaP.size(); i++){
            csv << "," << res->thetaP(i);
        }
        csv << "\n";
    }
    taus << "Taus\n";
    for (auto& tau: ppRes.taus){
        taus<< tau << "\n";
    }
    csv.close();
    taus.close();

    std::cout << "Done!" << std::endl;
}


int main()
{
    // 12734
    // 12754
    // 12755

    std::vector<double> events;
//    io::CSVReader<1> in("../DATA/synthetic/LogNormalSynthetic.csv");
    io::CSVReader<1> in("../DATA/synthetic/GaussianSynthetic.csv");
    in.read_header(io::ignore_extra_column, "RR");
    double rr;
    while(in.read_row(rr)){
        events.push_back(rr);
    }
    auto td = getTestData();

    auto pip = RegressionPipeline(PointProcessDistributions::InverseGaussian, 9, true);
    auto fullRes = pip.fullRegression(td.testEvents,60.0,0.005,true,1000, WeightsProducer(0.98));

    // produceCSV(fullRes);
}
