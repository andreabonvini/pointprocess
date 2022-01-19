//
// Created by Andrea Bonvini on 31/12/21.
//

#ifndef POINTPROCESS_BENCHMARKNEWTON_H
#define POINTPROCESS_BENCHMARKNEWTON_H
#include "benchmark/benchmark.h"
#include "../../src/pointprocess/OptimizersFactory.h"
#include "../data/testData.h"
#include "../../src/pointprocess/PointProcessUtils.h"


static void BM_optimizeNewton(benchmark::State& state) {
    // Perform setup here
    auto td = getTestData();
    bool rightCensoring = true;
    bool hasTheta0 = true;
    unsigned char AR_ORDER = 9;
    double windowLength = 60.0;
    double delta = 0.005;
    unsigned long maxIter = 1000;
    auto weightsProducer = WeightsProducer(0.02);

    // We want the first event to be at time 0 (events = events_times - events_times[0])
    std::vector<double> events = td.testEvents;
    double t0 = td.testEvents[0];
    std::transform(td.testEvents.begin(), td.testEvents.end(), events.begin(),[&](auto& value){ return value - t0;});

    auto setup = pp::utils::getPipelineSetup(events, hasTheta0, AR_ORDER, windowLength, delta, weightsProducer);
    auto datasetBuffer = DatasetBuffer(setup);

    unsigned int numberOfParams = AR_ORDER + hasTheta0 + 1; // AR_ORDER + hasTheta0 + 1
    // Declare some useful variables that will be used in the updateNewton() routine.
    auto gradient = Eigen::VectorXd(numberOfParams);
    auto hessian =  Eigen::MatrixXd(numberOfParams,numberOfParams);
    auto rcGradient = Eigen::VectorXd(numberOfParams);
    auto rcHessian = Eigen::MatrixXd(numberOfParams,numberOfParams);
    auto xold = Eigen::VectorXd(numberOfParams);
    auto alpha = Eigen::VectorXd(numberOfParams);
    auto vars = pp::TmpVars(
            gradient,
            hessian,
            rcGradient,
            rcHessian,
            xold,
            alpha
    );

    auto factory = OptimizersFactory();
    auto igOptimizer = factory.create(PointProcessDistributions::InverseGaussian);
    auto x = Eigen::VectorXd(numberOfParams);
    x << 0.0188767,0.0229873, 1.05662, -0.327049, 0.769105, -0.644893, -0.0240915, 0.140929;

    std::vector<std::shared_ptr<pp::RegressionResult>> results;
    results.reserve(datasetBuffer.size());
    int time_step = 1;


    for (auto _ : state) {
        for (auto [currentTime,eventHappened, resetParameters,dataset] : datasetBuffer) {
            pp::utils::logging::printProgress(currentTime, (double) time_step/ (double) datasetBuffer.size());
            time_step++;
            state.ResumeTiming();
            //  =========================== This code gets timed ======================
            auto result = igOptimizer->optimizeNewton(dataset, rightCensoring, maxIter, x, vars);
            //  =======================================================================
            state.PauseTiming();
            result->time = currentTime; // FIXME: Not really good practice (should be in the constructor of RegressionResult)
            results.push_back(result);
        }
    }

    // Compute some statistics
    unsigned long minNumberOfIterations = maxIter;
    unsigned long maxNumberOfIterations = 0;
    double meanNumberOfIterations = 0.0;
    double totalLikelihood = 0.0;
    double percentageOfConvergence = 0.0;
    double percentageOfCDFisOne = 0.0;
    unsigned int i = 0;
    for (auto& res: results){
        totalLikelihood += res->likelihood;
        if(std::isinf(res->likelihood)){
            i++;
            std::cout << i << std::endl;
        }
        percentageOfConvergence += (double) res->converged;
        percentageOfCDFisOne += (double) res->cdfIsOne;
        meanNumberOfIterations += (double) res->nIter;
        minNumberOfIterations = res->nIter < minNumberOfIterations? res->nIter : minNumberOfIterations;
        maxNumberOfIterations = res->nIter > maxNumberOfIterations? res->nIter : maxNumberOfIterations;
    }
    meanNumberOfIterations = meanNumberOfIterations / (double) results.size();
    percentageOfConvergence = percentageOfConvergence / (double) results.size() * 100.0;
    percentageOfCDFisOne = percentageOfCDFisOne / (double) results.size() * 100.0;
    std::cout << "\n=========================================================================" << std::endl;
    std::cout << "Total Likelihood             |         " << totalLikelihood << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Percentage of convergence    |         " << percentageOfConvergence <<  " %\n";
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Percentage of CDF == 1.0     |         " << percentageOfCDFisOne <<  " %\n";
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Minimum number of iterations |         " << minNumberOfIterations << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Maximum number of iterations |         " << maxNumberOfIterations << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Mean number of iterations    |         " << meanNumberOfIterations << std::endl;
    std::cout << "=========================================================================" << std::endl;
}
// Register the function as a benchmark
BENCHMARK(BM_optimizeNewton)->Iterations(1);


#endif //POINTPROCESS_BENCHMARKNEWTON_H
