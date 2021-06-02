//
// Created by Andrea Bonvini on 07/05/21.
//

#ifndef POINTPROCESS_TESTOPTIMIZERS_H
#define POINTPROCESS_TESTOPTIMIZERS_H

#include "testData.h"
#include "../pointprocess/InterEventDistributions.h"
#include "../pointprocess/optimizers/GaussianOptimizer.h"

bool gradientRcIsRight(const PointProcessDataset& dataset, PointProcessDistributions dist){

    // FIXME: Change precision in .isApprox() for the first parameter (Kappa or Sigma), the gradients for this values
    //  are in the order of 1e-6. (Maybe only for the Inverse Gaussian, it is enough to change x[0] in the test.)

    Eigen::VectorXd gradientRc(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd approxGradientRc(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd xStart(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd xRight(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd xLeft(dataset.AR_ORDER + dataset.hasTheta0 + 1);

    // retrieve dummy optimizer
    std::shared_ptr<BaseOptimizer> dummyOpt;
    std::vector<double> events = {};
    auto wp = WeightsProducer(1.0);
    auto set = PipelineSetup(0.005,events,true,true,8,10,100,10,10000,wp);

    if (dist == PointProcessDistributions::Gaussian){
        dummyOpt = std::make_shared<GaussianOptimizer>(GaussianOptimizer());
    }
    else if(dist == PointProcessDistributions::InverseGaussian){
        dummyOpt = std::make_shared<InverseGaussianOptimizer>(InverseGaussianOptimizer());
    }
    else if(dist == PointProcessDistributions::LogNormal){
        dummyOpt = std::make_shared<LogNormalOptimizer>(LogNormalOptimizer());
    }

    dummyOpt->populateStartingPoint(xStart);

    if(dist == PointProcessDistributions::Gaussian){
        xStart[0] = 0.1;
    }

    // Compute exact gradient
    dummyOpt->updateGradientRc(xStart,dataset,gradientRc);

    // Approximate gradient
    double machinePrecision = std::cbrt(1e-15);
    double step;
    double negloglikelRcRight;
    double negloglikelRcLeft;

    xRight = xStart;
    xLeft = xStart;
    for(long i = 0; i < xStart.size(); i++){
        step = machinePrecision * xStart[i];
        // Move right and left for parameter i
        xRight[i] += step;
        xLeft[i] -= step;
        // Compute tmp negloglikel and update gradient through finite difference approximation.
        negloglikelRcRight = dummyOpt->computeLikelRc(xRight, dataset);
        negloglikelRcLeft = dummyOpt->computeLikelRc(xLeft, dataset);
        approxGradientRc[i] = (negloglikelRcRight - negloglikelRcLeft) / (2.0 * step);
        // Reset...
        xRight[i] = xStart[i];
        xLeft[i] = xStart[i];
    }
    // Compare
    return approxGradientRc.isApprox(gradientRc,1e-5);
}

bool gradientIsRight(const PointProcessDataset& dataset, PointProcessDistributions dist){

    Eigen::VectorXd gradient(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd approxGradient(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd xStart(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd xRight(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd xLeft(dataset.AR_ORDER + dataset.hasTheta0 + 1);

    // retrieve dummy optimizer
    std::shared_ptr<BaseOptimizer> dummyOpt;
    std::vector<double> events = {};
    auto wp = WeightsProducer(1.0);
    auto set = PipelineSetup(0.005,events,true,true,8,10,100,10,10000,wp);
    if (dist == PointProcessDistributions::Gaussian){
        dummyOpt = std::make_shared<GaussianOptimizer>(GaussianOptimizer());
    }
    else if(dist == PointProcessDistributions::InverseGaussian){
        dummyOpt = std::make_shared<InverseGaussianOptimizer>(InverseGaussianOptimizer());
    }
    else if(dist == PointProcessDistributions::LogNormal){
        dummyOpt = std::make_shared<LogNormalOptimizer>(LogNormalOptimizer());
    }

    dummyOpt->populateStartingPoint(xStart);

    // Compute exact gradient
    dummyOpt->updateGradient(xStart,dataset,gradient);

    // Approximate gradient
    double machinePrecision = std::cbrt(1e-15);
    double step;
    double negloglikelRight;
    double negloglikelLeft;

    xRight = xStart;
    xLeft = xStart;
    for(long i = 0; i < xStart.size(); i++){
        step = machinePrecision * xStart[i];
        // Move right and left for parameter i
        xRight[i] += step;
        xLeft[i] -= step;
        // Compute tmp negloglikel and update gradient through finite difference approximation.
        negloglikelRight = dummyOpt->computeLikel(xRight, dataset);
        negloglikelLeft = dummyOpt->computeLikel(xLeft, dataset);
        approxGradient[i] = (negloglikelRight - negloglikelLeft) / (2.0 * step);
        // Reset...
        xRight[i] = xStart[i];
        xLeft[i] = xStart[i];
    }
    // Compare
    return approxGradient.isApprox(gradient,1e-5);
}

bool hessianIsRight(const PointProcessDataset& dataset, PointProcessDistributions dist){

    Eigen::MatrixXd hessian(dataset.AR_ORDER + dataset.hasTheta0 + 1, dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::MatrixXd approxHessian(dataset.AR_ORDER + dataset.hasTheta0 + 1, dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd xStart(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    auto gradRight = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    auto gradLeft = Eigen::VectorXd(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd xRight(dataset.AR_ORDER + dataset.hasTheta0 + 1);
    Eigen::VectorXd xLeft(dataset.AR_ORDER + dataset.hasTheta0 + 1);

    // retrieve dummy optimizer
    std::shared_ptr<BaseOptimizer> dummyOpt;
    std::vector<double> events = {};
    auto wp = WeightsProducer(1.0);
    auto set = PipelineSetup(0.005,events,true,true,8,10,100,10,10000,wp);
    if (dist == PointProcessDistributions::Gaussian){
        dummyOpt = std::make_shared<GaussianOptimizer>(GaussianOptimizer());
    }
    else if(dist == PointProcessDistributions::InverseGaussian){
        dummyOpt = std::make_shared<InverseGaussianOptimizer>(InverseGaussianOptimizer());
    }
    else if(dist == PointProcessDistributions::LogNormal){
        dummyOpt = std::make_shared<LogNormalOptimizer>(LogNormalOptimizer());
    }

    dummyOpt->populateStartingPoint(xStart);

    // Compute exact gradient
    dummyOpt->updateHessian(xStart,dataset,hessian);

    // Approximate gradient
    double machinePrecision = std::cbrt(1e-15);
    double step;

    xRight = xStart;
    xLeft = xStart;
    for (long i = 0; i < xStart.size(); i++){
        step = machinePrecision * xStart[i];
        // Move right and left for parameter i
        xRight[i] += step;
        xLeft[i] -= step;
        // Compute tmp gradient and update hessianRc through finite difference approximation.
        dummyOpt->updateGradient(xRight, dataset,gradRight);
        dummyOpt->updateGradient(xLeft, dataset,gradLeft);
        for (long j = 0; j < xStart.size() ; j++){
            approxHessian(i,j) = (gradRight[j] - gradLeft[j]) / (2.0 * step);
        }
        // Reset...
        xRight[i] = xStart[i];
        xLeft[i] = xStart[i];
    }
    // Compare
    return approxHessian.isApprox(hessian,1e-5);
}


bool testInverseGaussianGradient(){
    auto dataset = getTestDataset();
    return gradientIsRight(dataset, PointProcessDistributions::InverseGaussian);
}

bool testGaussianGradient(){
    auto dataset = getTestDataset();
    return gradientIsRight(dataset, PointProcessDistributions::Gaussian);
}

bool testLogNormalGradient(){
    auto dataset = getTestDataset();
    return gradientIsRight(dataset, PointProcessDistributions::LogNormal);
}

bool testInverseGaussianGradientRc(){
    auto dataset = getTestDataset();
    return gradientRcIsRight(dataset, PointProcessDistributions::InverseGaussian);
}

bool testGaussianGradientRc(){
    auto dataset = getTestDataset();
    return gradientRcIsRight(dataset, PointProcessDistributions::Gaussian);
}

bool testLogNormalGradientRc(){
    auto dataset = getTestDataset();
    return gradientRcIsRight(dataset, PointProcessDistributions::LogNormal);
}

bool testInverseGaussianHessian(){
    auto dataset = getTestDataset();
    return hessianIsRight(dataset, PointProcessDistributions::InverseGaussian);
}

bool testGaussianHessian(){
    auto dataset = getTestDataset();
    return hessianIsRight(dataset, PointProcessDistributions::Gaussian);
}

bool testLogNormalHessian(){
    auto dataset = getTestDataset();
    return hessianIsRight(dataset, PointProcessDistributions::LogNormal);
}


#endif //POINTPROCESS_TESTOPTIMIZERS_H
