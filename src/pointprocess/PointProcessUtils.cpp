//
// Created by Andrea Bonvini on 23/12/21.
//

#include "PointProcessUtils.h"


pp::Stats pp::utils::computeStats(std::vector<double> &taus) {
    std::vector<double> rescaledTimes;
    Eigen::VectorXd z(taus.size());
    for (long i = 0; i < taus.size(); i++) {
        z[i] = taus[i];
    }
    z = -z;
    z = 1.0 - z.array().exp();
    std::sort(z.data(), z.data() + z.size());
    auto lin = Eigen::VectorXd::LinSpaced(z.size(), 0.0, 1.0);
    auto lu = Eigen::VectorXd::LinSpaced(z.size(), 1.36 / sqrt(z.size()), 1.0 + 1.36 / sqrt(z.size()));
    auto ll = Eigen::VectorXd::LinSpaced(z.size(), -1.36 / sqrt(z.size()), 1.0 - 1.36 / sqrt(z.size()));
    double ksDistance = (z.array() - lin.array()).abs().maxCoeff() / sqrt(2.0);
    double percOut = 0.0;
    double autoCorr = 0.0; // TODO: IMPLEMENT!
    for (long i = 0; i < z.size(); i++) {
        percOut += (double) z[i] < ll[i] || z[i] > lu[i];
    }
    percOut = percOut / ((double) z.size());
    return Stats(ksDistance, percOut, autoCorr);

}


pp::PipelineSetup
pp::utils::getPipelineSetup(const std::vector<double> &events, bool hasTheta0_, unsigned char AR_ORDER_,
                 double windowLength, double delta, WeightsProducer weightsProducer) {

    // Consistency check
    if (events[events.size() - 1] < windowLength) {
        throw std::invalid_argument("The window length is too wide.");
    }
    if (fabs(std::remainder(windowLength,delta)) > 1e-10){
        // no better way to check if windowLength is a multiple of delta...
        throw std::invalid_argument("windowLength should be a multiple of delta.");
    }
    // Find the index of the last event within window_length
    unsigned long last_event_index = 0;
    for (unsigned long index = 0; index < events.size(); index++) {
        if (events[index] > windowLength) {
            last_event_index = index - 1;
            break;
        }
    }
    // Find total number of time bins
    auto bins = (unsigned long) std::ceil(events[events.size() - 1] / delta);
    auto bins_in_window = (unsigned long) (windowLength / delta);

    return pp::PipelineSetup(delta, events, hasTheta0_, AR_ORDER_, last_event_index, bins, bins_in_window,
            weightsProducer);
}

void pp::utils::computeTaus(std::vector<double> &taus, const std::vector<double> &lambdas, const PipelineSetup &setup) {

    unsigned long last_event_index = setup.last_event_index;
    double intL = 0.0;
    double pr = 0.0;
    bool wait = true;

    unsigned long offset = setup.bins_in_window;
    for (unsigned long bin_index = offset; bin_index <= setup.bins; bin_index++) {
        double currentTime = (double) bin_index * setup.delta;
        bool eventHappened = setup.events[last_event_index + 1] <= currentTime;
        if (eventHappened) {
            last_event_index++;
            double dt = setup.events[last_event_index] - (currentTime - setup.delta);
            intL = intL + dt * pr;
            if (!wait) {
                taus.push_back(intL);
            }
            // we set wait to false when we observe the FIRST event after the starting timeWindow, this way we'll
            // be able to compute the integral of lambda (i.e. tau) only for fully observed inter-event intervals.
            /* e.g.
             *              1  2    3   4    5    6  7
             * events:    --|--|----|---|----|----|--|
             *   bins:    ----------------------------
             * timeWindow:  \_____________/XX       <- WRONG
             *                 \_____________/      <- RIGHT
             *  We have to wait the first observable event (5) in order to compute a reliable estimate
             *  of the first tau (integral of lambda)
             */
            wait = false;
            intL = 0.0;
        } else {
            intL = intL + setup.delta * lambdas[bin_index - offset];
        }
        if (bin_index <= lambdas.size()) {
            pr = lambdas[bin_index - offset];
        }
    }
}



void pp::utils::serialize::ppResData2csv(Result &ppRes, const std::string &outputResultsName) {
    // The main Regression Results will be saved at outputResultsName
    std::ofstream csv(outputResultsName.c_str());
    // Define header...
    csv << ((ppRes.distribution == PointProcessDistributions::InverseGaussian)
            ? "Time,Mu,Sigma,Kappa,Lambda,eventHappened,nIter,Likelihood,maxGrad,meanInterval,Theta0"
            : "Time,Mu,Sigma,Lambda,eventHappened,nIter,Likelihood,maxGrad,meanInterval,Theta0");
    for (long i = 0; i < ppRes.AR_ORDER; i++) {
        csv << "," << "Theta" << std::to_string(i + 1);
    }
    csv << "\n";

    if (ppRes.distribution == PointProcessDistributions::InverseGaussian) {
        // In case we are serializing an InverseGaussian Regression Result we have to save also the Kappa parameter for
        // each time step.
        for (auto &res: ppRes.results) {
            auto tmp = dynamic_cast<IGRegressionResult *>(res.get());
            csv << ppRes.t0 + tmp->time << "," << tmp->mu << "," << tmp->sigma << "," << tmp->kappa << ","
                << tmp->lambda << ","
                << tmp->eventHappened << "," << tmp->nIter << "," << tmp->likelihood << "," << tmp->maxGrad
                << "," << tmp->meanInterval << "," << tmp->theta0;
            for (long i = 0; i < tmp->thetaP.size(); i++) {
                csv << "," << tmp->thetaP(i);
            }
            csv << "\n";
        }
    } else {
        for (auto &res: ppRes.results) {
            csv << ppRes.t0 + res->time << "," << res->mu << "," << res->sigma << "," << res->lambda << ","
                << res->eventHappened << "," << res->nIter << "," << res->likelihood << "," << res->maxGrad
                << res->meanInterval << "," << res->theta0;
            for (long i = 0; i < res->thetaP.size(); i++) {
                csv << "," << res->thetaP(i);
            }
            csv << "\n";
        }
    }

    csv.close();
}


void pp::utils::serialize::ppResTaus2csv(Result &ppRes, const std::string &outputTausName) {
    std::ofstream taus(outputTausName.c_str());
    taus << "Taus\n";
    for (auto &tau: ppRes.taus) {
        taus << tau << "\n";
    }
    taus.close();
}


void pp::utils::logging::printProgress(double currentTime, double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s] (Current time: %f s)", val, lpad, PBSTR, rpad, "", currentTime);
    fflush(stdout);
}
