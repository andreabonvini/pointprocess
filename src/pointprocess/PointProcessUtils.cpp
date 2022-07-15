//
// Created by Andrea Bonvini on 23/12/21.
//

#include "PointProcessUtils.h"

#include <utility>
#include <algorithm>

// LCOV_EXCL_START
pointprocess::Stats::Stats(double ksDistance, double percOut, double autoCorr) : ksDistance(ksDistance), percOut(percOut),
                                                                    autoCorr(autoCorr) {}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
pointprocess::TmpVars::TmpVars(
    Eigen::VectorXd gradient,
    Eigen::MatrixXd hessian,
    Eigen::VectorXd rcGradient,
    Eigen::MatrixXd rcHessian,
    Eigen::VectorXd xold,
    Eigen::VectorXd alpha
        ) :
    gradient(std::move(gradient)),
    hessian(std::move(hessian)),
    rcGradient(std::move(rcGradient)),
    rcHessian(std::move(rcHessian)),
    xold(std::move(xold)),
    alpha(std::move(alpha)){}
// LCOV_EXCL_STOP

pointprocess::PipelineSetup::PipelineSetup(
            double delta,
            std::vector<double> events,
            bool hasTheta0,
            unsigned char AR_ORDER,
            unsigned long last_event_index,
            unsigned long b,
            unsigned long biw,
            WeightsProducer weightsProducer
        ) :
            delta(delta),
            events(std::move(events)),
            hasTheta0(hasTheta0),
            AR_ORDER(AR_ORDER),
            last_event_index(last_event_index),
            bins(b),
            bins_in_window(biw),
            weightsProducer(weightsProducer) {}




// LCOV_EXCL_START
pointprocess::RegressionResult::RegressionResult(
        double theta0,
        Eigen::VectorXd thetaP,
        double mu,
        double sigma,
        double lambda,
        double meanInterval,
        unsigned long nIter,
        double likelihood,
        double maxGrad,
        bool converged,
        bool cdfIsOne,
        bool eventHappened,
        double time
    ) :
        theta0(theta0),
        thetaP(std::move(thetaP)),
        mu(mu),
        sigma(sigma),
        lambda(lambda),
        meanInterval(meanInterval),
        nIter(nIter),
        likelihood(likelihood),
        maxGrad(maxGrad),
        converged(converged),
        cdfIsOne(cdfIsOne),
        eventHappened(eventHappened),
        time(time){}

void pointprocess::RegressionResult::computeHRVIndices() {
    double variance = pow(sigma, 2.0); // The standard deviation (sigma) of a random variable, sample,
    // statistical population, data set, or probability distribution is the square root of its variance.
    double var = variance * static_cast<double>(1e6); // from [s^2] to [ms^2]
    auto poles = spectral::computePoles(thetaP, var, 1.0 / meanInterval);
    hrvIndices = spectral::computeHeartRateVariabilityIndices(poles);
}

// I declare a virtual destructor just to have run-time type information (RTTI), which is needed
// to guarantee polymorphic behaviour.
pointprocess::RegressionResult::~RegressionResult() = default;
// LCOV_EXCL_STOP

// LCOV_EXCL_START
pointprocess::IGRegressionResult::IGRegressionResult(
        double theta0_,
        Eigen::VectorXd thetaP_,
        double mu,
        double sigma,
        double lambda,
        double meanInterval,
        unsigned long nIter,
        double likelihood,
        double maxGrad,
        double kappa,
        bool converged,
        bool cdfIsOne,
        bool eventHappened,
        double time)
    :
    kappa(kappa),
    pointprocess::RegressionResult(theta0_, std::move(thetaP_), mu, sigma, lambda, meanInterval, nIter, likelihood, maxGrad, converged,cdfIsOne,eventHappened, time) {}
// LCOV_EXCL_STOP


// LCOV_EXCL_START
pointprocess::Result::Result(
        std::vector<std::shared_ptr<pointprocess::RegressionResult>> results,
        std::vector<double> taus,
        Distributions distribution,
        unsigned char AR_ORDER,
        bool hasTheta0,
        double windowLength,
        double delta,
        double t0,
        Stats stats
    ) :
    results(std::move(results)),
    taus(std::move(taus)),
    distribution(distribution),
    AR_ORDER(AR_ORDER),
    hasTheta0(hasTheta0),
    windowLength(windowLength),
    delta(delta),
    t0(t0),
    stats(stats) {}
// LCOV_EXCL_STOP


//TODO: test (maybe?)
// LCOV_EXCL_START
void pointprocess::Result::computeHRVIndices(){


    // Original MATLAB code:
    /*
            % warn==1 is OK
            % warn~=1 warning:
            %   mod(warn, 2) is the scale factor used to shrink the poles in case of instability (slightly less than 1 is fine)
            %   bitand(floor(warn), 2) if powLF was negative (increase AR order?)
            %   bitand(floor(warn), 4) if powHF was negative (increase AR order?)
            warn = NaN(1,J);

            for i = 1:J
                if isnan(Thetap(1,i))
                    continue;
                end
                [tot,comps,compsp,f,pole_freq,pole_pow,pole_res,poles,mod_scale] = spectral(Thetap(:,i), Var(i), fsamp(i), [], 0);
                warn(i) = mod_scale;
                pf = abs(pole_freq);
                powLF(i) = sum(pole_pow((pf>0.04) & (pf<0.15)));
                if powLF(i) <= 0
                    powLF(i) = powLF(i-1);
                    warn(i) = warn(i) + 2;
                end
                powHF(i) = sum(pole_pow((pf>0.15) & (pf<0.45)));
                if powHF(i) <= 0
                    powHF(i) = powHF(i-1);
                    warn(i) = warn(i) + 4;
                end
                powVLF(i) = sum(pole_pow(pf<0.04));
                powTot(i) = sum(pole_pow);
            end

            b = ceil(J/10);
            if b > 1
                b = hamming(min(21, b));
                b = b/sum(b);
                powLF = fliplr(filter(b,1,fliplr(filter(b,1,powLF))));
                powHF = fliplr(filter(b,1,fliplr(filter(b,1,powHF))));
                powVLF = fliplr(filter(b,1,fliplr(filter(b,1,powVLF))));
                powTot = fliplr(filter(b,1,fliplr(filter(b,1,powTot))));
            end

            if nargin > 3 && dwsample ~= 1
                powLF = decimate(powLF, dwsample, 'fir');
                powHF = decimate(powHF, dwsample, 'fir');
                powVLF = decimate(powVLF, dwsample, 'fir');
                powTot = decimate(powTot, dwsample, 'fir');
            end

            bal = powLF ./ powHF;

     */
    // Hide cursor
    indicators::show_console_cursor(false);
    indicators::ProgressBar bar{
            indicators::option::BarWidth{65},
            indicators::option::Start{"["},
            indicators::option::Fill{"■"},
            indicators::option::Lead{"■"},
            indicators::option::Remainder{"-"},
            indicators::option::End{" ]"},
            indicators::option::PrefixText("Computing HRV indices: "),
            indicators::option::ShowElapsedTime{true},
            indicators::option::ShowRemainingTime{true},
            indicators::option::ShowPercentage(true),
            indicators::option::ForegroundColor{indicators::Color::green},
            indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}
    };



    Eigen::VectorXd powLFVector(results.size());
    Eigen::VectorXd powHFVector(results.size());

    double lastValidPowLF = 1e-10;
    double lastValidPowHF = 1e-10;
    for(int i = 0; i < results.size(); i++){
        results[i]->computeHRVIndices();
        if (results[i]->hrvIndices.powLF <= 0){
            results[i]->hrvIndices.powLF = lastValidPowLF;
        }
        lastValidPowLF = results[i]->hrvIndices.powLF;

        powHFVector(i) = results[i]->hrvIndices.powLF;

        if (results[i]->hrvIndices.powHF <= 0){
            results[i]->hrvIndices.powHF = lastValidPowHF;
        }
        lastValidPowHF = results[i]->hrvIndices.powHF;

        powHFVector(i) = results[i]->hrvIndices.powHF;


        if (i % static_cast<int>(static_cast<double>(results.size()) / 100.0) == 0 || i == (results.size() - 1)){
            auto prog = static_cast<size_t>(static_cast<double>(i + 1) / static_cast<double>(results.size()) * 100 );
            bar.set_progress(prog);
        }
    }
    /*
    b = ceil(J/10);
    if b > 1
    b = hamming(min(21, b));
    b = b/sum(b);
    powLF = fliplr(filter(b,1,fliplr(filter(b,1,powLF))));
    powHF = fliplr(filter(b,1,fliplr(filter(b,1,powHF))));
    powVLF = fliplr(filter(b,1,fliplr(filter(b,1,powVLF))));
    powTot = fliplr(filter(b,1,fliplr(filter(b,1,powTot))));
    end

    if nargin > 3 && dwsample ~= 1
    powLF = decimate(powLF, dwsample, 'fir');
    powHF = decimate(powHF, dwsample, 'fir');
    powVLF = decimate(powVLF, dwsample, 'fir');
    powTot = decimate(powTot, dwsample, 'fir');
    end
    */
    unsigned int b = std::ceil(static_cast<double>(results.size()) / 10.0);
    Eigen::VectorXd hamming_window;
    if (b>1){
        hamming_window = pointprocess::spectral::hamming(std::min(static_cast<unsigned int>(21),b));
        hamming_window = hamming_window.array() / hamming_window.array().sum();
    }

    hrvIndicesComputed = true;
}

// LCOV_EXCL_STOP


//todo: test (maybe?)
// LCOV_EXCL_START
std::map<std::string, Eigen::MatrixXd> pointprocess::Result::toDict(){
    assert(!results.empty());
    std::map<std::string, Eigen::MatrixXd> out;
    Eigen::MatrixXd eventHappenedVector(results.size(),1);
    Eigen::MatrixXd meanIntervalVector(results.size(),1);
    Eigen::MatrixXd thetaPMatrix(results.size(), results[0]->thetaP.size());
    Eigen::MatrixXd theta0Vector(results.size(),1);
    Eigen::MatrixXd sigmaVector(results.size(),1);
    Eigen::MatrixXd lambdaVector(results.size(),1);
    Eigen::MatrixXd cdfIsOneVector(results.size(),1);
    Eigen::MatrixXd convergedVector(results.size(),1);
    Eigen::MatrixXd nIterVector(results.size(),1);
    Eigen::MatrixXd likelihoodVector(results.size(),1);
    Eigen::MatrixXd maxGradVector(results.size(),1);
    Eigen::MatrixXd timeVector(results.size(),1);
    Eigen::MatrixXd muVector(results.size(),1);
    Eigen::MatrixXd powVLFVector(results.size(),1);
    Eigen::MatrixXd powLFVector(results.size(),1);
    Eigen::MatrixXd powHFVector(results.size(),1);


    for(int i = 0; i < results.size(); i++){

        nIterVector(i) = static_cast<double>(results[i]->nIter);
        eventHappenedVector(i) = results[i]->eventHappened;
        convergedVector(i) = results[i]->converged;
        cdfIsOneVector(i) = results[i]->cdfIsOne;
        timeVector(i) = t0 + results[i]->time;
        muVector(i) = results[i]->mu;
        likelihoodVector(i) = results[i]->likelihood;
        sigmaVector(i) = results[i]->sigma;
        lambdaVector(i) = results[i]->lambda;
        meanIntervalVector(i) = results[i]->meanInterval;
        maxGradVector(i) = results[i]->maxGrad;

        for(int j = 0; j < results[i]->thetaP.size(); j++){
            thetaPMatrix(i, j) = results[i]->thetaP(j);
        }
        if(hasTheta0) {
            theta0Vector(i) = results[i]->lambda;
        }

        if(hrvIndicesComputed){
            powVLFVector(i) = results[i]->hrvIndices.powVLF;
            powLFVector(i) = results[i]->hrvIndices.powLF;
            powHFVector(i) = results[i]->hrvIndices.powHF;
        }

    }

    out["cdf_is_one"] = cdfIsOneVector;
    out["sigma"] = sigmaVector;
    out["lambda"] = lambdaVector;
    out["event_happened"] = eventHappenedVector;
    out["mean_interval"] = meanIntervalVector;
    out["thetap"] = thetaPMatrix;
    out["Likelihood"] = likelihoodVector;
    out["n_iter"] = nIterVector;
    out["converged"] = convergedVector;
    out["max_grad"] = maxGradVector;
    out["Time"] = timeVector;
    out["Mu"] = muVector;

    if (hasTheta0){
        out["theta0"] = theta0Vector;
    }

    if(hrvIndicesComputed){
        out["powVLF"] = powVLFVector;
        out["powLF"] = powLFVector;
        out["powHF"] = powHFVector;
    }
    return out;
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
pointprocess::Stats pointprocess::utils::computeStats(std::vector<double> &taus) {
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
    return {ksDistance, percOut, autoCorr};

}
// LCOV_EXCL_STOP


pointprocess::PipelineSetup pointprocess::utils::getPipelineSetup(const std::vector<double> &events, bool hasTheta0_, unsigned char AR_ORDER_,
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

    return {delta, events, hasTheta0_, AR_ORDER_, last_event_index, bins, bins_in_window,
            weightsProducer};
}

void pointprocess::utils::computeTaus(std::vector<double> &taus, const std::vector<double> &lambdas, const PipelineSetup &setup) {

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



