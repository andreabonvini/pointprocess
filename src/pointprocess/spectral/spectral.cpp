//
// Created by Andrea Bonvini on 11/06/22.
//

#include "spectral.h"
#include <unsupported/Eigen/Polynomials>
#include <boost/math/special_functions/erf.hpp>
#include <iostream>
#include <algorithm>
#include <numbers>
#include <utility>

static constexpr int N_SAMPLES = 2048;

// LCOV_EXCL_START
pointprocess::spectral::Pole::Pole(std::complex<double> pos_, double frequency_, double power_, std::complex<double> residual_){
    pos = pos_;
    frequency = frequency_;
    power = power_;
    residual = residual_;
}
// LCOV_EXCL_STOP


// LCOV_EXCL_START
pointprocess::spectral::SpectralAnalysis::SpectralAnalysis(
        Eigen::VectorXd frequencies_,
        Eigen::VectorXd powers_,
        std::vector<pointprocess::spectral::Pole>& poles_,
        std::vector<Eigen::VectorXcd> comps_
) : frequencies(std::move(frequencies_)), powers(std::move(powers_)), poles(poles_), comps(std::move(comps_)){}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
pointprocess::spectral::HeartRateVariabilityIndices::HeartRateVariabilityIndices(
        double powVLF_,
        double powLF_,
        double powHF_
){
    powVLF = powVLF_;
    powLF = powLF_;
    powHF = powHF_;
}
// LCOV_EXCL_STOP


// LCOV_EXCL_START
pointprocess::spectral::HeartRateVariabilityIndices::HeartRateVariabilityIndices() {
    powLF = -1.0;
    powVLF = -1.0;
    powHF = -1.0;

}
// LCOV_EXCL_STOP


pointprocess::spectral::HeartRateVariabilityIndices
pointprocess::spectral::computeHeartRateVariabilityIndices(
        std::vector<Pole>& poles
){
    double powVLF = 0.0;
    double powLF = 0.0;
    double powHF = 0.0;
    for (auto& pole : poles){
        if (std::abs(pole.frequency) <= 0.04 && pole.power > 0.0){
            powVLF += pole.power;
        }
        else if  (std::abs(pole.frequency) > 0.04 && std::abs(pole.frequency) <= 0.15 && pole.power > 0.0){
            powLF += pole.power;
        }
        else if (std::abs(pole.frequency) > 0.15 && std::abs(pole.frequency) <= 0.45 && pole.power > 0.0){
            powHF += pole.power;
        } else if (pole.power <= 0.0){
            // FIXME: What should we do in this case?
            // std::cout << "WARNING: negative power detected.\n";
        }
    }

    return {powVLF,powLF,powHF};
}


bool pointprocess::spectral::areConj(
        const std::complex<double>& c1,
        const std::complex<double>& c2,
        double tol
        ){
    return(
            (std::abs(std::real(c1) - std::real(c2)) < tol)
            &&
            (std::abs(std::imag(c1) + std::imag(c2)) < tol)
    );
}



// TODO: test!
// LCOV_EXCL_START
std::vector<pointprocess::spectral::Pole> pointprocess::spectral::computePoles(
        Eigen::VectorXd& thetaP, double var, double fSamp
        ){

    auto one = Eigen::VectorXd::Ones((long) 1);
    Eigen::VectorXd ar(thetaP.size() + 1);
    ar << one, - thetaP;  // ar <- [1, -θ1, -θ2, ..., θp]

    // ================ Compute poles complex values ================
    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    solver.compute(ar.reverse());
    const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &r = solver.roots();

    std::vector<std::complex<double>> polesValues;
    for(int i = 0; i < r.size(); i++){
        polesValues.push_back(r[i]);
    }

    // Sort the roots by absolute phase angle.
    std::sort(polesValues.begin(), polesValues.end(),
              [](
                      const std::complex<double>& a,
                      const std::complex<double>& b
              ) {
                  return std::abs(std::arg(a)) < std::abs(std::arg(b));
              });

    // Fix AR models that might have become slightly unstable due to the estimation process
    // using an exponential decay (see Stoica and Moses, Signal Processing 26(1) 1992).
    // Matches the MATLAB reference implementation by Barbieri & Citi.
    // TODO: Review — validate numerically against MATLAB reference for edge cases.
    double maxAbsPole = 0.0;
    for (const auto& root : polesValues){
        maxAbsPole = std::max(maxAbsPole, std::abs(root));
    }
    double modScale = std::min(0.99 / maxAbsPole, 1.0);
    for(auto& root : polesValues){
        root *= modScale;
    }
    // Apply cumulative product scaling to thetaP: thetaP[i] *= modScale^(i+1)
    double cumProd = 1.0;
    for(Eigen::Index i = 0; i < thetaP.size(); i++){
        cumProd *= modScale;
        thetaP(i) *= cumProd;
    }

    // Compute poles residuals.

    std::vector<std::complex<double>> polesResiduals(polesValues.size());
    for (size_t i = 0; i < polesValues.size(); i++){
        std::complex<double> prod1(1.0,0.0);
        std::complex<double> prod2(1.0,0.0);
        for (size_t j = 0; j < polesValues.size(); j++){
            if (j!=i)
                prod1 *= polesValues[i] - polesValues[j];
            prod2 *= 1.0/polesValues[i] - std::conj(polesValues[j]);
        }

        polesResiduals[i] = 1.0 / (
                polesValues[i]
                * prod1 * prod2
        );
    }

    // Compute poles frequencies.

    std::vector<double> polesFrequencies(polesValues.size());
    for(size_t i = 0; i < polesFrequencies.size(); i++){
        polesFrequencies[i] = std::arg(polesValues[i]) / (2.0 * std::numbers::pi) * fSamp;
    }

    // Compute poles powers.

    std::vector<double> polesPowers(polesValues.size());
    for(size_t i = 0; i < polesPowers.size(); i++){
        polesPowers[i] = var * std::real(polesResiduals[i]);
    }

    std::vector<pointprocess::spectral::Pole> poles;
    poles.reserve(polesValues.size());
    for (size_t i = 0; i < polesValues.size(); i++){
        poles.emplace_back(polesValues[i], polesFrequencies[i], polesPowers[i], polesResiduals[i]);
    }
    return poles;
}
// LCOV_EXCL_STOP


// TODO: test!
// LCOV_EXCL_START
pointprocess::spectral::SpectralAnalysis
pointprocess::spectral::computeSpectralAnalysis(
        Eigen::VectorXd& thetaP,
        double meanInterval,
        double variance,
        bool aggregate) {

    // Work on a local copy: computePoles may modify thetaP via the stability fix,
    // and we don't want to mutate the caller's data.
    Eigen::VectorXd localThetaP = thetaP;

    double var = variance * static_cast<double>(1e6); // from [s^2] to [ms^2]
    double fSamp = 1.0 / meanInterval;

    auto poles = computePoles(localThetaP, var, fSamp);

    // Build AR polynomial from (possibly stabilized) thetaP
    auto one = Eigen::VectorXd::Ones((long) 1);
    Eigen::VectorXd ar(localThetaP.size() + 1);
    ar << one, - localThetaP;

    auto fs = Eigen::VectorXcd ::LinSpaced(N_SAMPLES, -0.5, 0.5);
    // z = e^(-2πfT)
    // z: unit delay operator
    auto z = (fs.array() * (2.0 * std::complex(0.0,1.0) *  std::numbers::pi)).exp();
    // σ^2 : Sample variance
    // T: sampling interval
    // Power(z) = (σ^2*T)/ |1+θ1*z^(-1)+...+θp*z^(-p)|^2
    // powers = (var / fsamp) / abs(np.polyval(ar, np.conj(z))) ** 2
    Eigen::VectorXd powers = (var/fSamp) /  z.unaryExpr(
            [&ar](const std::complex<double>& z_el){
                return std::abs(Eigen::poly_eval(ar, std::conj(z_el)));
            }
    ).array().pow(2.0);

    // We also save the spectral components for each frequency value for each pole))
    // std::vector<std::array<std::complex<double>, N_SAMPLES>> polesComps(polesValues.size());
    std::vector<Eigen::VectorXcd> polesComps(poles.size());
    std::vector<std::complex<double>> refPoles(poles.size());
    for(size_t i = 0; i < refPoles.size(); i++){
        refPoles[i] = 1.0 / std::conj(poles[i].pos);
    }

    for(size_t i = 0; i < poles.size(); i++){
        std::vector<std::complex<double>> pp(z.size());
        for(Eigen::Index j = 0; j < z.size(); j++) {
            pp[j] = poles[i].residual * poles[i].pos / (z[j] -  poles[i].pos);
        }
        std::vector<std::complex<double>> refpp(z.size());
        for(Eigen::Index j = 0; j < z.size(); j++) {
            refpp[j] = -std::conj(poles[i].residual) * refPoles[i] / (z[j] -  refPoles[i]);
        }

        polesComps[i].resize(N_SAMPLES);
        for(Eigen::Index j = 0; j < z.size(); j++) {
            polesComps[i](j) = var / fSamp * (pp[j] + refpp[j]);
        }
    }

    // Aggregate complex conjugate poles in poles_comps_agg
    std::vector<Eigen::VectorXcd> polesCompsAgg;
    polesCompsAgg.push_back(polesComps[0]);

    size_t current_size_ = 0;
    for(size_t i = 1; i < poles.size(); i++){
        if(areConj(poles[i].pos, poles[i-1].pos, 1e-5)){
            for(Eigen::Index j = 0; j < polesCompsAgg[current_size_].size(); j++){
                polesCompsAgg[current_size_](j) += polesComps[i](j);
            }
        }
        else{
            polesCompsAgg.push_back(polesComps[i]);
            current_size_++;
        }
    }

    Eigen::VectorXd frequencies(fs.size());
    frequencies = fs.array().real() * fSamp;

    return {
            frequencies,
            powers,
            poles,
            (aggregate ? polesCompsAgg : polesComps)
    };
}
// LCOV_EXCL_STOP


Eigen::VectorXd pointprocess::spectral::hamming(int n){
    Eigen::VectorXd w(n+1);
    for (int i = 0; i <= n; i++){
        w(i) = 0.54 - 0.46*std::cos(2.0*std::numbers::pi*static_cast<float>(i)/n);
    }
    return w;
}


Eigen::VectorXd pointprocess::spectral::filter1D(Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& a){
    // TODO: probably not so efficient for huge x vectors, optimize.
    // TODO: only the case where a = [1] has been implemented so far.
    // https://it.mathworks.com/help/matlab/data_analysis/filtering-data.html
    // a(1)y(n)=b(1)x(n)+b(2)x(n−1)+…+b(N_b)x(n−N_b + 1)−a(2)y(n−1)−…−a(N_a)y(n−N_a + 1)
    /* e.g.
     * b = (1,2,3)
     * x = (10,20,30,40,50,60)
     * y = (
     *      1*10               = 10
     *      1*20 + 2*10        = 40
     *      1*30 + 2*20 + 3*10 = 100
     *      1*40 + 2*30 + 3*20 = 160
     *      1*50 + 2*40 + 3*30 = 220
     *      1*60 + 2*50 + 3*40 = 280
     *
     *
     *  X = 10 00 00
     *      20 10 00
     *      30 20 10
     *      40 30 20
     *      50 40 30
     *      60 50 40
     *
            )
    */
    Eigen::VectorXd result(x.size());
    if (a.size() == 1){
        Eigen::MatrixXd X = Eigen::MatrixXd::Zero(x.size(), b.size());
        for(int i =0; i < x.size(); i++){
            for(int j=0; j < b.size(); j++){
                if((i-j)>=0)
                    X(i,j) = x(i-j);
            }
        }
        result = X * b;
    }

    return result;

}


Eigen::VectorXd pointprocess::spectral::computeAutoCorrelation(std::vector<double>& taus, int maxlag){
    // Transform taus to uniform on (0,1]: Z = 1 - exp(-taus)
    // Then map to standard normal via inverse error function: N = erfinv(2*Z - 1)
    // Then compute autocorrelation of (N - mean(N)) for lags 1..maxlag.

    const auto n = static_cast<Eigen::Index>(taus.size());
    constexpr double small = 1e-5;

    Eigen::VectorXd N(n);
    for (Eigen::Index i = 0; i < n; i++){
        double z = 1.0 - std::exp(-taus[static_cast<size_t>(i)]);
        z = std::clamp(z, small, 1.0 - small);
        N(i) = boost::math::erf_inv(2.0 * z - 1.0);
    }

    double mean = N.mean();
    Eigen::VectorXd centered = N.array() - mean;
    double variance = centered.squaredNorm(); // sum of squares (unnormalized)

    // Direct autocorrelation computation: O(n * maxlag), efficient for typical maxlag << n.
    Eigen::VectorXd autoCorr(maxlag);
    for (int lag = 1; lag <= maxlag; lag++){
        double sum = 0.0;
        for (Eigen::Index i = 0; i < n - lag; i++){
            sum += centered(i) * centered(i + lag);
        }
        autoCorr(lag - 1) = sum / variance;
    }

    return autoCorr;
}

