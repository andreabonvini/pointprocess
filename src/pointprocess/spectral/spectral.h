//
// Created by Andrea Bonvini on 11/06/22.
//

#ifndef POINTPROCESS_SPECTRAL_H
#define POINTPROCESS_SPECTRAL_H

#include <complex>
#include <vector>
#include <Eigen/Core>

namespace pointprocess::spectral{
    struct Pole{
        std::complex<double> pos;
        double frequency;
        double power;
        std::complex<double> residual;
        Pole(std::complex<double> pos_, double frequency_, double power_, std::complex<double> residual_);
    };

    struct SpectralAnalysis{
        Eigen::VectorXd frequencies;
        Eigen::VectorXd powers;
        std::vector<Pole> poles;
        std::vector<Eigen::VectorXcd> comps;
        SpectralAnalysis(
                Eigen::VectorXd frequencies_,
                Eigen::VectorXd powers_,
                std::vector<Pole>& poles_,
                std::vector<Eigen::VectorXcd> comps_
        );
    };

    struct HeartRateVariabilityIndices{
        double powVLF = -1.0;
        double powLF = -1.0;
        double powHF = -1.0;
        HeartRateVariabilityIndices();
        HeartRateVariabilityIndices(
                double powVLF_,
                double powLF_,
                double powHF_
        );
    };

    bool areConj(
            const std::complex<double>& c1,
            const std::complex<double>& c2,
            double tol
    );

    std::vector<pointprocess::spectral::Pole> computePoles(
            Eigen::VectorXd& thetaP, double var, double fSamp
    );

    SpectralAnalysis computeSpectralAnalysis(
            Eigen::VectorXd& thetaP,
            double meanInterval,
            double variance,
            bool aggregate
    );

    HeartRateVariabilityIndices computeHeartRateVariabilityIndices(
            std::vector<Pole>& poles
    );

    Eigen::VectorXd hamming(int n);


    Eigen::VectorXd filter1D(Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& a);
}




#endif //POINTPROCESS_SPECTRAL_H
