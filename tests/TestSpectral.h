//
// Created by Andrea Bonvini on 12/07/22.
//

#ifndef POINTPROCESS_TESTSPECTRAL_H
#define POINTPROCESS_TESTSPECTRAL_H


#include "../src/pointprocess/spectral/spectral.h"
#include <Eigen/Core>
#include <gtest/gtest.h>


namespace pointprocess {
    namespace tests {
        TEST(TestSpectral,TestHamming){
            Eigen::VectorXd expected_hamming_5(6);
            expected_hamming_5 << 0.08, 0.39785218, 0.91214782, 0.91214782, 0.39785218, 0.08;
            Eigen::VectorXd expected_hamming_6(7);
            expected_hamming_6 << 0.08, 0.31, 0.77, 1.  , 0.77, 0.31, 0.08;
            auto result_5 = pointprocess::spectral::hamming(5);
            EXPECT_TRUE(result_5.isApprox(expected_hamming_5,1e-5));
            auto result_6 = pointprocess::spectral::hamming(6);
            EXPECT_TRUE(result_6.isApprox(expected_hamming_6,1e-5));
        }


        TEST(TestSpectral,TestComputeHeartRateVariabilityIndices){
            auto dummyComp = std::complex<double>(0.5,0.5);
            std::vector<pointprocess::spectral::Pole> dummyPoles = {
                    pointprocess::spectral::Pole(dummyComp,0.03,1.0,dummyComp),
                    pointprocess::spectral::Pole(dummyComp,0.04,1.0,dummyComp),
                    pointprocess::spectral::Pole(dummyComp,0.08,2.0,dummyComp),
                    pointprocess::spectral::Pole(dummyComp,0.09,2.0,dummyComp),
                    pointprocess::spectral::Pole(dummyComp,0.20,4.0,dummyComp),
                    pointprocess::spectral::Pole(dummyComp,0.21,4.0,dummyComp)
            };
            auto result = pointprocess::spectral::computeHeartRateVariabilityIndices(dummyPoles);
            EXPECT_EQ(result.powVLF, 2.0);
            EXPECT_EQ(result.powLF, 4.0);
            EXPECT_EQ(result.powHF, 8.0);

        }

        TEST(TestSpectral,TestAreConjTrue1){
            double tol = 1e-3;
            std::complex<double> v1 = {0.500,-0.5005};
            std::complex<double> v2 = {0.500, 0.5000};
            auto result = pointprocess::spectral::areConj(v1,v2,tol);
            EXPECT_TRUE(result);
        }
        TEST(TestSpectral,TestAreConjTrue2){
            double tol = 1e-3;
            std::complex<double> v1 = {0.500,-0.5005};
            std::complex<double> v2 = {0.500, 0.5000};
            auto result = pointprocess::spectral::areConj(v1,v2,tol);
            EXPECT_TRUE(result);
        }
        TEST(TestSpectral,TestAreConjFalse1){
            double tol = 1e-5;
            std::complex<double> v1 = {0.500,-0.5005};
            std::complex<double> v2 = {0.500, 0.5000};
            auto result = pointprocess::spectral::areConj(v1,v2,tol);
            EXPECT_FALSE(result);
        }

        TEST(TestSpectral,TestAreConjFalse2){
        double tol = 1e-3;
            std::complex<double> v1 = {-0.5,-0.5005};
            std::complex<double> v2 = {0.5, 0.5005};
            auto result = pointprocess::spectral::areConj(v1,v2,tol);
            EXPECT_FALSE(result);
        }
        /*
        TEST(TestSpectral,TestComputePoles){
            Eigen::VectorXd thetaP(6);
            double var = ...;
            double fSamp = ...;
            auto result = pointprocess::spectral::computePoles(thetaP, var, fSamp);
        }
        */

    }
}

#endif //POINTPROCESS_TESTSPECTRAL_H
