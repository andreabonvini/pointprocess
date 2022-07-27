//
// Created by Andrea Bonvini on 18/05/21.
//

#include "pointprocess/PointProcessUtils.h"
#include "pointprocess/RegressionPipeline.h"
#include "pointprocess/spectral/spectral.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // Some automatic conversions are optional and require extra
// headers to be included when compiling your pybind11 module.
#include <pybind11/eigen.h> // Needed for Eigen conversions.
#include <pybind11/iostream.h> // Needed to redirect stdout to Python

namespace py = pybind11;


pointprocess::Result computeFullRegression(
        std::vector<double>& events,
        double windowLength,
        double delta,
        unsigned char AR_ORDER,
        bool hasTheta0,
        bool rightCensoring,
        double alpha,
        pointprocess::Distributions distribution,
        unsigned int maxIter
){

    auto pip = RegressionPipeline(distribution, AR_ORDER, hasTheta0);
    return pip.fullRegression(
            events,
            windowLength,
            delta,
            rightCensoring,
            maxIter,
            alpha
            );
}


std::shared_ptr<pointprocess::RegressionResult> computeSingleRegression(
        std::deque<double>& events,
        unsigned char AR_ORDER,
        bool hasTheta0,
        bool rightCensoring,
        double alpha,
        pointprocess::Distributions distribution,
        unsigned int maxIter
){

    auto wp = WeightsProducer(alpha);

    auto dataset = PointProcessDataset::load(
            events, // event_times
            AR_ORDER, // AR_ORDER
            hasTheta0, // hasTheta0
            wp
    );

    auto factory = OptimizersFactory();
    auto optimizer = factory.create(distribution);

    return optimizer->singleRegression(
            dataset,
            rightCensoring,
            maxIter
    );
}


// Create Python Bindings


PYBIND11_MODULE(pointprocess, m) {
    m.doc() = "PointProcess module - https://github.com/andreabonvini/pointprocess";

    py::enum_<pointprocess::Distributions>(m, "Distributions")
            .value("InverseGaussian", pointprocess::Distributions::InverseGaussian)
            .value("LogNormal", pointprocess::Distributions::LogNormal)
            .value("Gaussian", pointprocess::Distributions::Gaussian);

    py::class_<pointprocess::spectral::Pole>(m, "Pole")
            .def(py::init<std::complex<double> , double, double, std::complex<double>>(),
                    py::arg("pos"),
                    py::arg("frequency"),
                    py::arg("power"),
                    py::arg("residual")
                    )
            .def_readwrite("pos", &pointprocess::spectral::Pole::pos)
            .def_readwrite("frequency", &pointprocess::spectral::Pole::frequency)
            .def_readwrite("power", &pointprocess::spectral::Pole::power)
            .def_readwrite("residual", &pointprocess::spectral::Pole::residual);


    py::class_<pointprocess::spectral::SpectralAnalysis>(m, "SpectralAnalysis")
            .def(py::init<Eigen::VectorXd , Eigen::VectorXd, std::vector<pointprocess::spectral::Pole>&, std::vector<Eigen::VectorXcd>>(),
                 py::arg("frequencies"),
                 py::arg("powers"),
                 py::arg("poles"),
                 py::arg("comps")
            )
            .def_readwrite("frequencies", &pointprocess::spectral::SpectralAnalysis::frequencies)
            .def_readwrite("powers", &pointprocess::spectral::SpectralAnalysis::powers)
            .def_readwrite("poles", &pointprocess::spectral::SpectralAnalysis::poles)
            .def_readwrite("comps", &pointprocess::spectral::SpectralAnalysis::comps);


    py::class_<pointprocess::spectral::HeartRateVariabilityIndices>(m, "HeartRateVariabilityIndices")
            .def(py::init<double,double,double>(),
                 py::arg("powVLF"),
                 py::arg("powLF"),
                 py::arg("powHF")
            )
            .def_readwrite("powVLF", &pointprocess::spectral::HeartRateVariabilityIndices::powVLF)
            .def_readwrite("powLF", &pointprocess::spectral::HeartRateVariabilityIndices::powLF)
            .def_readwrite("powHF", &pointprocess::spectral::HeartRateVariabilityIndices::powHF);



    py::class_<pointprocess::RegressionResult, std::shared_ptr<pointprocess::RegressionResult>>(m, "RegressionResult")
            .def(
                py::init<double, Eigen::VectorXd, double, double, double, double, unsigned long, double, double, bool, bool, bool, double>(),
                py::arg("theta0"),
                py::arg("thetap"),
                py::arg("mu"),
                py::arg("sigma"),
                py::arg("lambda_"),
                py::arg("mean_interval"),
                py::arg("n_iter"),
                py::arg("likelihood"),
                py::arg("max_grad"),
                py::arg("converged"),
                py::arg("cdf_is_one"),
                py::arg("event_happened"),
                py::arg("time")
             )
            .def_readwrite("theta0", &pointprocess::RegressionResult::theta0)
            .def_readwrite("thetap", &pointprocess::RegressionResult::thetaP)
            .def_readwrite("mu", &pointprocess::RegressionResult::mu)
            .def_readwrite("sigma", &pointprocess::RegressionResult::sigma)
            .def_readwrite("lambda_", &pointprocess::RegressionResult::lambda)
            .def_readwrite("mean_interval", &pointprocess::RegressionResult::meanInterval)
            .def_readwrite("n_iter", &pointprocess::RegressionResult::nIter)
            .def_readwrite("likelihood", &pointprocess::RegressionResult::likelihood)
            .def_readwrite("max_grad", &pointprocess::RegressionResult::maxGrad)
            .def_readwrite("converged", &pointprocess::RegressionResult::converged)
            .def_readwrite("cdf_is_one", &pointprocess::RegressionResult::cdfIsOne)
            .def_readwrite("time", &pointprocess::RegressionResult::time)
            .def_readwrite("hrv_indices", &pointprocess::RegressionResult::hrvIndices)
            .def("compute_hrv_indices",
                 [](pointprocess::RegressionResult &self) {
                     py::scoped_ostream_redirect stream(
                             std::cout,                               // std::ostream&
                             py::module_::import("sys").attr("stdout") // Python output
                     );
                     self.computeHRVIndices();
                 }
            );

    py::class_<pointprocess::IGRegressionResult, std::shared_ptr<pointprocess::IGRegressionResult>, pointprocess::RegressionResult>(m, "IGRegressionResult")
            .def(
                    py::init<double, Eigen::VectorXd, double, double, double, double, unsigned long, double, double, double, bool, bool, bool, double>(),
                    py::arg("theta0"),
                    py::arg("thetap"),
                    py::arg("mu"),
                    py::arg("sigma"),
                    py::arg("lambda_"),
                    py::arg("mean_interval"),
                    py::arg("n_iter"),
                    py::arg("likelihood"),
                    py::arg("max_grad"),
                    py::arg("kappa"),
                    py::arg("converged"),
                    py::arg("cdf_is_one"),
                    py::arg("event_happened"),
                    py::arg("time")
            )
            .def_readwrite("theta0", &pointprocess::IGRegressionResult::theta0)
            .def_readwrite("thetap", &pointprocess::IGRegressionResult::thetaP)
            .def_readwrite("mu", &pointprocess::IGRegressionResult::mu)
            .def_readwrite("sigma", &pointprocess::IGRegressionResult::sigma)
            .def_readwrite("lambda_", &pointprocess::IGRegressionResult::lambda)
            .def_readwrite("mean_interval", &pointprocess::IGRegressionResult::meanInterval)
            .def_readwrite("n_iter", &pointprocess::IGRegressionResult::nIter)
            .def_readwrite("likelihood", &pointprocess::IGRegressionResult::likelihood)
            .def_readwrite("max_grad", &pointprocess::IGRegressionResult::maxGrad)
            .def_readwrite("kappa", &pointprocess::IGRegressionResult::kappa)
            .def_readwrite("converged", &pointprocess::IGRegressionResult::converged)
            .def_readwrite("cdf_is_one", &pointprocess::IGRegressionResult::cdfIsOne)
            .def_readwrite("event_happened", &pointprocess::IGRegressionResult::eventHappened)
            .def_readwrite("time", &pointprocess::IGRegressionResult::time);


    py::class_<pointprocess::Stats>(m, "Stats")
            .def(
                    py::init<double, double, double>()
            )
            .def_readwrite("ks_distance", &pointprocess::Stats::ksDistance)
            .def_readwrite("perc_out", &pointprocess::Stats::percOut)
            .def_readwrite("auto_corr", &pointprocess::Stats::autoCorr);

    py::class_<pointprocess::KsCoords>(m,"KsCoords")
            .def(
                    py::init<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd,Eigen::VectorXd>(),
                            py::arg("z"),
                            py::arg("inner"),
                            py::arg("upper"),
                            py::arg("lower")
            )
            .def_readwrite("z", &pointprocess::KsCoords::z)
            .def_readwrite("inner", &pointprocess::KsCoords::lin)
            .def_readwrite("upper", &pointprocess::KsCoords::lu)
            .def_readwrite("lower", &pointprocess::KsCoords::ll);



    py::class_<pointprocess::Result>(m, "Result")
            .def(
                    py::init<
                            std::vector<std::shared_ptr<pointprocess::RegressionResult>>,
                            std::vector<double>,
                            pointprocess::Distributions,
                            unsigned char,
                            bool,
                            double,
                            double,
                            double,
                            pointprocess::Stats
                            >(),
                            py::arg("results"),
                            py::arg("taus"),
                            py::arg("distribution"),
                            py::arg("ar_order"),
                            py::arg("has_theta0"),
                            py::arg("window_length"),
                            py::arg("delta"),
                            py::arg("t0"),
                            py::arg("stats")

            )
            .def("compute_hrv_indices",
                 [](pointprocess::Result &self) {
                     py::scoped_ostream_redirect stream(
                             std::cout,                               // std::ostream&
                             py::module_::import("sys").attr("stdout") // Python output
                     );
                     self.computeHRVIndices();
                 }
                 )
             .def("to_dict",
                  [](pointprocess::Result &self) {
                      py::scoped_ostream_redirect stream(
                              std::cout,                               // std::ostream&
                              py::module_::import("sys").attr("stdout") // Python output
                      );
                      return self.toDict();
                  }
                 )
            .def_readwrite("results", &pointprocess::Result::results)
            .def_readwrite("taus", &pointprocess::Result::taus)
            .def_readwrite("distribution", &pointprocess::Result::distribution)
            .def_readwrite("ar_order", &pointprocess::Result::AR_ORDER)
            .def_readwrite("has_theta0", &pointprocess::Result::hasTheta0)
            .def_readwrite("windows_length", &pointprocess::Result::windowLength)
            .def_readwrite("delta", &pointprocess::Result::delta)
            .def_readwrite("t0", &pointprocess::Result::t0)
            .def_readwrite("stats", &pointprocess::Result::stats);

    m.def(
            "compute_single_regression",
            &computeSingleRegression,
            "This function performs a single point process regression on a series of events.",
            py::arg("events"),
            py::arg("ar_order"),
            py::arg("has_theta0"),
            py::arg("right_censoring"),
            py::arg("alpha"),
            py::arg("distribution"),
            py::arg("max_iter")
    );

    m.def(
            "get_ks_coords",
            &pointprocess::utils::getKsCoords,
            "This function returns the transformed taus and coordinates useful to create a KS plot.",
            py::arg("taus")
            );

    m.def(
            "compute_full_regression",
            []( std::vector<double>& events,
                double windowLength,
                double delta,
                unsigned char AR_ORDER,
                bool hasTheta0,
                bool rightCensoring,
                double alpha,
                pointprocess::Distributions distribution,
                unsigned int maxIter
                ) {
                py::scoped_ostream_redirect stream(
                        std::cout,                               // std::ostream&
                        py::module_::import("sys").attr("stdout") // Python output
                );
                return computeFullRegression(events, windowLength, delta, AR_ORDER, hasTheta0, rightCensoring, alpha,distribution,maxIter);
            },
            "This function performs a full point process regression on a series of events.",
            py::arg("events"),
            py::arg("window_length"),
            py::arg("delta"),
            py::arg("ar_order"),
            py::arg("has_theta0"),
            py::arg("right_censoring"),
            py::arg("alpha"),
            py::arg("distribution"),
            py::arg("max_iter")
    );

    m.def(
            "compute_spectral_analysis",
            []( Eigen::VectorXd& thetaP, //TODO: DOES THIS COPY thetap each time!?
                double meanInterval,
                double variance,
                bool aggregate
            ) {
                py::scoped_ostream_redirect stream(
                        std::cout,                               // std::ostream&
                        py::module_::import("sys").attr("stdout") // Python output
                );
                return pointprocess::spectral::computeSpectralAnalysis(
                        thetaP,meanInterval,variance,aggregate
                        );
            },
            "This function performs a full point process regression on a series of events.",
            py::arg("thetap"),
            py::arg("mean_interval"),
            py::arg("variance"),
            py::arg("aggregate")
    );
}


