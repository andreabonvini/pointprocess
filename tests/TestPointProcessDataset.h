//
// Created by Andrea Bonvini on 21/04/21.
//

#ifndef POINTPROCESS_TESTPOINTPROCESSDATASET_H
#define POINTPROCESS_TESTPOINTPROCESSDATASET_H
#include "../pointprocess/PointProcessDataset.h"
#include <Eigen/Core>

bool testToeplitz1(){
    std::vector<double> col {1.0,2.0,3.0};
    std::vector<double> row {1.0,5.0,6.0};
    Eigen::MatrixXd expected(col.size(),row.size());
    expected << 1.0,5.0,6.0,
                2.0,1.0,5.0,
                3.0,2.0,1.0;

    Eigen::MatrixXd res(col.size(),row.size());
    toeplitz(col,row,res);
    return res.isApprox(expected);
}

bool testToeplitz2(){
    std::vector<double> col {1.0,2.0,3.0,4.0};
    std::vector<double> row {1.0,6.0,7.0};
    Eigen::MatrixXd expected(col.size(),row.size());
    expected << 1.0,6.0,7.0,
                2.0,1.0,6.0,
                3.0,2.0,1.0,
                4.0,3.0,2.0;
    Eigen::MatrixXd res(col.size(),row.size());
    toeplitz(col,row,res);

    return res.isApprox(expected);
}

bool testToeplitz3(){
    std::vector<double> col {5.0,6.0,7.0};
    std::vector<double> row {5.0,2.0,3.0,4.0};
    Eigen::MatrixXd expected(col.size(),row.size());
    expected << 5.0,2.0,3.0,4.0,
                6.0,5.0,2.0,3.0,
                7.0,6.0,5.0,2.0;
    Eigen::MatrixXd res(col.size(),row.size());
    toeplitz(col,row,res);
    return res.isApprox(expected);
}


bool testPointProcessDataset(){
    // TODO: split in n different tests...
    std::deque<double> events = {
            727.391,
            728.297,
            729.188,
            730.062,
            730.984,
            731.93,
            732.875,
            733.828,
            734.781,
            735.711,
    };
    unsigned char ar_order = 3;
    unsigned char n_samples = 6;
    bool hasTheta0 = false;
    Eigen::VectorXd expected_wn(n_samples);
    Eigen::VectorXd expected_eta(n_samples);
    Eigen::VectorXd expected_xt(ar_order + hasTheta0);
    Eigen::MatrixXd expected_xn(n_samples, ar_order + hasTheta0);

    // note that diff: 0.906, 0.891, 0.874, 0.922, 0.946, 0.945, 0.953, 0.953, 0.93

    expected_eta <<  1.000, 1.000, 1.000, 1.000, 1.000, 1.000;

    expected_wn << 0.922, 0.946, 0.945, 0.953, 0.953, 0.930;

    expected_xt << 0.930, 0.953, 0.953;

    expected_xn << 0.874, 0.891, 0.906,
                   0.922, 0.874, 0.891,
                   0.946, 0.922, 0.874,
                   0.945, 0.946, 0.922,
                   0.953, 0.945, 0.946,
                   0.953, 0.953, 0.945;

    auto wp = WeightsProducer(0.0);
    double current_time = 736.0;
    double expected_wt = current_time - events[events.size() - 1];

    auto dataset = PointProcessDataset::load(events,ar_order,hasTheta0,wp,current_time);

    return dataset.wn.isApprox(expected_wn) && dataset.xn.isApprox(expected_xn) && dataset.xt.isApprox(expected_xt) && dataset.eta.isApprox(expected_eta) && dataset.wt == expected_wt;
}

bool testPointProcessDataset_hasTheta0(){
    // TODO: split in n different tests...
    std::deque<double> events = {
            727.391,
            728.297,
            729.188,
            730.062,
            730.984,
            731.93,
            732.875,
            733.828,
            734.781,
            735.711,
    };
    unsigned char ar_order = 3;
    unsigned char n_samples = 6;
    bool hasTheta0 = true;
    Eigen::VectorXd expected_wn(n_samples);
    Eigen::VectorXd expected_eta(n_samples);
    Eigen::VectorXd expected_xt(ar_order + hasTheta0);
    Eigen::MatrixXd expected_xn(n_samples, ar_order + hasTheta0);

    // note that diff: 0.906, 0.891, 0.874, 0.922, 0.946, 0.945, 0.953, 0.953, 0.93

    expected_eta <<  1.000, 1.000, 1.000, 1.000, 1.000, 1.000;

    expected_wn << 0.922, 0.946, 0.945, 0.953, 0.953, 0.930;

    expected_xt << 1.000, 0.930, 0.953, 0.953;

    expected_xn << 1.000, 0.874, 0.891, 0.906,
                   1.000, 0.922, 0.874, 0.891,
                   1.000, 0.946, 0.922, 0.874,
                   1.000, 0.945, 0.946, 0.922,
                   1.000, 0.953, 0.945, 0.946,
                   1.000, 0.953, 0.953, 0.945;

    auto wp = WeightsProducer(0.0);
    double current_time = 736.0;
    double expected_wt = current_time - events[events.size() - 1];

    auto dataset = PointProcessDataset::load(events,ar_order,hasTheta0,wp,current_time);

    return dataset.wn.isApprox(expected_wn) && dataset.xn.isApprox(expected_xn) && dataset.xt.isApprox(expected_xt) && dataset.eta.isApprox(expected_eta) && dataset.wt == expected_wt;
}

#endif //POINTPROCESS_TESTPOINTPROCESSDATASET_H
