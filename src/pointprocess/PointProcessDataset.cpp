//
// Created by Andrea Bonvini on 20/04/21.
//

#include "PointProcessDataset.h"

PointProcessDataset::PointProcessDataset(
        unsigned char N_SAMPLES_,
        unsigned char AR_ORDER_,
        bool hasTheta0_,
        Eigen::MatrixXd xn_,
        Eigen::VectorXd wn_,
        Eigen::VectorXd eta_,
        Eigen::VectorXd xt_,
        double wt_
){
    N_SAMPLES = N_SAMPLES_;
    AR_ORDER = AR_ORDER_;
    hasTheta0 = hasTheta0_;
    xn = std::move(xn_);
    wn = std::move(wn_);
    eta = std::move(eta_);
    xt = std::move(xt_);
    wt = wt_;
}

PointProcessDataset PointProcessDataset::load(
        std::deque<double> events_times,
        unsigned char AR_ORDER_,
        bool hasTheta0_,
        WeightsProducer& weightsProducer,
        double current_time
){
    // In order to build a dataset we need at least 1 sample (i.e. AR_ORDER_ + 2 events)
    // e.g if we have an AR_ORDER = 3 and 5 events, we then have 4 inter arrival times from which we can
    // build 1 sample.
    assert (events_times.size() >= AR_ORDER_ + 2);
    std::deque<double> inter_events_times;
    inter_events_times = events_times;
    std::adjacent_difference(inter_events_times.begin(),inter_events_times.end(),inter_events_times.begin());
    inter_events_times.pop_front();
    // We now define wn, i.e. the intervals we have to predict once we build our AR model.
    std::vector<double> wn_v = std::vector<double>(inter_events_times.begin() + AR_ORDER_, inter_events_times.end());
    // Let's copy it into an Eigen Vector object.
    Eigen::VectorXd wn_ = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(wn_v.data(), (long) wn_v.size());
    // We now have to build a matrix xn s.t. for i = 0, ..., len(inter_events_times)-p-1 the i_th element of xn will be
    // xn[i] = [1, inter_events_times[i + p - 1], inter_events_times[i + p - 2], ..., rr[i]]

    /*
     * Here we have to distinguish between the case in which (AR_ORDER == 0 && hasTheta == true) and the other cases.
     * Note that when AR_ORDER = 0 and hasTheta0 = true we are directly optimizing for the mean of the distribution
     * (since we would have mu = theta0) without trying to learn any autoregressive coefficient (which are mainly
     * used for parametric spectral estimation when dealing with HRV data).
     */
    Eigen::MatrixXd xn_(events_times.size() - 1 - AR_ORDER_ , hasTheta0_ + AR_ORDER_);
    assert(!(AR_ORDER_ == 0 && !hasTheta0_)); // This is the only forbidden case.
    if (AR_ORDER_ == 0 && hasTheta0_){
        xn_ << Eigen::MatrixXd::Ones( (long) wn_.size(),1);
    }
    else {
        // a = inter_events_times[p - 1 : -1]
        std::vector<double> a = std::vector<double>(inter_events_times.begin() + AR_ORDER_ - 1,inter_events_times.end() - 1);
        // b = inter_events_times[p - 1 :: -1]
        std::vector<double> b = std::vector<double>(inter_events_times.begin(),inter_events_times.begin() + AR_ORDER_);
        std::reverse(b.begin(), b.end());
        // xn = toeplitz(a, b)
        Eigen::MatrixXd xn_tmp(a.size(), b.size());
        toeplitz(a, b, xn_tmp);
        // Note that the 1 at the beginning of each row is added only if the hasTheta0 parameter is set to True.
        if (hasTheta0_){
            auto ones = Eigen::MatrixXd::Ones( (long) xn_tmp.rows(),1);
            xn_ << ones, xn_tmp;
        }
        else{
            xn_ = xn_tmp;
        }
    }

    // xt = inter_events_times[-p:][::-1]
    std::vector<double> xt_v = std::vector<double>(inter_events_times.end() - AR_ORDER_, inter_events_times.end());
    std::reverse(xt_v.begin(), xt_v.end());
    Eigen::VectorXd xt_tmp = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(xt_v.data(), (long) xt_v.size());
    Eigen::VectorXd xt_(AR_ORDER_ + hasTheta0_);
    if (hasTheta0_) {
        auto one = Eigen::VectorXd::Ones((long) 1);
        xt_ << one, xt_tmp;
    }
    else
        xt_ = xt_tmp;
    if (current_time == 0.0)
        current_time = events_times[events_times.size() - 1];

    assert(current_time >= events_times[events_times.size() - 1]);
    double wt_ = current_time - events_times[events_times.size() - 1];
    // target_distances = current_time - event_times[p + 1 :]
    std::vector<double> target_distances;
    for(auto it = events_times.begin() + AR_ORDER_ + 1; it != events_times.end(); ++it){
        target_distances.push_back( current_time - *it);
    }
    // eta = weights_producer(current_time - uk)
    Eigen::VectorXd eta_ = weightsProducer.produce(target_distances);
    unsigned char N_SAMPLES = xn_.cols();
    return PointProcessDataset(N_SAMPLES , AR_ORDER_, hasTheta0_, xn_, wn_, eta_, xt_, wt_);
}

void PointProcessDataset::toeplitz(const std::vector<double>& col, const std::vector<double>& row, Eigen::MatrixXd& toep){
    assert(col[0] == row[0]);
    for(long i = 0; i != col.size(); i++) {
        for(long j = 0; j != row.size(); j++) {
            if (j >= i){
                // Upper diagonal
                toep(i,j) = row[j - i];
            }
            else{
                // Lower diagonal
                toep(i,j) = col[i - j];
            }
        }
    }
}
