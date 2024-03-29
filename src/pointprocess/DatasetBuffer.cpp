//
// Created by Andrea Bonvini on 28/05/22.
//

#include "DatasetBuffer.h"


DatasetBuffer::Iterator::Iterator(
        DatasetBuffer* ptr,
        unsigned long bin_index
        ) :
        m_ptr(ptr) ,
        bin_index(bin_index),
        values(
                std::move(
                        std::make_tuple(
                                m_ptr->currentTime_ ,
                                m_ptr->eventHappened_,
                                m_ptr->resetParameters_,
                                m_ptr->current_
                        )
                )
        ) {}

// TODO: what is this const?
DatasetBuffer::Iterator::value_type DatasetBuffer::Iterator::operator*() const {return values;}
// reference operator*() const { return *m_ptr; }

// pointer operator->() { return m_ptr; }
DatasetBuffer::Iterator::pointer DatasetBuffer::Iterator::operator->() { return &values; }

// Prefix increment
DatasetBuffer::Iterator& DatasetBuffer::Iterator::operator++() {
    bin_index++;
    m_ptr->update(bin_index);
    std::get<0>(values) = m_ptr->currentTime_;
    std::get<1>(values) = m_ptr->eventHappened_;
    std::get<2>(values) = m_ptr->resetParameters_;
    // FIXME: is this a copy?
    std::get<3>(values) = m_ptr->current_;

    return *this;
}



// Postfix increment // FIXME: is this right?
// Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }

bool operator==(const DatasetBuffer::Iterator& a, const DatasetBuffer::Iterator& b) { return a.bin_index == b.bin_index;}
bool operator!=(const DatasetBuffer::Iterator& a, const DatasetBuffer::Iterator& b) { return a.bin_index != b.bin_index;}

    // ============================ End Iterator definition ==================================

DatasetBuffer::DatasetBuffer(pointprocess::PipelineSetup& setup)
            :
            AR_ORDER_(setup.AR_ORDER),
            hasTheta0_(setup.hasTheta0),
            weightsProducer_(setup.weightsProducer),
            delta_(setup.delta),
            events_(setup.events), // DO NOT std::move!
            observed_events_(std::deque<double>(events_.begin(), events_.begin() + (long) setup.last_event_index + 1)),
            windowLength_(setup.delta * (double) (setup.bins_in_window)),
            bins_in_window_(setup.bins_in_window),
            last_bin_index_(setup.bins),
            last_event_index_(setup.last_event_index),
            currentTime_((double) setup.bins_in_window * setup.delta),
            current_(
                    PointProcessDataset::load(
                            observed_events_,
                            setup.AR_ORDER,
                            setup.hasTheta0,
                            setup.weightsProducer,
                            (double) setup.bins_in_window * setup.delta)
            ) {
    }

unsigned int DatasetBuffer::getNumberOfRegressionParameters() const{
        return AR_ORDER_ + (unsigned int) hasTheta0_;
    }

DatasetBuffer::Iterator DatasetBuffer::begin() {
        return DatasetBuffer::Iterator(this, bins_in_window_); }

DatasetBuffer::Iterator DatasetBuffer::end()   {
        return DatasetBuffer::Iterator(this, last_bin_index_ + 1); }

unsigned long DatasetBuffer::size() const{
    return last_bin_index_ - bins_in_window_ + 1;
}

void DatasetBuffer::update(unsigned long bin_index){
    if (bin_index <= last_bin_index_) {
        currentTime_ = (double) bin_index * delta_;

        // If in the previous time step we reset our distribution parameters (since the observed events changed), we
        // can set this flag to false.
        if (resetParameters_) {
            resetParameters_ = false;
        }

        /* If the first element of observed_events happened before the
         * time window between (current_time - window_length) and (current_time)
         * we can discard it since it will not be part of the current optimization process.
         */
        if (!observed_events_.empty() && observed_events_[0] < currentTime_ - windowLength_) {
            observed_events_.pop_front();
            // Force re-evaluation of starting point for theta and kappa.
            resetParameters_ = true;
        }
        // We check whether an event happened in ((bin_index - 1) * delta, bin_index * delta]
        eventHappened_ = events_[last_event_index_ + 1] <= currentTime_;

        if (eventHappened_) {
            last_event_index_++;
            observed_events_.push_back(events_[last_event_index_]);
            resetParameters_ = true;
        }

        // We create a PointProcessDataset for the current time bin
        current_ = PointProcessDataset::load(observed_events_, AR_ORDER_, hasTheta0_, weightsProducer_,
                                             currentTime_);
    }
}