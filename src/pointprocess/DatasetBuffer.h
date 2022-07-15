//
// Created by Andrea Bonvini on 31/12/21.
//

#ifndef POINTPROCESS_DATASETBUFFER_H
#define POINTPROCESS_DATASETBUFFER_H

#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t
#include "PointProcessDataset.h"
#include "PointProcessUtils.h"


class DatasetBuffer {
public:
    // ============================ Start Iterator definition ==================================
    struct Iterator
    {
        using value_type = std::tuple<double, bool, bool, PointProcessDataset>;
        using pointer  = value_type*;
        using reference = value_type&;

        Iterator(DatasetBuffer* ptr, unsigned long bin_index);

        // TODO: what is this const?
        value_type operator*() const;
        // reference operator*() const { return *m_ptr; }

        // pointer operator->() { return m_ptr; }
        pointer operator->();

        // Prefix increment
        Iterator& operator++();

        // Postfix increment // FIXME: is this right?
        // Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }

        friend bool operator== (const Iterator& a, const Iterator& b);
        friend bool operator!= (const Iterator& a, const Iterator& b);
    private:
        DatasetBuffer* m_ptr;
        value_type values;
        unsigned long bin_index;
    };
    // ============================ End Iterator definition ==================================

    explicit DatasetBuffer(pointprocess::PipelineSetup& setup);

    [[nodiscard]] unsigned int getNumberOfRegressionParameters() const;

    Iterator begin();
    Iterator end();

    [[nodiscard]] unsigned long size() const;

private:
    unsigned char AR_ORDER_;
    bool hasTheta0_;
    WeightsProducer weightsProducer_;
    double delta_;
    std::vector<double> events_;
    /* observed_events here is the subset of events observed during the first window, this std::deque will keep track
     * of the events used for local regression at each time bin, discarding old events and adding new ones.
     * It works as a buffer for our regression pipeline.
     */
    std::deque<double> observed_events_;
    double windowLength_;
    unsigned long bins_in_window_;
    unsigned long last_bin_index_;
    unsigned long last_event_index_;
    double currentTime_;
    bool eventHappened_ = false;
    bool resetParameters_ = true;
    PointProcessDataset current_;

    void update(unsigned long bin_index);
};

#endif //POINTPROCESS_DATASETBUFFER_H
