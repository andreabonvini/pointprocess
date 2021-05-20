Optimizers
===========

A :code:`BaseOptimizer` object is all you need in order to perform a :code:`singleRegression` on a set of events.

Based on the specific distribution you want to use to fit your data you can instantiate one of three :code:`BaseOptimizer` objects:

- :code:`auto opt = InverseGaussianOptimizer()`
- :code:`auto opt = GaussianOptimizer()`
- :code:`auto opt = LogNormalOptimizer()`

Once you have instantiated your optimizer, you have to build a :code:`WeightsProducer` instance, a :ref:`PointProcessDataset` and call the :code:`.singleRegression()` method as you
can see in the example below:

.. code-block:: cpp

    #include "pointprocess/WeightsProducer.h"
    #include "pointprocess/PointProcessDataset.h"
    #include "pointprocess/optimizers/LogNormalOptimizer.h"


    int main()
    {
        // Retrieve some data...
        auto td = getTestData();
        // Since td.testEvents is a std::vector<double>, in this case we need to copy it inside a std::deque<double>
        // in order to build a PointProcessDataset.
        std::deque<double> events;
        for (int i = 75; i < 300; i++){
            events.push_back(td.testEvents[i]);
        }

        auto wp = WeightsProducer(1.0);

        auto dataset = PointProcessDataset::load(
                events, // event_times
                9, // AR_ORDER
                true, // hasTheta0
                wp
                );

        auto optimizer = LogNormalOptimizer();
        auto result = optimizer.singleRegression(
                dataset,
                false, // rightCensoring
                10000 // maxIter
                );
    }

The output of the :code:`singleRegression()` method is a :code:`RegressionResult`:

.. code-block:: cpp

    struct RegressionResult{
        double theta0;
        Eigen::VectorXd thetaP;
        double mu;
        double sigma;
        double lambda;
        double meanInterval;
        double time = 0.0;
        unsigned long nIter;
        double likelihood;
        bool eventHappened = false;
    }