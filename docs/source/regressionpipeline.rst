RegressionPipeline
====================

A :code:`RegressionPipeline` object is all you need in order to perform a :code:`fullRegression` on a set of events.

Based on the specific distribution you want to use to fit your data you can instantiate one of three :code:`BaseOptimizer` objects:

- :code:`auto opt = InverseGaussianOptimizer()`
- :code:`auto opt = GaussianOptimizer()`
- :code:`auto opt = LogNormalOptimizer()`

Once you have instatiated your optimizer, you have to build a :code:`WeightsProducer` instance, a :code:`PointProcessDataset` and call the :code:`.singleRegression()` method as you
can see in the example below:

.. code-block:: cpp

    #include "pointprocess/InterEventDistributions.h"
    #include "pointprocess/WeightsProducer.h"
    #include "pointprocess/PointProcessDataset.h"
    #include "pointprocess/RegressionPipeline.h"


    int main()
    {
        // Retrieve some data...
        auto td = getTestData();
        auto wp = WeightsProducer(0.98);

        // Initialize Pipelined...
        auto pip = RegressionPipeline(
            PointProcessDistributions::InverseGaussian,
            9,   // AR_ORDER
            true // hasTheta0
        );

        // Compute regression...
        auto fullRes = pip.fullRegression(
            td.testEvents, // std::vector<double> (events times)
            60.0,          // windowLength
            0.005,         // delta
            true,          // rightCensoring
            1000,          // maxIter
            wp
        );
        // Serialize...
        ppRes2csv(fullRes, std::string("../notebooks/myData.csv"), std::string("../notebooks/myTaus.csv"));
    }

The output of the :code:`fullRegression()` method is a :code:`PointProcessResult`:

.. code-block:: cpp

    struct PointProcessResult{
        std::vector<std::shared_ptr<RegressionResult>> results;
        std::vector<double> taus;
        PointProcessDistributions distribution;
        double percOut;
        double ksDistance;
        double t0;
        double autoCorr;
        unsigned char AR_ORDER;
        bool hasTheta0;
        double windowLength;
        double delta;
    // ...
    }