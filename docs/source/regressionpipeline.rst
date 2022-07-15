RegressionPipeline
====================

A :code:`RegressionPipeline` object is all you need in order to perform a :code:`fullRegression` on a set of events.

.. code-block:: cpp

    #include "../../pointprocess/InterEventDistributions.h"
    #include "../../pointprocess/RegressionPipeline.h"
    #include "../../pointprocess/WeightsProducer.h"
    #include "../../tests/testData.h"

    #include <vector>
    #include "../../pointprocess/serialize.h"



    int main() {
        auto td = getTestData();
        auto pip = RegressionPipeline(
                Distributions::InverseGaussian, // distribution
                9, // AR_ORDER
                true // hasTheta0
                );
        auto ppRes = pip.fullRegression(
                td.testEvents, // event times
                60.0,  // windowLength
                0.005, // delta
                true,  // rightCensoring
                1000,  // maxIter
                WeightsProducer(0.98)
                );
        // Serialization...
        ppResData2csv(ppRes, std::string("myData.csv"));
        ppResTaus2csv(ppRes, std::string("myTaus.csv"));
    }

The output of the :code:`fullRegression()` method is a :code:`PointProcessResult`:

.. code-block:: cpp

    struct PointProcessResult{
        std::vector<std::shared_ptr<RegressionResult>> results;
        std::vector<double> taus;
        Distributions distribution;
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