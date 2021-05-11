.. image:: ../images/ppbig.png


A C++ Library for Point Process analysis.
==========================================

.. toctree::
   :maxdepth: 2
   :caption: Basic info:

   license
   help

   
.. toctree::
   :maxdepth: 2
   :caption: Scientific docs:

   Gradients-and-Hessians
   inhomogenouspoisson
   hazardfunction


Quick Tour
==================

.. code-block:: c

    #include "pointprocess/InterEventDistributions.h"
    #include "pointprocess/RegressionPipeline.h"
    #include "pointprocess/WeightsProducer.h"
    #include "tests/testData.h"
    #include <vector>

    int main(){
        std::vector<double> td = getTestData();
        auto pip = RegressionPipeline(PointProcessDistributions::InverseGaussian, 9, true);
        auto fullRes = pip.fullRegression(td.testEvents, 60.0, 0.005, true, 1000, WeightsProducer(0.98));
        // Serialize...
	}
