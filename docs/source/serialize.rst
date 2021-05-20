Serialization
==============

Once you have a :code:`PointProcessResult` (output of :code:`RegressionPipeline::fullRegression()` you can use the following functions to serialize
both the :code:`regression parameters` at each time step (**ALERT**: The :code:`.csv` file generated tends to be quite large, e.g. :code:`1 GB` for :code:`3` hours of data)
and the :code:`Taus` (integral of hazard-rate for each observed event) in order to assess goodness of fit.

.. code-block:: cpp

    // ppRes is a PointProcessResult obtained by running RegressionPipeline::fullRegression()
    ppResData2csv(ppRes, std::string("myData.csv"));
    ppResTaus2csv(ppRes, std::string("myTaus.csv"));

The header of these files will be:

.. csv-table:: myData.csv
   :file: _static/myDataExample.csv
   :widths: 8,8,8,8,14,8,14,8,8,8,8
   :header-rows: 1

.. csv-table:: myTaus.csv
   :file: _static/myTausExample.csv
   :widths: 100
   :header-rows: 1




