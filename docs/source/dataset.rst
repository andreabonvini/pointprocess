PointProcessDataset
====================

A :code:`PointProcessDataset` is needed in order to compute a :code:`singleRegression` (see :ref:`Optimizers`).

The easiest way to construct a :code:`PointProcessDataset` is by using the  static method :code:`load()`, you just have to supply:

- :code:`std::deque<double> events_times`
- :code:`unsigned char AR_ORDER`
- :code:`bool hasTheta0`
- :code:`WeightsProducer& weightsProducer`
- :code:`double current_time`

:code:`current_time` is the time at which we're evaluating our model, should be greater or equal than the last element of :code:`events_times` (which is the default value). It is used only in case :code:`rightCensoring` is applied.

Example:

.. code-block:: cpp

    auto dataset = PointProcessDataset::load(
            events_times,
            AR_ORDER,
            hasTheta0,
            weightsProducer
        )

