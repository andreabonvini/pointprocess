Spectral Analysis
==================

The spectral analysis module computes power spectral density (PSD) and heart rate variability
(HRV) indices from the autoregressive (AR) model parameters estimated by the point process
regression pipeline.

All spectral computations are implemented in C++ for maximum efficiency and exposed to Python
via pybind11 bindings.

Functions
---------

``compute_spectral_analysis(thetap, mean_interval, variance, aggregate=True)``
    Compute the power spectral density from AR model coefficients.

    **Parameters:**

    - ``thetap`` — AR coefficients (numpy array)
    - ``mean_interval`` — mean inter-event interval in seconds
    - ``variance`` — model variance (sigma²) in seconds²
    - ``aggregate`` — if True (default), aggregate complex conjugate pole components

    **Returns:** a ``SpectralAnalysis`` object with attributes:

    - ``frequencies`` — frequency axis (Hz)
    - ``powers`` — power spectral density (ms²/Hz)
    - ``poles`` — list of ``Pole`` objects (position, frequency, power, residual)
    - ``comps`` — spectral components per pole (or aggregated conjugate pair)

    **Example:**

    .. code-block:: python

        from pointprocess import compute_spectral_analysis

        analysis = compute_spectral_analysis(
            thetap=result.thetap,
            mean_interval=result.mean_interval,
            variance=result.sigma ** 2
        )

``compute_auto_correlation(taus, maxlag=60)``
    Compute autocorrelation of the quantile-normal-transformed rescaled times for
    model goodness-of-fit assessment. Under a correct model, the rescaled times
    should be independent (uncorrelated at all lags).

    **Parameters:**

    - ``taus`` — rescaled inter-event intervals (list or array of floats)
    - ``maxlag`` — maximum lag to compute (default: 60)

    **Returns:** numpy array of autocorrelation coefficients for lags 1 through ``maxlag``.

HRV Indices
-----------

HRV frequency band power is computed by aggregating pole powers into standard bands:

- **VLF** (Very Low Frequency): ≤ 0.04 Hz
- **LF** (Low Frequency): 0.04–0.15 Hz
- **HF** (High Frequency): 0.15–0.45 Hz

These bands follow the definitions in `Barbieri et al. (2005)
<https://pubmed.ncbi.nlm.nih.gov/15374824/>`_ and the MATLAB reference implementation.

To compute HRV indices after a full regression:

.. code-block:: python

    result = compute_full_regression(...)
    result.compute_hrv_indices()
    d = result.to_dict()
    # d["powVLF"], d["powLF"], d["powHF"] are now available

Implementation Notes
--------------------

- The PSD is computed on a 2048-point frequency grid.
- A pole stability fix (Stoica and Moses, *Signal Processing* 26(1), 1992) is applied
  to prevent slightly unstable AR models from producing invalid spectral estimates.
- Hamming-window low-pass filtering is applied to the time-varying HRV power traces
  to reduce estimation noise, matching the MATLAB reference.