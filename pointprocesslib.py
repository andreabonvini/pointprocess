import ctypes
import numpy as np
from enum import Enum
import os
import matplotlib.pyplot as plt
from scipy.special import erfinv
from copy import deepcopy
from dataclasses import dataclass
from typing import List
from sys import platform

# TODO: Ideally all the functions in this file (except from regr_likel) should be implemented in C++ to obtain a
#  considerable gain in efficiency and to directly serialize all the needed information in a .csv file.
#  (especially the spectral functions are extremely slow in Python when compared to a pure C++ implementation.

# ================================================ START SPECTRAL ======================================

@dataclass
class Pole:
    """
    Pole or pair of conjugate poles,
    """

    # TODO add Documentation
    pos: complex  # position on the complex plane
    frequency: float  # pragma: no cover
    power: float  # pragma: no cover
    residual: float  # pragma: no cover


@dataclass
class SpectralAnalysis:
    frequencies: np.array  # Hz
    powers: np.array  # ms^2 / Hz
    poles: List[Pole]
    comps: List[List[complex]]


@dataclass
class HeartRateVariabilityIndices:
    powVLF: float
    powLF: float
    powHF: float


def hrv_indices(analysis: SpectralAnalysis) -> HeartRateVariabilityIndices:
    powVLF = sum(
        p.power for p in analysis.poles if np.abs(p.frequency) <= 0.04 and p.power > 0.0
    )
    powLF = sum(
        p.power
        for p in analysis.poles
        if 0.04 < np.abs(p.frequency) <= 0.15 and p.power > 0.0
    )
    powHF = sum(
        p.power
        for p in analysis.poles
        if 0.15 < np.abs(p.frequency) < 0.45 and p.power > 0.0
    )
    return HeartRateVariabilityIndices(powVLF, powLF, powHF)


def compute_psd(
        thetap: np.ndarray, mean_interval: float, variance: float, aggregate=True
) -> SpectralAnalysis:  # pragma: no cover
    var = 1e6 * variance  # from [s^2] to [ms^2]
    fsamp = 1 / mean_interval
    ar = np.r_[1, -thetap]  # [1, -θ1, -θ2, ..., θp]
    # Compute poles' complex values
    poles_values = np.roots(ar)
    # Order them by absolute angle
    poles_values = np.array(sorted(poles_values, key=lambda x: abs(np.angle(x))))

    # Fix AR models that might have become slightly unstable due to the estimation process
    # using an exponential decay (see Stoica and Moses, Signal Processing 26(1) 1992)
    mod_scale = min(0.99 / max(np.abs(poles_values)), 1)
    poles_values = list(mod_scale * poles_values)

    thetap = thetap * np.cumprod(np.ones(thetap.shape) * mod_scale)

    fs = np.linspace(-0.5, 0.5, 2048)
    # z = e^(-2πfT)
    # z: unit delay operator
    z = np.exp(2j * np.pi * fs)
    # P(z) = (σ^2*T)/ |1+θ1*z^(-1)+...+θp*z^(-p)|^2
    # σ^2 : Sample variance
    # T: sampling interval
    powers = (var / fsamp) / abs(np.polyval(np.r_[1, -thetap], np.conj(z))) ** 2
    frequencies = fs * fsamp
    poles_residuals = [
        1
        / (
                p
                * np.prod(p - np.array([val for val in poles_values if val is not p]))
                * np.prod(1 / p - np.conj(poles_values))
        )
        for p in poles_values
    ]

    poles_frequencies = [np.angle(p) / (2 * np.pi) * fsamp for p in poles_values]
    poles_powers = [var * np.real(p) for p in poles_residuals]
    # We also save the spectral components for each frequency value for each pole
    poles_comps = []
    ref_poles = 1 / np.conj(poles_values)
    for i in range(len(poles_values)):
        pp = poles_residuals[i] * poles_values[i] / (z - poles_values[i])
        refpp = -np.conj(poles_residuals[i]) * ref_poles[i] / (z - ref_poles[i])
        poles_comps.append(var / fsamp * (pp + refpp))

    poles_comps_agg = [deepcopy(poles_comps[0])]

    # Aggregate complex conjugate poles in poles_comps_agg

    for i in range(1, len(poles_values)):
        if np.isclose(poles_values[i], np.conj(poles_values[i - 1])):
            poles_comps_agg[-1] += poles_comps[i]
        else:
            poles_comps_agg.append(deepcopy(poles_comps[i]))

    poles = [
        Pole(pos, freq, power, res)
        for pos, freq, power, res in zip(
            poles_values, poles_frequencies, poles_powers, poles_residuals,
        )
    ]

    return SpectralAnalysis(
        frequencies,
        powers,
        poles,
        poles_comps_agg if aggregate else poles_comps,
    )

# ================================================ END SPECTRAL =======================================


def ks_distance(taus: np.ndarray, plot: bool = False):  # pragma: no cover
    """
    Compute KS-distance through the Time-Rescaling theorem
    """
    z = 1 - np.exp(-taus)
    z = sorted(z)
    d = len(z)
    lin = np.linspace(0, 1, d)
    if plot:
        plt.figure(figsize=(10, 10))
        lu = np.linspace(1.36 / np.sqrt(d), 1 + 1.36 / np.sqrt(d), d)
        ll = np.linspace(-1.36 / np.sqrt(d), 1 - 1.36 / np.sqrt(d), d)
        plt.plot(z, lin)
        plt.plot(lin, lin)
        plt.plot(lu, lin)
        plt.plot(ll, lin)
    KSdistance = max(abs(z - lin)) / np.sqrt(2.0)
    return KSdistance


def check_corr(taus: np.ndarray, maxlag: int = 60, plot: bool = True):
    Z = 1-np.exp(-taus)
    small = 0.00001
    Z = np.minimum(np.maximum(Z,small),1-small)
    N = erfinv(Z * 2 - 1)
    Nf = np.fft.fft(N - np.mean(N), 2 ** int(np.ceil(np.log2(len(N)))))
    Ns = np.abs(Nf)**2
    Ncorr = np.real(np.fft.ifft(Ns))
    Ncorr = Ncorr[1:maxlag+1] / Ncorr[0]
    c = Ncorr
    th = 2.0 / np.sqrt(len(Z))
    if plot:
        plt.plot(c, "k")
        plt.plot([0,maxlag],[th,th],"r--")
        plt.plot([0,maxlag],[-th,-th],"r--")
        plt.xlim([0,maxlag])
        plt.ylim([-0.6,0.6])


class InterEventDistribution(Enum):
    InverseGaussian = 0
    LogNormal = 1
    Gaussian = 2


def regrlikel(
        events: np.array,
        window_length: float,
        delta: float,
        ar_order: int,
        has_theta0: bool = True,
        right_censoring: bool = True,
        alpha: float = 0.02,
        distribution: InterEventDistribution = InterEventDistribution.InverseGaussian,
        max_iter: int = 1000,
        serialize_data: bool = True,
        output_data_path: str = "Data.csv",
        output_taus_path: str = "Taus.csv"
) -> None:
    assert len(events.shape) == 1
    n_events = len(events)
    c_events_pointer = events.astype(np.double).ctypes.data_as(c_double_p)
    # Compute result...
    cdll.regrlikel(
        n_events,
        c_events_pointer,
        window_length,
        delta,
        ar_order,
        has_theta0,
        right_censoring,
        alpha,
        distribution.value,
        max_iter,
        serialize_data,
        output_data_path.encode('utf-8'),
        output_taus_path.encode('utf-8')
    )


# COMMENT MARKED STUFF...
# def regrlikel_marked(
#         events: np.array,
#         amplitudes: np.array,
#         windowLength: float,
#         delta: float,
#         AR_ORDER: int,
#         hasTheta0: bool = True,
#         alpha: float = 0.02,
#         beta: float = 0.02,
#         distribution: InterEventDistribution = InterEventDistribution.Gaussian,
#         maxIter: int = 1000,
#         serializeData: bool = False,
#         outputDataName: str = "Data.csv",
#         outputTausName: str = "NOSENSE.csv"
# ) -> None:
#
#     assert len(events.shape) == 1
#     assert len(amplitudes.shape) == 1
#     assert len(events) == len(amplitudes)
#     n_events = len(events)
#     c_events_pointer = events.astype(np.double).ctypes.data_as(c_double_p)
#     c_amplitudes_pointer = amplitudes.astype(np.double).ctypes.data_as(c_double_p)
#     # Compute result...
#     cdll.regrlikel_marked(
#         n_events,
#         c_events_pointer,
#         c_amplitudes_pointer,
#         windowLength,
#         delta,
#         AR_ORDER,
#         hasTheta0,
#         alpha,
#         beta,
#         distribution.value,
#         maxIter,
#         serializeData,
#         outputDataName.encode('utf-8'),
#         outputTausName.encode('utf-8')
#     )

def _get_extension() -> str:

    # Checking whether the platform is MacOs, Windows or Linux.
    if platform == "darwin":
        return ".dylib"
    elif platform == "win32":
        return ".DLL"
    elif platform == "linux" or platform == "linux2":
        return ".so"
    else:
        raise Exception("_get_extension() function is broken.")


current_path = os.path.dirname(os.path.realpath(__file__))
LIB_NAME = "libpointprocess" + _get_extension()
lib_file = os.path.join(current_path, "build/src", LIB_NAME)
cdll = ctypes.cdll.LoadLibrary(lib_file)
c_double_p = ctypes.POINTER(ctypes.c_double)
cdll.regrlikel.argtypes = [
    ctypes.c_uint,    # n_events
    c_double_p,       # c_events_pointer
    ctypes.c_double,  # windowLength
    ctypes.c_double,  # delta
    ctypes.c_ubyte,   # AR_ORDER
    ctypes.c_bool,    # hasTheta0
    ctypes.c_bool,    # rightCensoring
    ctypes.c_double,  # alpha
    ctypes.c_uint,    # distribution
    ctypes.c_uint,    # maxIter
    ctypes.c_bool,    # serializeData
    ctypes.c_char_p,  # outputDataName
    ctypes.c_char_p   # outputTausName
 ]

# cdll.regrlikel_marked.argtypes = [
#     ctypes.c_uint,    # n_events
#     c_double_p,       # c_events_pointer
#     c_double_p,       # c_amplitudes_pointer
#     ctypes.c_double,  # windowLength
#     ctypes.c_double,  # delta
#     ctypes.c_ubyte,   # AR_ORDER
#     ctypes.c_bool,    # hasTheta0
#     ctypes.c_double,  # alpha
#     ctypes.c_double,  # beta
#     ctypes.c_uint,    # distribution
#     ctypes.c_uint,    # maxIter
#     ctypes.c_bool,    # serializeData
#     ctypes.c_char_p,  # outputDataName
#     ctypes.c_char_p   # outputTausName
# ]

