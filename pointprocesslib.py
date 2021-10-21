import ctypes
import numpy as np
from enum import Enum
import os


class InterEventDistribution(Enum):
    InverseGaussian = 0
    LogNormal = 1
    Gaussian = 2


def regrlikel(
        events: np.array,
        windowLength: float,
        delta: float,
        AR_ORDER: int,
        hasTheta0: bool = True,
        rightCensoring: bool = True,
        alpha: float = 0.02,
        distribution: InterEventDistribution = InterEventDistribution.InverseGaussian,
        maxIter: int = 1000,
        serializeData: bool = False,
        outputDataName: str = "Data.csv",
        outputTausName: str = "Taus.csv"
) -> None:
    assert len(events.shape) == 1
    n_events = len(events)
    c_events_pointer = events.astype(np.double).ctypes.data_as(c_double_p)
    # Compute result...
    cdll.regrlikel(
        n_events,
        c_events_pointer,
        windowLength,
        delta,
        AR_ORDER,
        hasTheta0,
        rightCensoring,
        alpha,
        distribution.value,
        maxIter,
        serializeData,
        outputDataName.encode('utf-8'),
        outputTausName.encode('utf-8')
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


current_path = os.path.dirname(os.path.realpath(__file__))
lib_file = os.path.join(current_path, "build", "libpointprocess.dylib")
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