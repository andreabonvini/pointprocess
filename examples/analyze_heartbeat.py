"""
Standalone example demonstrating the pointprocess library for heart rate
variability (HRV) analysis.

The script fits a point-process regression model to R-peak event times
shipped with the library (``events.npy``), computes spectral analysis,
and plots the results.

Usage:
    python analyze_heartbeat.py

Requirements (installed by the companion README):
    pointprocess, numpy, scipy, matplotlib, loguru
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from loguru import logger

from pointprocess import (
    Distributions,
    compute_full_regression,
    compute_single_regression,
    compute_spectral_analysis,
)

SCRIPT_DIR = Path(__file__).resolve().parent
EVENTS_NPY = SCRIPT_DIR / "events.npy"
OUTPUT_DIR = SCRIPT_DIR / "output"


# ---------------------------------------------------------------------------
# 1.  Point-process regression on a window of R-R events
# ---------------------------------------------------------------------------

def run_single_regression(events: np.ndarray) -> None:
    """Fit a single Inverse-Gaussian point-process model to a data window."""

    window = events[75:301]
    logger.info(
        "Single regression: {} events, span {:.1f}--{:.1f} s",
        len(window), window[0], window[-1],
    )

    t0 = time.perf_counter()
    result = compute_single_regression(
        events=window,
        ar_order=9,
        has_theta0=True,
        right_censoring=False,
        alpha=0.02,
        distribution=Distributions.InverseGaussian,
        max_iter=10000,
    )
    elapsed = time.perf_counter() - t0
    logger.info("Single regression completed in {:.4f} s", elapsed)
    logger.info("  theta_0={:.6f}  kappa={:.4f}  log-lik={:.2f}  converged={}",
                result.theta0, result.kappa, result.likelihood, result.converged)

    # -- Spectral analysis on the fitted model --------------------------------
    t0 = time.perf_counter()
    spectral = compute_spectral_analysis(
        thetap=result.thetap,
        mean_interval=result.mean_interval,
        variance=result.sigma ** 2,
    )
    elapsed = time.perf_counter() - t0
    logger.info("Spectral analysis completed in {:.4f} s", elapsed)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Left: RR intervals in the window
    ax1.plot(window[1:], 1000 * np.diff(window), color="steelblue", linewidth=0.8)
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("RR interval [ms]")
    ax1.set_title("RR Intervals (analysis window)")

    # Right: Power spectral density
    ax2.plot(spectral.frequencies, spectral.powers, color="black", linewidth=0.8)
    for pole in spectral.poles:
        ax2.axvline(pole.frequency, color="red", linestyle="--", alpha=0.3)
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Power [ms^2/Hz]")
    ax2.set_xlim(0, 0.5)
    ax2.set_title("Power Spectral Density (single window)")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "single_regression.png", dpi=150)
    logger.info("Saved single-regression figure")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 2.  Full (time-varying) regression over a longer segment
# ---------------------------------------------------------------------------

def run_full_regression(events: np.ndarray) -> None:
    """Run sliding-window regression and plot time-varying HRV indices."""

    logger.info(
        "Full regression: {} events, span {:.1f}--{:.1f} s",
        len(events), events[0], events[-1],
    )

    t0 = time.perf_counter()
    result = compute_full_regression(
        events=events.tolist(),
        window_length=60.0,
        delta=0.005,
        ar_order=9,
        has_theta0=True,
        right_censoring=True,
        alpha=0.02,
        distribution=Distributions.InverseGaussian,
        max_iter=1000,
    )
    elapsed = time.perf_counter() - t0
    logger.info("Full regression completed in {:.4f} s", elapsed)

    t0 = time.perf_counter()
    result.compute_hrv_indices()
    elapsed = time.perf_counter() - t0
    logger.info("HRV indices computed in {:.4f} s", elapsed)

    d = result.to_dict()

    fig, axes = plt.subplots(3, 1, figsize=(13, 10), sharex=True)

    # Mean interval (mu)
    axes[0].plot(d["Time"], 1000 * np.array(d["Mu"]), color="steelblue", linewidth=0.5)
    axes[0].set_ylabel("Mean RR [ms]")
    axes[0].set_title("Time-varying Mean Interval")

    # Hazard rate (lambda)
    axes[1].plot(d["Time"], d["lambda"], color="darkgreen", linewidth=0.5)
    axes[1].set_ylabel("Hazard rate")
    axes[1].set_title("Instantaneous Hazard Rate")

    # LF / HF ratio
    lf = np.array(d["powLF"])
    hf = np.array(d["powHF"])
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(hf > 0, lf / hf, np.nan)
    axes[2].plot(d["Time"], ratio, color="crimson", linewidth=0.5)
    axes[2].set_ylabel("LF / HF")
    axes[2].set_xlabel("Time [s]")
    axes[2].set_title("Sympatho-vagal Balance (LF/HF)")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "full_regression.png", dpi=150)
    logger.info("Saved full-regression figure")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)

    total_t0 = time.perf_counter()

    if not EVENTS_NPY.exists():
        logger.error("events.npy not found at {}. Cannot proceed.", EVENTS_NPY)
        return

    events = np.load(str(EVENTS_NPY))
    logger.info("Loaded {} R-peak event times from events.npy", len(events))

    run_single_regression(events)
    run_full_regression(events)

    total_elapsed = time.perf_counter() - total_t0
    logger.info("Total processing time: {:.3f} s", total_elapsed)
    logger.info("Output saved to {}", OUTPUT_DIR)


if __name__ == "__main__":
    main()
