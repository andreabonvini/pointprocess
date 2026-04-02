# Examples

A standalone Python script that demonstrates the `pointprocess` library on
R-peak event times bundled with the repository (`events.npy`).

## Quick start

### 1. Create a virtual environment

```bash
cd examples
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
```

### 2. Install dependencies

```bash
pip install pointprocess numpy scipy matplotlib loguru
```

If you are developing locally and want to use the library from source, replace
`pip install pointprocess` with a local install from the repository root:

```bash
pip install ..   # from the examples/ directory
```

### 3. Run the example

```bash
python analyze_heartbeat.py
```

Figures are saved to `examples/output/`.

## What the script does

1. **Single-window regression** -- fits an Inverse-Gaussian point-process model
   (AR order 9) to a window of R-R event times and computes the power spectral
   density.
2. **Full (time-varying) regression** -- slides a 60 s window across the full
   recording, estimates HRV indices, and plots the mean RR interval, hazard
   rate, and LF/HF ratio over time.

All steps are profiled and logged with `loguru`.
