<!--
Author(s): Shibaji Chakraborty
-->

# Cross-Language Demo

This demo processes the same reference signal through all 11 language
implementations and verifies that every output matches the Python reference
within a tolerance of \(10^{-9}\).

---

## Overview

The six-panel figure below summarises every verification step: input signal,
CC forward magnitude, CD spectrum spikes, round-trip reconstruction, round-trip
error (below \(10^{-9}\)), and the cross-language max-error bar chart:

![Cross-language demo overview](../assets/figures/cross_language_demo.png)

---

## Signal definition

All languages use the same signal parameters:

| Parameter | Value |
|-----------|-------|
| \(N\) | 128 samples |
| \(\Delta t\) | 1 ms (fs = 1 kHz) |
| Frequency \(f_1\) | 60 Hz, amplitude 2.5 |
| LP impulse cut-off \(f_c\) | 1/3600 Hz (1-hour period) |
| Impulse response \(N_2\) | 256 samples, \(\Delta t=60\) s |

---

## Step 1 — Generate reference vectors (Python)

```bash
cd UniversalFFT
python tests/cross_language/generate_reference.py
```

This writes three `.npz` files to `tests/data/`:

- `cosine_vectors.npz` — CC forward/inverse of the cosine signal
- `filter_vectors.npz` — LP filter output
- `impulse_response_vectors.npz` — CC inverse of boxcar TF

---

## Step 2 — Run each language demo

=== "Python"

    ```bash
    python examples/cross_language_demo/run_all.py
    ```

=== "C"

    ```bash
    make c-test
    # writes: tests/data/c_X_cc_cosine.csv
    #         tests/data/c_x_cc_rec_cosine.csv
    #         tests/data/c_h_impulse.csv
    ```

=== "C++"

    ```bash
    make cpp-test
    # writes: tests/data/cpp_X_cc_cosine.csv
    #         tests/data/cpp_x_cc_rec_cosine.csv
    #         tests/data/cpp_h_impulse.csv
    ```

=== "Fortran"

    ```bash
    make fortran-test
    # writes: tests/data/fortran_X_cc_cosine.csv
    #         tests/data/fortran_x_cc_rec_cosine.csv
    #         tests/data/fortran_h_impulse.csv
    ```

=== "Julia"

    ```bash
    make julia-test
    # writes: tests/data/julia_X_cc_cosine.csv
    #         tests/data/julia_x_cc_rec_cosine.csv
    #         tests/data/julia_h_impulse.csv
    ```

=== "Rust"

    ```bash
    make rust-test
    # writes: tests/data/rust_X_cc_cosine.csv
    #         tests/data/rust_x_cc_rec_cosine.csv
    #         tests/data/rust_h_impulse.csv
    ```

=== "MATLAB"

    ```matlab
    % From MATLAB command window (add matlab/ to path first)
    addpath('matlab/')
    ufft_demo('../tests/data')
    % writes tests/data/matlab_*.csv
    ```

=== "R"

    ```bash
    Rscript r/ufft_demo.R
    # writes tests/data/r_*.csv
    ```

=== "JavaScript"

    ```bash
    make js-test
    # writes: tests/data/js_X_cc_cosine.csv
    #         tests/data/js_x_cc_rec_cosine.csv
    #         tests/data/js_h_impulse.csv
    ```

=== "Octave"

    ```bash
    make octave-test
    # writes: tests/data/octave_X_cc_cosine.csv
    #         tests/data/octave_x_cc_rec_cosine.csv
    #         tests/data/octave_h_impulse.csv
    ```

=== "CUDA / HIP"

    ```bash
    make -C cuda test           # NVIDIA CUDA
    make -C cuda HIP=1 test     # AMD HIP
    # writes: tests/data/cuda_X_cc_cosine.csv
    #         tests/data/cuda_x_cc_rec_cosine.csv
    #         tests/data/cuda_h_impulse.csv
    ```

=== "IDL / GDL"

    ```bash
    # GDL (free runtime — install: sudo apt install gnudatalanguage)
    gdl idl/universalfft.pro
    # Commercial IDL:
    # idl -e "@idl/universalfft.pro"
    ```

---

## Step 3 — Validate

```bash
python tests/cross_language/validate_all.py
```

Expected output (after all demos have run):

```
=================================================================
UniversalFFT cross-language validation
Tolerance: 1e-09
=================================================================

[1] CC forward transform of cosine signal
  [PASS] C:          max_abs_error = 1.93e-16  (tol=1e-09)
  [PASS] Fortran:    max_abs_error = 1.93e-16  (tol=1e-09)
  [PASS] C++:        max_abs_error = 1.93e-16  (tol=1e-09)
  [PASS] JavaScript: max_abs_error = 1.93e-16  (tol=1e-09)
  [PASS] Rust:       max_abs_error = 1.93e-16  (tol=1e-09)
  [PASS] Julia:      max_abs_error = 1.93e-16  (tol=1e-09)
  [PASS] MATLAB:     max_abs_error = 1.93e-16  (tol=1e-09)
  [PASS] R:          max_abs_error = 1.93e-16  (tol=1e-09)
  [PASS] Octave:     max_abs_error = 1.93e-16  (tol=1e-09)
  [PASS] CUDA:       max_abs_error = 1.93e-16  (tol=1e-09)

[2] CC round-trip (forward → inverse) on cosine signal
  ...

[3] Impulse response (CC inverse of boxcar TF)
  ...

=================================================================
OVERALL: PASS
=================================================================
```

!!! note "Missing language CSVs are SKIPped, not FAILed"
    The validator prints `[SKIP]` for any CSV that doesn't exist yet.
    The overall result is still `PASS` unless a *present* file exceeds the tolerance.

!!! note "MATLAB/R RNG differs from Python"
    MATLAB's `rng(0)` and Python's `default_rng(0)` use different RNG algorithms,
    so the random filter test signals differ between languages. The validation
    compares each language against the Python-generated `.npz` reference, which
    is the authoritative ground truth.

---

## Normalisation map

Each language requires a different correction to map its native FFT output
to the Boteler CC/CD conventions:

| Language | Native forward | CC fwd | CC inv | CD fwd | CD inv |
|----------|---------------|--------|--------|--------|--------|
| Python (NumPy) | raw sum | `× Δt` | `× N·Δf` | `/ N` | `× N` |
| C (self-contained) | raw sum | `× Δt` | `× Δf` | `/ N` | raw |
| C++ (self-contained) | raw sum | `× Δt` | `× Δf` | `/ N` | raw |
| Fortran (self-contained) | raw sum | `× Δt` | `× Δf` | `/ N` | raw |
| Julia (FFTW.jl) | raw sum | `× Δt` | `× N·Δf` | `/ N` | `× N` |
| Rust (rustfft) | raw sum | `× Δt` | `× Δf` | `/ N` | raw |
| MATLAB | raw sum | `× Δt` | `× N·Δf` | `/ N` | `× N` |
| R | raw sum (inv: **also raw**) | `× Δt` | `× Δf` | `/ N` | raw |
| JavaScript (self-contained) | raw sum | `× Δt` | `× Δf` | `/ N` | raw |
| Octave | raw sum | `× Δt` | `× N·Δf` | `/ N` | `× N` |
| CUDA/HIP (cuFFT/hipFFT) | raw sum | `× Δt` | `× Δf` | `/ N` | raw |
| IDL / GDL | `(1/N)·raw` ← **reversed!** | `× N·Δt` | `× Δf` | raw | raw |

!!! warning "IDL / GDL forward convention is reversed"
    IDL's `FFT(x,-1)` already divides by N in the **forward** direction —
    the opposite of every other library. `universalfft.pro` corrects for this
    transparently.

!!! info "R inverse is uniquely raw"
    R's `fft(inverse=TRUE)` returns the raw sum with **no** `1/N` — unlike
    NumPy, MATLAB, Julia, and Octave which all divide by N on inverse.
    The R CC inverse therefore multiplies only by `Δf` (not `N·Δf`).
