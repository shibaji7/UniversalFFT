#!/usr/bin/env python3
"""
Language-Agnostic Demo — Python driver.

Processes the reference signal through all four operations and prints
a comparison table.  Call from the repo root:

    python examples/cross_language_demo/run_all.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../python"))

import numpy as np
from universalfft import fft_cc, ifft_cc, fft_cd, freqs, fft_filter
from universalfft.utils import low_pass_response, sinc_integral_check
from reference_signal import N, dt, f1, f3, A1, A3, fc, make_signal

t, x = make_signal()
f    = freqs(N, dt)
H_lp = low_pass_response(f, fc).astype(complex)

# --- 1. CC forward transform ---
X_cc = fft_cc(x, dt)

# --- 2. CD spectrum (amplitudes at f1, f3) ---
X_cd = fft_cd(x)
idx1 = np.argmin(np.abs(f - f1))
idx3 = np.argmin(np.abs(f - f3))

# --- 3. LP filter ---
y_lp = fft_filter(x, H_lp, dt)

# --- 4. Impulse response (CC inverse of LP TF) ---
H_ir = low_pass_response(f, fc).astype(complex)
h    = ifft_cc(H_ir, dt)
ok   = sinc_integral_check(h, dt)

# --- Report ---
print("=" * 62)
print("UniversalFFT — Cross-language demo (Python reference)")
print("=" * 62)
print(f"Signal: {A1}·cos(2π·{f1}·t) + {A3}·cos(2π·{f3}·t)")
print(f"N={N}, dt={dt} s, fs={1/dt:.0f} Hz")
print()
print("1. CC forward: X_CC[0] (DC) =", f"{X_cc[0].real:.6e}")
print()
print("2. CD spectrum (Boteler §4.3):")
print(f"   |X_CD[f={f1:.0f} Hz]| = {abs(X_cd[idx1]):.6f}  (expect {A1/2:.4f} = A1/2)")
print(f"   |X_CD[f={f3:.0f} Hz]| = {abs(X_cd[idx3]):.6f}  (expect {A3/2:.4f} = A3/2)")
print()
print("3. LP-filtered signal:")
print(f"   max |y_lp - x| = {np.max(np.abs(np.real(y_lp) - x)):.2e}")
print(f"   (fc={fc} Hz passes 60 Hz but removes 180 Hz; error reflects 3rd harmonic removal + Gibbs ringing)")
print()
print("4. Impulse response (CC inverse):")
print(f"   Integral = {np.sum(np.real(h)) * dt:.8f}  (expect 1.0)")
print(f"   sinc_integral_check passed: {ok}")
print("=" * 62)
print()
print("To validate C/MATLAB/R outputs against this reference:")
print("  python tests/cross_language/generate_reference.py")
print("  make -C c test")
print("  # matlab -batch ufft_demo   (from matlab/)")
print("  # Rscript r/ufft_demo.R")
print("  python tests/cross_language/validate_all.py")
