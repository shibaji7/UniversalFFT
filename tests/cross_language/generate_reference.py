#!/usr/bin/env python3
"""
Generate reference test vectors for cross-language validation.

Writes NumPy .npz files to tests/data/ that the C, MATLAB, and R wrappers
must reproduce within the agreed tolerance of 1e-9.

Usage
-----
    cd UniversalFFT
    python tests/cross_language/generate_reference.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../python"))

import numpy as np
from universalfft import fft_cc, ifft_cc, fft_cd, ifft_cd, freqs
from universalfft.utils import low_pass_response

DATA_DIR = os.path.join(os.path.dirname(__file__), "../data")


def generate_cosine_vectors():
    """60 Hz cosine, N=128, dt=1e-3 s."""
    N, dt, f1, A = 128, 1e-3, 60.0, 2.5
    t = np.arange(N) * dt
    x = A * np.cos(2 * np.pi * f1 * t)
    f = freqs(N, dt)

    X_cc = fft_cc(x, dt)
    x_cc_rec = ifft_cc(X_cc, dt)
    X_cd = fft_cd(x, dt)
    x_cd_rec = ifft_cd(X_cd, dt)

    np.savez(
        os.path.join(DATA_DIR, "cosine_vectors.npz"),
        x=x,
        t=t,
        f=f,
        dt=np.array([dt]),
        N=np.array([N]),
        X_cc_real=np.real(X_cc),
        X_cc_imag=np.imag(X_cc),
        x_cc_rec_real=np.real(x_cc_rec),
        x_cc_rec_imag=np.imag(x_cc_rec),
        X_cd_real=np.real(X_cd),
        X_cd_imag=np.imag(X_cd),
        x_cd_rec_real=np.real(x_cd_rec),
        x_cd_rec_imag=np.imag(x_cd_rec),
    )
    print(f"[OK] cosine_vectors.npz  (N={N}, dt={dt}, f1={f1} Hz, A={A})")


def generate_filter_vectors():
    """Magnetic-storm LP filter, N=256, dt=60 s."""
    rng = np.random.default_rng(0)
    N, dt = 256, 60.0
    x = rng.standard_normal(N)
    f = freqs(N, dt)
    fc = 1.0 / 3600.0          # 1 hr
    H = low_pass_response(f, fc).astype(complex)

    X_cc = fft_cc(x, dt)
    Y_cc = X_cc * H
    y_filtered = ifft_cc(Y_cc, dt)

    np.savez(
        os.path.join(DATA_DIR, "filter_vectors.npz"),
        x=x,
        f=f,
        H_real=np.real(H),
        H_imag=np.imag(H),
        dt=np.array([dt]),
        N=np.array([N]),
        X_cc_real=np.real(X_cc),
        X_cc_imag=np.imag(X_cc),
        y_filtered_real=np.real(y_filtered),
        y_filtered_imag=np.imag(y_filtered),
    )
    print(f"[OK] filter_vectors.npz  (N={N}, dt={dt}, fc={fc:.6f} Hz)")


def generate_impulse_response_vectors():
    """Sinc impulse response of boxcar LP filter, N=256, dt=60 s."""
    N, dt = 256, 60.0
    f = freqs(N, dt)
    fc = 1.0 / 3600.0
    H = low_pass_response(f, fc).astype(complex)
    h = ifft_cc(H, dt)          # CC inverse — correct choice (Boteler §4.2)

    np.savez(
        os.path.join(DATA_DIR, "impulse_response_vectors.npz"),
        f=f,
        H_real=np.real(H),
        H_imag=np.imag(H),
        dt=np.array([dt]),
        N=np.array([N]),
        h_real=np.real(h),
        h_imag=np.imag(h),
    )
    integral = np.sum(np.real(h)) * dt
    print(f"[OK] impulse_response_vectors.npz  (integral={integral:.8f}, expect≈1.0)")


if __name__ == "__main__":
    os.makedirs(DATA_DIR, exist_ok=True)
    generate_cosine_vectors()
    generate_filter_vectors()
    generate_impulse_response_vectors()
    print("\nAll reference vectors written to tests/data/")
