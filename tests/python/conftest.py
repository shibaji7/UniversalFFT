"""Shared fixtures for UniversalFFT tests."""

import numpy as np
import pytest

@pytest.fixture
def simple_cosine():
    """A 60 Hz cosine sampled at 1 kHz for 0.1 s — N=100 samples."""
    N = 100
    dt = 1e-3          # 1 ms  →  fs = 1000 Hz
    f1 = 60.0          # Hz
    A = 3.0
    t = np.arange(N) * dt
    x = A * np.cos(2 * np.pi * f1 * t)
    return {"x": x, "dt": dt, "N": N, "f1": f1, "A": A, "t": t}


@pytest.fixture
def magnetic_storm_like():
    """Synthetic band-limited 'magnetic storm' time series, N=2048, dt=60 s."""
    rng = np.random.default_rng(42)
    N = 2048
    dt = 60.0
    # white noise then low-pass by zeroing upper half of spectrum
    raw = rng.standard_normal(N)
    F = np.fft.fft(raw)
    F[N // 4:3 * N // 4] = 0          # crude LP
    x = np.real(np.fft.ifft(F))
    return {"x": x, "dt": dt, "N": N}
