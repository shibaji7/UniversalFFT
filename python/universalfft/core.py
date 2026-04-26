"""
Core FFT/IFFT wrappers implementing Boteler (2012) conventions.

NumPy reference conventions (what we wrap)
------------------------------------------
np.fft.fft(x)  → X_raw[k] = Σ_{n=0}^{N-1} x[n] e^{-i2πkn/N}   (no scaling)
np.fft.ifft(X) → x_raw[n] = (1/N) Σ_{k=0}^{N-1} X[k] e^{+i2πkn/N}

Boteler CC forward  = np.fft.fft(x)  * dt
Boteler CC inverse  = np.fft.ifft(X) * N * df   (= ifft(X) / dt,  since N*df=1/dt)
Boteler CD forward  = np.fft.fft(x)  / N        (= ifft-scaled forward)
Boteler CD inverse  = np.fft.ifft(X) * N        (= raw summation, no 1/N)
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def freqs(N: int, dt: float) -> NDArray[np.float64]:
    """Return the FFT frequency bin array (Hz) for N samples spaced dt seconds.

    Follows NumPy ordering: [0, df, 2df, …, f_Nyquist, -f_Nyquist+df, …, -df].
    Use ``np.fft.fftshift`` to reorder to [-f_Nyquist, …, f_Nyquist].

    Parameters
    ----------
    N:
        Number of samples.
    dt:
        Sampling interval in seconds.

    Returns
    -------
    f : ndarray, shape (N,)
        Frequency bins in Hz.
    """
    return np.fft.fftfreq(N, d=dt)


# ---------------------------------------------------------------------------
# Continuous–Continuous (CC) pair  —  Boteler Eqs. (21a)/(21b)
# ---------------------------------------------------------------------------

def fft_cc(
    x: ArrayLike,
    dt: float,
) -> NDArray[np.complexfloating]:
    """Continuous–Continuous forward DFT (Boteler Eq. 21a).

    X(k) = Σ_{n=0}^{N-1} x(n) e^{-i2πkn/N} Δt

    The factor Δt makes X(k) approximate the Fourier integral:
        F(f) = ∫ f(t) e^{-i2πft} dt

    Parameters
    ----------
    x:
        Time-domain samples (real or complex), shape (N,).
    dt:
        Sampling interval in seconds (Δt).

    Returns
    -------
    X : ndarray, shape (N,), complex
        Frequency-domain values in the FFT-output order
        (0, +f, …, f_Nyquist, -f_Nyquist+df, …, -df).
    """
    x = np.asarray(x, dtype=complex)
    return np.fft.fft(x) * dt


def ifft_cc(
    X: ArrayLike,
    dt: float,
) -> NDArray[np.complexfloating]:
    """Continuous–Continuous inverse DFT (Boteler Eq. 21b).

    x(n) = Σ_{k=0}^{N-1} X(k) e^{+i2πkn/N} Δf

    where Δf = 1/(N Δt), so the scaling factor is Δf = 1/(N dt).

    Parameters
    ----------
    X:
        Frequency-domain values (FFT-output order), shape (N,).
    dt:
        Sampling interval in seconds of the *original* time series.

    Returns
    -------
    x : ndarray, shape (N,), complex
        Reconstructed time-domain samples.
    """
    X = np.asarray(X, dtype=complex)
    N = X.shape[0]
    df = 1.0 / (N * dt)
    # np.fft.ifft multiplies by 1/N; we need to multiply by Δf instead,
    # so we undo the 1/N and apply Δf: factor = N * df = 1/dt
    return np.fft.ifft(X) * N * df


# ---------------------------------------------------------------------------
# Continuous–Discrete (CD) pair  —  Boteler Eqs. (22a)/(22b)
# ---------------------------------------------------------------------------

def fft_cd(
    x: ArrayLike,
    dt: float | None = None,   # dt unused for CD forward but kept for API symmetry
) -> NDArray[np.complexfloating]:
    """Continuous–Discrete forward DFT (Boteler Eq. 22a).

    X(k) = (1/N) Σ_{n=0}^{N-1} x(n) e^{-i2πkn/N}

    Appropriate when the time-domain signal is continuous but the
    frequency-domain representation is a set of discrete components
    (e.g. Fourier harmonics of a periodic waveform).  Taking the DFT
    of a cosine of amplitude A yields spikes of amplitude A/2 at ±f₁,
    which combine via Euler's formula to recover amplitude A.

    Parameters
    ----------
    x:
        Time-domain samples, shape (N,).
    dt:
        Unused; present for API symmetry with ``fft_cc``.

    Returns
    -------
    X : ndarray, shape (N,), complex
    """
    x = np.asarray(x, dtype=complex)
    N = x.shape[0]
    return np.fft.fft(x) / N


def ifft_cd(
    X: ArrayLike,
    dt: float | None = None,
) -> NDArray[np.complexfloating]:
    """Continuous–Discrete inverse DFT (Boteler Eq. 22b).

    x(n) = Σ_{k=0}^{N-1} X(k) e^{+i2πkn/N}

    Parameters
    ----------
    X:
        Frequency-domain values (CD-scaled), shape (N,).
    dt:
        Unused; present for API symmetry with ``ifft_cc``.

    Returns
    -------
    x : ndarray, shape (N,), complex
    """
    X = np.asarray(X, dtype=complex)
    N = X.shape[0]
    # np.fft.ifft applies 1/N; undo it to get raw summation
    return np.fft.ifft(X) * N


# ---------------------------------------------------------------------------
# Convenience: filter in frequency domain (either convention works)
# ---------------------------------------------------------------------------

def fft_filter(
    x: ArrayLike,
    H: ArrayLike,
    dt: float,
) -> NDArray[np.complexfloating]:
    """Apply a linear filter H(f) to time series x via frequency-domain multiplication.

    Uses the CC convention internally, but because the combined scaling
    of a forward+inverse pair is Δt·Δf = 1/N regardless of convention
    (Boteler Section 4.1), either convention produces identical results.

    Parameters
    ----------
    x:
        Input time-domain samples, shape (N,).
    H:
        Transfer function evaluated at the FFT frequency bins (same order
        as returned by ``freqs``), shape (N,).
    dt:
        Sampling interval in seconds.

    Returns
    -------
    y : ndarray, shape (N,), complex
        Filtered time-domain signal.
    """
    x = np.asarray(x, dtype=complex)
    H = np.asarray(H, dtype=complex)
    if x.shape != H.shape:
        raise ValueError(
            f"x and H must have the same shape; got {x.shape} vs {H.shape}"
        )
    X = fft_cc(x, dt)
    Y = X * H
    return ifft_cc(Y, dt)
