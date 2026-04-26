"""Utility helpers for UniversalFFT."""

from __future__ import annotations

import math
import numpy as np
from numpy.typing import ArrayLike, NDArray


def next_power_of_two(n: int) -> int:
    """Return the smallest power of two >= n."""
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    return 1 << (n - 1).bit_length()


def pad_to_length(
    x: ArrayLike,
    length: int,
    mode: str = "zero",
) -> NDArray:
    """Zero-pad or reflect-pad x to the given length along axis 0.

    Parameters
    ----------
    x:
        Input array.
    length:
        Target length (must be >= len(x)).
    mode:
        ``"zero"`` (default) or ``"reflect"``.

    Returns
    -------
    x_padded : ndarray, shape (length, ...)
    """
    x = np.asarray(x)
    n = x.shape[0]
    if length < n:
        raise ValueError(
            f"Target length {length} is shorter than input length {n}."
        )
    if length == n:
        return x.copy()

    pad_width = [(0, length - n)] + [(0, 0)] * (x.ndim - 1)
    if mode == "zero":
        return np.pad(x, pad_width, mode="constant", constant_values=0)
    elif mode == "reflect":
        return np.pad(x, pad_width, mode="reflect")
    else:
        raise ValueError(f"Unknown mode {mode!r}; choose 'zero' or 'reflect'.")


def low_pass_response(f: ArrayLike, fc: float) -> NDArray[np.float64]:
    """Build a brick-wall low-pass transfer function on frequency array f.

    H(f) = 1 for |f| <= fc, 0 otherwise.

    Parameters
    ----------
    f:
        Frequency bin array as returned by ``universalfft.freqs``.
    fc:
        Cut-off frequency in Hz.

    Returns
    -------
    H : ndarray, shape like f, real-valued (0.0 or 1.0)
    """
    f = np.asarray(f, dtype=float)
    return (np.abs(f) <= fc).astype(float)


def sinc_integral_check(h: ArrayLike, dt: float, tol: float = 1e-6) -> bool:
    """Verify that the impulse response integrates to 1.0 (Boteler Section 4.2).

    Parameters
    ----------
    h:
        Impulse-response values (samples of a continuous function).
    dt:
        Sampling interval in seconds.
    tol:
        Absolute tolerance for the integral check.

    Returns
    -------
    bool
        True if ``|Σ h[n]·Δt − 1| < tol``.
    """
    h = np.asarray(h)
    integral = np.sum(np.real(h)) * dt
    return bool(abs(integral - 1.0) < tol)
