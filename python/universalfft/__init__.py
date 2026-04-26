"""
UniversalFFT — Boteler (2012) compliant FFT/IFFT wrappers.

Two transform pairs are provided following the conventions in:
  Boteler, D.H. (2012). On Choosing Fourier Transforms for Practical
  Geoscience Applications. International Journal of Geosciences, 3, 952-959.

Convention summary
------------------
Forward transform uses e^{-i2πft}; inverse uses e^{+i2πft}.
Frequency variable is *f* (Hz), never angular frequency ω, so the
factor 1/(2π) never appears (Parseval's theorem holds with integration
over f).

Two pairs
---------
CC (Continuous–Continuous, Eq. 21 in paper)
    Both time-domain and frequency-domain values are treated as
    samples of continuous functions.

    X(k) = Σ x(n) e^{-i2πkn/N} Δt          (forward, scaled by Δt)
    x(n) = Σ X(k) e^{+i2πkn/N} Δf          (inverse, scaled by Δf)

    Use for: impulse-response recovery (single inverse transform).

CD (Continuous–Discrete, Eq. 22 in paper)
    Time domain is continuous; frequency domain is a set of discrete
    components (e.g. harmonics of a periodic signal).

    X(k) = (1/N) Σ x(n) e^{-i2πkn/N}       (forward, scaled by 1/N)
    x(n) =       Σ X(k) e^{+i2πkn/N}       (inverse, unscaled)

    Use for: spectrum determination of periodic signals.

For filter applications that use a *pair* of transforms (FFT →
multiply → IFFT) either convention gives identical results because
Δt·Δf = 1/N cancels the difference in scaling.
"""

from .core import (
    fft_cc,
    ifft_cc,
    fft_cd,
    ifft_cd,
    fft_filter,
    freqs,
)
from .utils import next_power_of_two, pad_to_length

__version__ = "1.0.0"
__author__ = "UniversalFFT contributors"

__all__ = [
    "fft_cc",
    "ifft_cc",
    "fft_cd",
    "ifft_cd",
    "fft_filter",
    "freqs",
    "next_power_of_two",
    "pad_to_length",
]
