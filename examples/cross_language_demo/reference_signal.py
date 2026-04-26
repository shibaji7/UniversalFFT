"""
Shared signal definition for the cross-language demo.

Every language reads/generates the same signal so results are comparable.
The signal is a 60 Hz cosine plus a 180 Hz (3rd harmonic) component,
simulating a GIC-distorted transformer magnetising current waveform.
"""

import numpy as np

# Signal parameters — identical across all language demos
N   = 256          # samples (power of two for efficiency)
dt  = 1e-3         # 1 ms  →  fs = 1 kHz
f1  = 60.0         # Hz  — fundamental
f3  = 180.0        # Hz  — 3rd harmonic
A1  = 1.0          # amplitude of fundamental
A3  = 0.3          # amplitude of harmonic (GIC distortion level)
fc  = 100.0        # LP filter cut-off (retains both components)


def make_signal():
    t = np.arange(N) * dt
    x = A1 * np.cos(2 * np.pi * f1 * t) + A3 * np.cos(2 * np.pi * f3 * t)
    return t, x
