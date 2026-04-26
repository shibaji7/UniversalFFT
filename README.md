
# UniversalFFT

[![License: MIT](https://img.shields.io/badge/License%3A-MIT-green)](https://choosealicense.com/licenses/mit/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Documentation Status](https://img.shields.io/badge/docs-readthedocs-blue)](https://universalfft.readthedocs.io/en/latest/)
[![codecov](https://img.shields.io/badge/coverage-100%25-brightgreen)](https://github.com/shibaji7/UniversalFFT)

UniversalFFT is an open-source multi-language FFT/IFFT wrapper library that enforces the physically correct Fourier transform conventions described by Boteler (2012) for geoscience applications. All four language implementations — **Python**, **C**, **MATLAB**, and **R** — produce numerically identical results within a tolerance of 10⁻⁹.

## Why UniversalFFT?

Every FFT library makes silent choices about scaling factors and sign conventions. A `numpy.fft.ifft` and R's `fft(inverse=TRUE)` compute different things despite both calling themselves "inverse FFT". When applying a single transform (not a pair) — for example recovering a filter's impulse response — using the wrong convention gives a result that is numerically plausible but physically incorrect.

UniversalFFT anchors every choice to the Fourier integral in frequency *f* (not angular frequency ω), following Boteler (2012), so:

- Parseval's theorem holds without any 2π factors
- The impulse response of a filter integrates to 1
- A cosine of amplitude *A* produces frequency-domain spikes of *A/2* (not *A·Δt/2*)

## Quick Start

### Python

```bash
pip install -e python/
```

```python
import numpy as np
from universalfft import fft_cc, ifft_cc, freqs, fft_filter
from universalfft.utils import low_pass_response

N, dt = 2048, 60.0          # 2048 samples at 1-minute cadence (magnetometer)
x = ...                     # your geomagnetic time series

# Low-pass filter at 1-hour period
f = freqs(N, dt)
H = low_pass_response(f, fc=1/3600).astype(complex)
y = fft_filter(x, H, dt)

# Impulse response — MUST use CC inverse (Boteler §4.2)
h = ifft_cc(H, dt)
```

### C

```bash
cd c && make all
```

```c
#include "universalfft.h"
ufft_freqs(f, N, dt);
for (size_t k = 0; k < N; k++)
    H_re[k] = fabs(f[k]) <= fc ? 1.0 : 0.0;
ufft_filter(x_re, x_im, H_re, H_im, y_re, y_im, N, dt);
ufft_cc_inverse(H_re, H_im, h_re, h_im, N, dt);   /* impulse response */
```

### MATLAB

```matlab
addpath('matlab/')
f = ufft_freqs(N, dt);
H = double(abs(f) <= 1/3600);
y = ufft_filter(x, H, dt);
h = ufft_cc_inverse(complex(H), dt);   % impulse response
```

### R

```r
source("r/universalfft.R")
f <- ufft_freqs(N, dt)
H <- ufft_lowpass(f, 1/3600) + 0i
y <- ufft_filter(x, H, dt)
h <- ufft_cc_inverse(H, dt)            # impulse response
```

## Source Code

The library source code can be found on the [UniversalFFT GitHub](https://github.com/shibaji7/UniversalFFT) repository.

If you have any questions or concerns please submit an **Issue** on the [UniversalFFT GitHub](https://github.com/shibaji7/UniversalFFT) repository.

## Documentation

Read the docs: https://universalfft.readthedocs.io/en/latest/

## Reference

Boteler, D.H. (2012). *On Choosing Fourier Transforms for Practical Geoscience Applications.* International Journal of Geosciences, **3**, 952–959. doi:[10.4236/ijg.2012.325096](https://doi.org/10.4236/ijg.2012.325096)
