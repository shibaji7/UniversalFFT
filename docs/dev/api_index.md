<!--
Author(s): Shibaji Chakraborty
-->

# API Documentation Overview

## Symbol legend

<span class="api-badge api-package">P</span> Package &nbsp;
<span class="api-badge api-class">C</span> Class &nbsp;
<span class="api-badge api-method">M</span> Function / Method

---

## Python package — `universalfft`

| Symbol | Type | Description |
|--------|:----:|-------------|
| `universalfft` | <span class="api-badge api-package">P</span> | Top-level package |
| `fft_cc` | <span class="api-badge api-method">M</span> | CC forward DFT (Boteler Eq. 21a) |
| `ifft_cc` | <span class="api-badge api-method">M</span> | CC inverse DFT (Boteler Eq. 21b) |
| `fft_cd` | <span class="api-badge api-method">M</span> | CD forward DFT (Boteler Eq. 22a) |
| `ifft_cd` | <span class="api-badge api-method">M</span> | CD inverse DFT (Boteler Eq. 22b) |
| `fft_filter` | <span class="api-badge api-method">M</span> | Frequency-domain filter via CC pair |
| `freqs` | <span class="api-badge api-method">M</span> | FFT frequency bin array (Hz) |
| `universalfft.utils` | <span class="api-badge api-package">P</span> | Utility helpers |
| `next_power_of_two` | <span class="api-badge api-method">M</span> | Smallest power of two ≥ n |
| `pad_to_length` | <span class="api-badge api-method">M</span> | Zero-pad or reflect-pad array |
| `low_pass_response` | <span class="api-badge api-method">M</span> | Brick-wall LP transfer function |
| `sinc_integral_check` | <span class="api-badge api-method">M</span> | Verify impulse response integrates to 1 |

---

## C library — `universalfft.h`

| Function | Description |
|----------|-------------|
| `ufft_inplace` | Core Cooley–Tukey radix-2 DIT FFT (in-place) |
| `ufft_cc_forward` | CC forward: \(X[k] = \sum x[n] e^{-i2\pi kn/N} \Delta t\) |
| `ufft_cc_inverse` | CC inverse: \(x[n] = \sum X[k] e^{+i2\pi kn/N} \Delta f\) |
| `ufft_cd_forward` | CD forward: \(X[k] = (1/N)\sum x[n] e^{-i2\pi kn/N}\) |
| `ufft_cd_inverse` | CD inverse: raw summation |
| `ufft_freqs` | Frequency bin array (Hz) |
| `ufft_filter` | LP/BP filter via CC pair |

---

## MATLAB functions

| Function | Description |
|----------|-------------|
| `ufft_cc_forward(x, dt)` | CC forward |
| `ufft_cc_inverse(X, dt)` | CC inverse |
| `ufft_cd_forward(x)` | CD forward |
| `ufft_cd_inverse(X)` | CD inverse |
| `ufft_freqs(N, dt)` | Frequency bin array |
| `ufft_filter(x, H, dt)` | Frequency-domain filter |

---

## R functions

| Function | Description |
|----------|-------------|
| `ufft_cc_forward(x, dt)` | CC forward |
| `ufft_cc_inverse(X, dt)` | CC inverse |
| `ufft_cd_forward(x)` | CD forward |
| `ufft_cd_inverse(X)` | CD inverse |
| `ufft_freqs(N, dt)` | Frequency bin array |
| `ufft_filter(x, H, dt)` | Frequency-domain filter |
| `ufft_lowpass(f, fc)` | Brick-wall LP transfer function |

---

## Julia — `UniversalFFT.jl`

| Function | Description |
|----------|-------------|
| `fft_cc(x, dt)` | CC forward |
| `ifft_cc(X, dt)` | CC inverse |
| `fft_cd(x)` | CD forward |
| `ifft_cd(X)` | CD inverse |
| `freqs(N, dt)` | Frequency bin array |
| `low_pass_response(f, fc)` | Brick-wall LP response |
| `fft_filter(x, H, dt)` | Frequency-domain filter |

---

## Fortran 90 module — `universalfft_mod`

| Subroutine | Description |
|------------|-------------|
| `ufft_cc_forward(x, N, dt, X)` | CC forward |
| `ufft_cc_inverse(X, N, dt, x)` | CC inverse |
| `ufft_cd_forward(x, N, X)` | CD forward |
| `ufft_cd_inverse(X, N, x)` | CD inverse |
| `ufft_freqs(N, dt, f)` | Frequency bin array |
| `ufft_filter(x, H, N, dt, y)` | Frequency-domain filter |

---

## C++ — `namespace ufft`

| Function | Description |
|----------|-------------|
| `fft_cc(x, dt)` / `fft_cc_real(x, dt)` | CC forward |
| `ifft_cc(X, dt)` | CC inverse |
| `fft_cd(x)` / `fft_cd_real(x)` | CD forward |
| `ifft_cd(X)` | CD inverse |
| `freqs(N, dt)` | Frequency bin array |
| `fft_filter(x, H, dt)` | Frequency-domain filter |

---

## Rust crate — `universalfft`

| Function | Description |
|----------|-------------|
| `fft_cc(x, dt)` | CC forward |
| `ifft_cc(X, dt)` | CC inverse |
| `fft_cd(x)` | CD forward |
| `ifft_cd(X)` | CD inverse |
| `freqs(n, dt)` | Frequency bin array |
| `low_pass_response(f, fc)` | Brick-wall LP response |
| `fft_filter(x, H, dt)` | Frequency-domain filter |
| `sinc_integral_check(h, dt, tol)` | Impulse response normalisation check |

---

## CUDA / HIP — host API

| Function | Description |
|----------|-------------|
| `ufft_cc_forward_host(…, dt)` | CC forward (GPU accelerated) |
| `ufft_cc_inverse_host(…, dt)` | CC inverse (GPU accelerated) |
| `ufft_cd_forward_host(…)` | CD forward (GPU accelerated) |
| `ufft_cd_inverse_host(…)` | CD inverse (GPU accelerated) |
| `ufft_freqs_cpu(f, N, dt)` | Frequency bin array (CPU side) |
| `ufft_filter_host(…, dt)` | Frequency-domain filter (GPU accelerated) |

---

## JavaScript ES module — `universalfft.js`

| Function | Description |
|----------|-------------|
| `fftCC(re, im, dt)` / `fftCCReal(x, dt)` | CC forward |
| `ifftCC(X, dt)` | CC inverse |
| `fftCD(re, im)` / `fftCDReal(x)` | CD forward |
| `ifftCD(X)` | CD inverse |
| `freqs(N, dt)` | Frequency bin array |
| `lowPassResponse(f, fc)` | Brick-wall LP response |
| `fftFilter(x, H, dt)` | Frequency-domain filter |

---

## Octave / MATLAB — `universalfft.m`

| Function | Description |
|----------|-------------|
| `ufft_cc_forward(x, dt)` | CC forward |
| `ufft_cc_inverse(X, dt)` | CC inverse |
| `ufft_cd_forward(x)` | CD forward |
| `ufft_cd_inverse(X)` | CD inverse |
| `ufft_freqs(N, dt)` | Frequency bin array |
| `ufft_low_pass(f, fc)` | Brick-wall LP response |
| `ufft_filter(x, H, dt)` | Frequency-domain filter |

---

## IDL / GDL — `universalfft.pro`

!!! note "GDL — free IDL runtime"
    Use `gdl universalfft.pro` for testing without a commercial IDL licence.
    See [IDL/GDL docs](idl/index.md) for installation instructions.

| Function | Description |
|----------|-------------|
| `ufft_cc_forward(x, dt)` | CC forward (corrects IDL's reversed 1/N convention) |
| `ufft_cc_inverse(X, dt)` | CC inverse |
| `ufft_cd_forward(x)` | CD forward |
| `ufft_cd_inverse(X)` | CD inverse |
| `ufft_freqs(N, dt)` | Frequency bin array |
| `ufft_low_pass(f, fc)` | Brick-wall LP response |
| `ufft_filter(x, H, dt)` | Frequency-domain filter |
