<!--
Author(s): Shibaji Chakraborty
-->

# MATLAB Functions

Six `.m` files in `matlab/` wrap MATLAB's native `fft`/`ifft` to implement the
Boteler (2012) conventions.

!!! note "MATLAB's `ifft` applies 1/N"
    MATLAB's `ifft(X)` returns `(1/N) * sum(X .* exp(+i2pi*(k-1)*(n-1)/N))`.
    All inverse functions in this library compensate for that factor.

---

## `ufft_cc_forward(x, dt)`

**File:** `matlab/ufft_cc_forward.m`

\[X[k] = \texttt{fft}(x) \cdot \Delta t\]

| Argument | Type | Description |
|----------|------|-------------|
| `x` | `double` vector, length N | Time-domain samples |
| `dt` | `double` scalar | Sampling interval (s) |
| **returns** | complex vector, length N | CC forward DFT |

---

## `ufft_cc_inverse(X, dt)`

**File:** `matlab/ufft_cc_inverse.m`

\[x[n] = \texttt{ifft}(X) \cdot N \cdot \Delta f \quad (= \texttt{ifft}(X) / \Delta t)\]

| Argument | Type | Description |
|----------|------|-------------|
| `X` | complex vector, length N | Frequency-domain values |
| `dt` | `double` scalar | Sampling interval of original time series (s) |
| **returns** | complex vector, length N | Reconstructed time-domain samples |

---

## `ufft_cd_forward(x)`

**File:** `matlab/ufft_cd_forward.m`

\[X[k] = \texttt{fft}(x) / N\]

`dt` accepted but unused (API symmetry).

---

## `ufft_cd_inverse(X)`

**File:** `matlab/ufft_cd_inverse.m`

\[x[n] = N \cdot \texttt{ifft}(X)\]

---

## `ufft_freqs(N, dt)`

**File:** `matlab/ufft_freqs.m`

Returns the FFT frequency bin array in Hz (FFT output order), matching
`numpy.fft.fftfreq(N, dt)`.

---

## `ufft_filter(x, H, dt)`

**File:** `matlab/ufft_filter.m`

Applies transfer function `H` using the CC pair:

```matlab
y = ufft_cc_inverse(ufft_cc_forward(x, dt) .* H(:), dt);
```
