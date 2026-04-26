<!--
Author(s): Shibaji Chakraborty
-->

# R Functions — `universalfft.R`

Source the single file to load all functions:

```r
source("/path/to/UniversalFFT/r/universalfft.R")
```

!!! warning "R's `fft(inverse=TRUE)` does NOT divide by N"
    Unlike NumPy and MATLAB, R's `stats::fft(inverse=TRUE)` returns the raw
    summation \(\sum X[k] e^{+i2\pi kn/N}\) with no \(1/N\) factor. This
    makes the R CD inverse a direct call — but it also means the CC inverse
    formula differs from the other languages (multiply by \(\Delta f\) alone,
    not by \(N \cdot \Delta f\)).

---

## `ufft_cc_forward(x, dt)`

\[X[k] = \texttt{fft}(x) \cdot \Delta t\]

---

## `ufft_cc_inverse(X, dt)`

\[x[n] = \texttt{fft}(X, \texttt{inverse=TRUE}) \cdot \Delta f\]

where \(\Delta f = 1/(N \Delta t)\). Because R's `fft(inverse=TRUE)` already
returns the raw sum (no \(1/N\)), multiplying by \(\Delta f\) is sufficient.

---

## `ufft_cd_forward(x, dt = NULL)`

\[X[k] = \texttt{fft}(x) / N\]

---

## `ufft_cd_inverse(X, dt = NULL)`

```r
stats::fft(X, inverse = TRUE)
```

R's `fft(inverse=TRUE)` is already the raw summation — no extra factor needed.

---

## `ufft_freqs(N, dt)`

Returns a numeric vector of length N with the FFT frequency bins in Hz
(FFT output order, matching `numpy.fft.fftfreq`).

---

## `ufft_filter(x, H, dt)`

```r
ufft_cc_inverse(ufft_cc_forward(x, dt) * H, dt)
```

---

## `ufft_lowpass(f, fc)`

Convenience: returns `as.numeric(abs(f) <= fc)` (0.0 or 1.0).
