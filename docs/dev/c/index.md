<!--
Author(s): Shibaji Chakraborty
-->

# C Library — `universalfft.h` / `universalfft.c`

A self-contained C99 implementation with no external dependencies. The FFT
algorithm is Cooley–Tukey radix-2 decimation-in-time (DIT); N must be a
power of two.

## Build

```bash
cd c/
make lib      # → libuniversalfft.so
make demo     # → ufft_demo
```

## Core function — `ufft_inplace`

```c
int ufft_inplace(double *x_re, double *x_im, size_t N, int inverse);
```

In-place Cooley–Tukey FFT. Sets the exponential sign:
- `inverse = 0` → forward (\(e^{-i2\pi kn/N}\))
- `inverse = 1` → inverse (\(e^{+i2\pi kn/N}\)) — **no** 1/N is applied

Returns `0` on success, `-1` if N is not a power of two.

---

## CC forward — `ufft_cc_forward`

```c
int ufft_cc_forward(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N, double dt);
```

\[X[k] = \sum_{n=0}^{N-1} x[n]\, e^{-i2\pi kn/N}\, \Delta t\]

---

## CC inverse — `ufft_cc_inverse`

```c
int ufft_cc_inverse(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N, double dt);
```

\[x[n] = \sum_{k=0}^{N-1} X[k]\, e^{+i2\pi kn/N}\, \Delta f, \quad \Delta f = \tfrac{1}{N\Delta t}\]

!!! note "Scaling"
    `ufft_inplace` with `inverse=1` returns the raw summation (no \(1/N\)).
    `ufft_cc_inverse` multiplies by \(\Delta f = 1/(N\Delta t)\) directly —
    no intermediate \(1/N\) is applied and removed.

---

## CD forward — `ufft_cd_forward`

```c
int ufft_cd_forward(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N);
```

\[X[k] = \frac{1}{N}\sum_{n=0}^{N-1} x[n]\, e^{-i2\pi kn/N}\]

---

## CD inverse — `ufft_cd_inverse`

```c
int ufft_cd_inverse(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N);
```

Raw summation: \(x[n] = \sum_{k=0}^{N-1} X[k]\, e^{+i2\pi kn/N}\)

---

## Frequency bins — `ufft_freqs`

```c
void ufft_freqs(double *f, size_t N, double dt);
```

Fills `f` with the FFT frequency bin array in Hz (FFT output order, matching
`numpy.fft.fftfreq`).

---

## Filter — `ufft_filter`

```c
int ufft_filter(
    const double *in_re, const double *in_im,
    const double *H_re,  const double *H_im,
    double *out_re, double *out_im,
    size_t N, double dt);
```

Applies transfer function H via the CC pair. Allocates two temporary arrays
of size N internally; returns `-2` on allocation failure.
