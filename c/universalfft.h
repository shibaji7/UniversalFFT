/**
 * universalfft.h — Boteler (2012) compliant FFT/IFFT wrappers in C.
 *
 * Implements the same two transform pairs as the Python package:
 *
 *   CC (Continuous–Continuous, Boteler Eq. 21):
 *     Forward:  X[k] = Σ x[n] e^{-i2πkn/N} Δt
 *     Inverse:  x[n] = Σ X[k] e^{+i2πkn/N} Δf
 *
 *   CD (Continuous–Discrete, Boteler Eq. 22):
 *     Forward:  X[k] = (1/N) Σ x[n] e^{-i2πkn/N}
 *     Inverse:  x[n] = Σ X[k] e^{+i2πkn/N}
 *
 * Implementation uses a self-contained Cooley–Tukey radix-2 DIT FFT so
 * the library has no external dependencies.  N must be a power of two.
 *
 * Build:
 *   cc -O2 -o libuniversalfft.so -shared -fPIC universalfft.c -lm
 *   cc -O2 -o ufft_demo demo.c universalfft.c -lm
 */

#ifndef UNIVERSALFFT_H
#define UNIVERSALFFT_H

#include <stddef.h>
#include <complex.h>   /* C99 complex numbers */

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------
 * Core in-place DFT (Cooley–Tukey radix-2 DIT).
 * x_re / x_im : real and imaginary parts, length N (modified in place).
 * inverse      : 0 → forward (e^{-i2πkn/N}), 1 → inverse (e^{+i2πkn/N}).
 * Returns 0 on success, -1 if N is not a power of two.
 * ----------------------------------------------------------------------- */
int ufft_inplace(double *x_re, double *x_im, size_t N, int inverse);

/* -----------------------------------------------------------------------
 * CC forward  (Boteler Eq. 21a)
 *   X[k] = Σ x[n] e^{-i2πkn/N} Δt
 *
 * in_re / in_im  : input real + imag (imag may be all zeros for real input)
 * out_re / out_im: output arrays (length N)
 * dt             : sampling interval (seconds)
 * ----------------------------------------------------------------------- */
int ufft_cc_forward(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N, double dt);

/* -----------------------------------------------------------------------
 * CC inverse  (Boteler Eq. 21b)
 *   x[n] = Σ X[k] e^{+i2πkn/N} Δf     where Δf = 1/(N·Δt)
 * ----------------------------------------------------------------------- */
int ufft_cc_inverse(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N, double dt);

/* -----------------------------------------------------------------------
 * CD forward  (Boteler Eq. 22a)
 *   X[k] = (1/N) Σ x[n] e^{-i2πkn/N}
 * ----------------------------------------------------------------------- */
int ufft_cd_forward(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N);

/* -----------------------------------------------------------------------
 * CD inverse  (Boteler Eq. 22b)
 *   x[n] = Σ X[k] e^{+i2πkn/N}   (raw summation, no 1/N)
 * ----------------------------------------------------------------------- */
int ufft_cd_inverse(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N);

/* -----------------------------------------------------------------------
 * Frequency bin array (Hz), length N, FFT output order.
 *   f[k] = k/(N·dt)  for k < N/2
 *   f[k] = (k-N)/(N·dt)  for k >= N/2
 * ----------------------------------------------------------------------- */
void ufft_freqs(double *f, size_t N, double dt);

/* -----------------------------------------------------------------------
 * Apply a frequency-domain filter H in-place using CC pair.
 *   in_re / in_im : time-domain input (length N)
 *   H_re / H_im   : transfer function at FFT frequency bins (length N)
 *   out_re / out_im: filtered time-domain output (length N)
 * ----------------------------------------------------------------------- */
int ufft_filter(
    const double *in_re, const double *in_im,
    const double *H_re,  const double *H_im,
    double *out_re, double *out_im,
    size_t N, double dt);

#ifdef __cplusplus
}
#endif

#endif /* UNIVERSALFFT_H */
