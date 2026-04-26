/**
 * universalfft.c — Boteler (2012) FFT conventions, self-contained C99.
 *
 * FFT algorithm: Cooley–Tukey radix-2 decimation-in-time.
 * Bit-reversal permutation then butterfly passes.
 */

#include "universalfft.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* -----------------------------------------------------------------------
 * Internal helpers
 * ----------------------------------------------------------------------- */

static int is_power_of_two(size_t n) {
    return n > 0 && (n & (n - 1)) == 0;
}

static void bit_reverse_permute(double *re, double *im, size_t N) {
    size_t j = 0;
    for (size_t i = 1; i < N; i++) {
        size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double tmp;
            tmp = re[i]; re[i] = re[j]; re[j] = tmp;
            tmp = im[i]; im[i] = im[j]; im[j] = tmp;
        }
    }
}

/* -----------------------------------------------------------------------
 * Core: in-place Cooley–Tukey radix-2 DIT FFT.
 * sign = -1 → forward (e^{-i2πkn/N}), +1 → inverse (before 1/N).
 * ----------------------------------------------------------------------- */
int ufft_inplace(double *x_re, double *x_im, size_t N, int inverse) {
    if (!is_power_of_two(N)) return -1;

    bit_reverse_permute(x_re, x_im, N);

    double sign = inverse ? +1.0 : -1.0;

    for (size_t len = 2; len <= N; len <<= 1) {
        double ang = sign * 2.0 * M_PI / (double)len;
        double wr = cos(ang);
        double wi = sin(ang);

        for (size_t i = 0; i < N; i += len) {
            double wre = 1.0, wim = 0.0;
            for (size_t j = 0; j < len / 2; j++) {
                size_t u = i + j;
                size_t v = i + j + len / 2;
                /* butterfly: t = w * x[v] */
                double tre = wre * x_re[v] - wim * x_im[v];
                double tim = wre * x_im[v] + wim * x_re[v];
                x_re[v] = x_re[u] - tre;
                x_im[v] = x_im[u] - tim;
                x_re[u] += tre;
                x_im[u] += tim;
                /* w *= (wr + i*wi) */
                double wre_new = wre * wr - wim * wi;
                wim = wre * wi + wim * wr;
                wre = wre_new;
            }
        }
    }
    return 0;
}

/* -----------------------------------------------------------------------
 * CC forward: X[k] = Σ x[n] e^{-i2πkn/N} * dt
 * ----------------------------------------------------------------------- */
int ufft_cc_forward(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N, double dt)
{
    memcpy(out_re, in_re, N * sizeof(double));
    memcpy(out_im, in_im, N * sizeof(double));
    int rc = ufft_inplace(out_re, out_im, N, 0);
    if (rc != 0) return rc;
    for (size_t k = 0; k < N; k++) {
        out_re[k] *= dt;
        out_im[k] *= dt;
    }
    return 0;
}

/* -----------------------------------------------------------------------
 * CC inverse: x[n] = Σ X[k] e^{+i2πkn/N} * Δf   (Δf = 1/(N*dt))
 * NumPy-ifft applies 1/N; we need Δf = 1/(N*dt), so factor = N * Δf = 1/dt.
 * ----------------------------------------------------------------------- */
int ufft_cc_inverse(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N, double dt)
{
    memcpy(out_re, in_re, N * sizeof(double));
    memcpy(out_im, in_im, N * sizeof(double));
    int rc = ufft_inplace(out_re, out_im, N, 1);   /* inverse twiddles */
    if (rc != 0) return rc;
    /* ufft_inplace does NOT apply 1/N — we scale by Δf = 1/(N*dt) */
    double df = 1.0 / ((double)N * dt);
    for (size_t n = 0; n < N; n++) {
        out_re[n] *= df;
        out_im[n] *= df;
    }
    return 0;
}

/* -----------------------------------------------------------------------
 * CD forward: X[k] = (1/N) Σ x[n] e^{-i2πkn/N}
 * ----------------------------------------------------------------------- */
int ufft_cd_forward(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N)
{
    memcpy(out_re, in_re, N * sizeof(double));
    memcpy(out_im, in_im, N * sizeof(double));
    int rc = ufft_inplace(out_re, out_im, N, 0);
    if (rc != 0) return rc;
    double inv_N = 1.0 / (double)N;
    for (size_t k = 0; k < N; k++) {
        out_re[k] *= inv_N;
        out_im[k] *= inv_N;
    }
    return 0;
}

/* -----------------------------------------------------------------------
 * CD inverse: x[n] = Σ X[k] e^{+i2πkn/N}  (raw summation — no 1/N)
 * ----------------------------------------------------------------------- */
int ufft_cd_inverse(
    const double *in_re, const double *in_im,
    double *out_re, double *out_im,
    size_t N)
{
    memcpy(out_re, in_re, N * sizeof(double));
    memcpy(out_im, in_im, N * sizeof(double));
    /* forward twiddles with conjugate sign = inverse, no 1/N */
    return ufft_inplace(out_re, out_im, N, 1);
}

/* -----------------------------------------------------------------------
 * Frequency bin array
 * ----------------------------------------------------------------------- */
void ufft_freqs(double *f, size_t N, double dt) {
    double df = 1.0 / ((double)N * dt);
    for (size_t k = 0; k < N; k++) {
        if (k < N / 2)
            f[k] = (double)k * df;
        else
            f[k] = ((double)k - (double)N) * df;
    }
}

/* -----------------------------------------------------------------------
 * Filter via CC pair
 * ----------------------------------------------------------------------- */
int ufft_filter(
    const double *in_re, const double *in_im,
    const double *H_re,  const double *H_im,
    double *out_re, double *out_im,
    size_t N, double dt)
{
    double *X_re = malloc(N * sizeof(double));
    double *X_im = malloc(N * sizeof(double));
    if (!X_re || !X_im) { free(X_re); free(X_im); return -2; }

    int rc = ufft_cc_forward(in_re, in_im, X_re, X_im, N, dt);
    if (rc != 0) { free(X_re); free(X_im); return rc; }

    /* complex multiply: Y = X * H */
    for (size_t k = 0; k < N; k++) {
        double yr = X_re[k] * H_re[k] - X_im[k] * H_im[k];
        double yi = X_re[k] * H_im[k] + X_im[k] * H_re[k];
        X_re[k] = yr;
        X_im[k] = yi;
    }

    rc = ufft_cc_inverse(X_re, X_im, out_re, out_im, N, dt);
    free(X_re);
    free(X_im);
    return rc;
}
