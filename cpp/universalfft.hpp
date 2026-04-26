/**
 * universalfft.hpp — Boteler (2012) compliant FFT/IFFT, C++17 header-only.
 *
 * Uses std::complex<double> and a self-contained Cooley–Tukey radix-2 DIT
 * algorithm — no external dependencies.  N must be a power of two.
 *
 * Conventions (identical to the C and Fortran implementations):
 *   CC forward : X[k] = Σ x[n] e^{-i2πkn/N} Δt
 *   CC inverse : x[n] = Σ X[k] e^{+i2πkn/N} Δf   (Δf = 1/(N Δt))
 *   CD forward : X[k] = (1/N) Σ x[n] e^{-i2πkn/N}
 *   CD inverse : x[n] = Σ X[k] e^{+i2πkn/N}        (raw summation)
 */
#pragma once

#include <complex>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>

namespace ufft {

using cdouble = std::complex<double>;
using cvec    = std::vector<cdouble>;
using rvec    = std::vector<double>;

static constexpr double PI = 3.14159265358979323846;

// ── Helpers ──────────────────────────────────────────────────────────────────

inline bool is_power_of_two(std::size_t n) {
    return n > 0 && (n & (n - 1)) == 0;
}

inline void bit_reverse(cvec& x) {
    std::size_t N = x.size(), j = 0;
    for (std::size_t i = 1; i < N; ++i) {
        std::size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(x[i], x[j]);
    }
}

// ── Core in-place Cooley–Tukey ────────────────────────────────────────────────
/**
 * In-place radix-2 DIT FFT.
 * @param x       Input/output complex vector (length must be power of two).
 * @param inverse false → forward (e^{-i2πkn/N}), true → backward (e^{+i2πkn/N}).
 *                Note: no 1/N scaling is applied in either direction.
 */
inline void fft_inplace(cvec& x, bool inverse) {
    std::size_t N = x.size();
    if (!is_power_of_two(N))
        throw std::invalid_argument("ufft: N must be a power of two");
    bit_reverse(x);
    double sign = inverse ? +1.0 : -1.0;
    for (std::size_t len = 2; len <= N; len <<= 1) {
        cdouble w_step = std::exp(cdouble(0.0, sign * 2.0 * PI / static_cast<double>(len)));
        for (std::size_t i = 0; i < N; i += len) {
            cdouble w(1.0, 0.0);
            for (std::size_t j = 0; j < len / 2; ++j) {
                cdouble u = x[i + j];
                cdouble t = w * x[i + j + len / 2];
                x[i + j]           = u + t;
                x[i + j + len / 2] = u - t;
                w *= w_step;
            }
        }
    }
}

// ── CC forward  (Boteler Eq. 21a) ────────────────────────────────────────────
/**
 * X[k] = Σ x[n] e^{-i2πkn/N} Δt
 */
inline cvec fft_cc(cvec x, double dt) {
    fft_inplace(x, false);
    for (auto& v : x) v *= dt;
    return x;
}

// Real overload
inline cvec fft_cc(const rvec& x, double dt) {
    cvec cx(x.begin(), x.end());
    return fft_cc(std::move(cx), dt);
}

// ── CC inverse  (Boteler Eq. 21b) ────────────────────────────────────────────
/**
 * x[n] = Σ X[k] e^{+i2πkn/N} Δf,   Δf = 1/(N Δt)
 *
 * fft_inplace(inverse=true) gives the raw summation (no 1/N).
 * Multiply by Δf to complete the approximation of the Fourier integral.
 */
inline cvec ifft_cc(cvec X, double dt) {
    std::size_t N = X.size();
    double df = 1.0 / (static_cast<double>(N) * dt);
    fft_inplace(X, true);          // raw summation
    for (auto& v : X) v *= df;
    return X;
}

// ── CD forward  (Boteler Eq. 22a) ────────────────────────────────────────────
/** X[k] = (1/N) Σ x[n] e^{-i2πkn/N} */
inline cvec fft_cd(cvec x) {
    std::size_t N = x.size();
    fft_inplace(x, false);
    double inv_N = 1.0 / static_cast<double>(N);
    for (auto& v : x) v *= inv_N;
    return x;
}

inline cvec fft_cd(const rvec& x) {
    cvec cx(x.begin(), x.end());
    return fft_cd(std::move(cx));
}

// ── CD inverse  (Boteler Eq. 22b) ────────────────────────────────────────────
/** x[n] = Σ X[k] e^{+i2πkn/N}  (raw summation) */
inline cvec ifft_cd(cvec X) {
    fft_inplace(X, true);
    return X;
}

// ── Frequency bin array ──────────────────────────────────────────────────────
/** FFT frequency bins in Hz, FFT-output order (matches numpy.fft.fftfreq). */
inline rvec freqs(std::size_t N, double dt) {
    double df = 1.0 / (static_cast<double>(N) * dt);
    rvec f(N);
    for (std::size_t k = 0; k < N; ++k) {
        if (k < N / 2)
            f[k] = static_cast<double>(k) * df;
        else
            f[k] = (static_cast<double>(k) - static_cast<double>(N)) * df;
    }
    return f;
}

// ── Low-pass brick-wall response ─────────────────────────────────────────────
inline cvec low_pass_response(const rvec& f, double fc) {
    cvec H(f.size());
    for (std::size_t k = 0; k < f.size(); ++k)
        H[k] = (std::abs(f[k]) <= fc) ? cdouble(1.0, 0.0) : cdouble(0.0, 0.0);
    return H;
}

// ── Filter via CC pair ────────────────────────────────────────────────────────
inline cvec fft_filter(const rvec& x, const cvec& H, double dt) {
    if (x.size() != H.size())
        throw std::invalid_argument("ufft::fft_filter: x and H must have the same size");
    cvec X = fft_cc(x, dt);
    for (std::size_t k = 0; k < X.size(); ++k) X[k] *= H[k];
    return ifft_cc(std::move(X), dt);
}

// ── Parseval / impulse check ─────────────────────────────────────────────────
inline bool sinc_integral_check(const cvec& h, double dt, double tol = 1e-4) {
    double integral = 0.0;
    for (const auto& v : h) integral += v.real();
    integral *= dt;
    return std::abs(integral - 1.0) < tol;
}

} // namespace ufft
