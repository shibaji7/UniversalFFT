//! UniversalFFT — Boteler (2012) compliant FFT/IFFT wrappers for Rust.
//!
//! Uses the `rustfft` crate. Convention notes:
//!
//! `rustfft` produces raw summations in both directions — no 1/N factor:
//!   - Forward: Σ x[n] e^{-i2πkn/N}
//!   - Inverse: Σ X[k] e^{+i2πkn/N}
//!
//! Boteler mapping:
//!   CC forward  = rustfft_forward * dt
//!   CC inverse  = rustfft_inverse * df          (Δf = 1/(N·Δt))
//!   CD forward  = rustfft_forward / N
//!   CD inverse  = rustfft_inverse               (raw sum — already correct)

use rustfft::{FftPlanner, num_complex::Complex};
use std::sync::Arc;

pub type C64 = Complex<f64>;

// ── Internal helper: run rustfft forward or inverse ─────────────────────────
fn run_fft(x: &[C64], inverse: bool) -> Vec<C64> {
    let n = x.len();
    let mut planner = FftPlanner::<f64>::new();
    let fft = if inverse {
        planner.plan_fft_inverse(n)
    } else {
        planner.plan_fft_forward(n)
    };
    let mut buf = x.to_vec();
    fft.process(&mut buf);
    buf
}

// ── CC forward  (Boteler Eq. 21a) ────────────────────────────────────────────
/// X[k] = Σ x[n] e^{-i2πkn/N} Δt
pub fn fft_cc(x: &[C64], dt: f64) -> Vec<C64> {
    run_fft(x, false)
        .into_iter()
        .map(|v| v * dt)
        .collect()
}

/// Real input overload.
pub fn fft_cc_real(x: &[f64], dt: f64) -> Vec<C64> {
    let cx: Vec<C64> = x.iter().map(|&v| C64::new(v, 0.0)).collect();
    fft_cc(&cx, dt)
}

// ── CC inverse  (Boteler Eq. 21b) ────────────────────────────────────────────
/// x[n] = Σ X[k] e^{+i2πkn/N} Δf,   Δf = 1/(N·Δt)
pub fn ifft_cc(X: &[C64], dt: f64) -> Vec<C64> {
    let n = X.len();
    let df = 1.0 / (n as f64 * dt);
    run_fft(X, true)
        .into_iter()
        .map(|v| v * df)
        .collect()
}

// ── CD forward  (Boteler Eq. 22a) ────────────────────────────────────────────
/// X[k] = (1/N) Σ x[n] e^{-i2πkn/N}
pub fn fft_cd(x: &[C64]) -> Vec<C64> {
    let n = x.len();
    run_fft(x, false)
        .into_iter()
        .map(|v| v / n as f64)
        .collect()
}

pub fn fft_cd_real(x: &[f64]) -> Vec<C64> {
    let cx: Vec<C64> = x.iter().map(|&v| C64::new(v, 0.0)).collect();
    fft_cd(&cx)
}

// ── CD inverse  (Boteler Eq. 22b) ────────────────────────────────────────────
/// x[n] = Σ X[k] e^{+i2πkn/N}  (raw summation — rustfft inverse is already this)
pub fn ifft_cd(X: &[C64]) -> Vec<C64> {
    run_fft(X, true)
}

// ── Frequency bin array ──────────────────────────────────────────────────────
/// FFT frequency bins in Hz, FFT-output order (matches numpy.fft.fftfreq).
pub fn freqs(n: usize, dt: f64) -> Vec<f64> {
    let df = 1.0 / (n as f64 * dt);
    (0..n)
        .map(|k| {
            if k < n / 2 { k as f64 * df }
            else          { (k as f64 - n as f64) * df }
        })
        .collect()
}

// ── Low-pass brick-wall response ─────────────────────────────────────────────
pub fn low_pass_response(f: &[f64], fc: f64) -> Vec<C64> {
    f.iter()
        .map(|&fk| if fk.abs() <= fc { C64::new(1.0, 0.0) } else { C64::new(0.0, 0.0) })
        .collect()
}

// ── Filter via CC pair ────────────────────────────────────────────────────────
pub fn fft_filter(x: &[f64], H: &[C64], dt: f64) -> Vec<C64> {
    assert_eq!(x.len(), H.len(), "x and H must have the same length");
    let X: Vec<C64> = fft_cc_real(x, dt)
        .into_iter()
        .zip(H.iter())
        .map(|(xk, &hk)| xk * hk)
        .collect();
    ifft_cc(&X, dt)
}

// ── Impulse integral check ────────────────────────────────────────────────────
pub fn sinc_integral_check(h: &[C64], dt: f64, tol: f64) -> bool {
    let integral: f64 = h.iter().map(|v| v.re).sum::<f64>() * dt;
    (integral - 1.0).abs() < tol
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    const TOL: f64 = 1e-9;
    const N: usize = 128;
    const DT: f64  = 1e-3;
    const F1: f64  = 60.0;
    const A: f64   = 2.5;

    fn cosine_signal() -> Vec<f64> {
        (0..N)
            .map(|n| A * (2.0 * std::f64::consts::PI * F1 * n as f64 * DT).cos())
            .collect()
    }

    #[test]
    fn cc_roundtrip() {
        let x  = cosine_signal();
        let cx: Vec<C64> = x.iter().map(|&v| C64::new(v, 0.0)).collect();
        let X  = fft_cc(&cx, DT);
        let xr = ifft_cc(&X, DT);
        for (a, b) in x.iter().zip(xr.iter()) {
            assert_abs_diff_eq!(b.re, *a, epsilon = TOL);
        }
    }

    #[test]
    fn cd_roundtrip() {
        let x  = cosine_signal();
        let cx: Vec<C64> = x.iter().map(|&v| C64::new(v, 0.0)).collect();
        let X  = fft_cd(&cx);
        let xr = ifft_cd(&X);
        for (a, b) in x.iter().zip(xr.iter()) {
            assert_abs_diff_eq!(b.re, *a, epsilon = TOL);
        }
    }

    #[test]
    fn cd_spectrum_spike() {
        // Use f1=62.5 Hz so k = f1*N*DT = 62.5*128*1e-3 = 8 (exact bin, no leakage).
        let f1_exact = 62.5f64;
        let x: Vec<f64> = (0..N)
            .map(|n| A * (2.0 * std::f64::consts::PI * f1_exact * n as f64 * DT).cos())
            .collect();
        let cx: Vec<C64> = x.iter().map(|&v| C64::new(v, 0.0)).collect();
        let X = fft_cd(&cx);
        let f = freqs(N, DT);
        let i1 = f.iter().enumerate()
            .min_by(|(_, a), (_, b)| ((*a - f1_exact).abs()).partial_cmp(&((*b - f1_exact).abs())).unwrap())
            .map(|(i, _)| i).unwrap();
        assert_abs_diff_eq!(X[i1].norm(), A / 2.0, epsilon = 1e-9);
    }

    #[test]
    fn impulse_response_integrates_to_one() {
        let n2 = 2048usize;
        let dt2 = 60.0f64;
        let f2  = freqs(n2, dt2);
        let H   = low_pass_response(&f2, 1.0 / 3600.0);
        let h   = ifft_cc(&H, dt2);
        assert!(sinc_integral_check(&h, dt2, 1e-4));
    }

    #[test]
    fn filter_cc_equals_cd() {
        let x: Vec<f64> = (0..2048).map(|i| (i as f64 * 0.01).sin()).collect();
        let dt2 = 60.0f64;
        let f2  = freqs(2048, dt2);
        let H   = low_pass_response(&f2, 1.0 / 3600.0);
        let cx: Vec<C64> = x.iter().map(|&v| C64::new(v, 0.0)).collect();
        let y_cc = {
            let X: Vec<C64> = fft_cc(&cx, dt2).into_iter().zip(H.iter()).map(|(a, &b)| a * b).collect();
            ifft_cc(&X, dt2)
        };
        let y_cd = {
            let X: Vec<C64> = fft_cd(&cx).into_iter().zip(H.iter()).map(|(a, &b)| a * b).collect();
            ifft_cd(&X)
        };
        for (a, b) in y_cc.iter().zip(y_cd.iter()) {
            assert_abs_diff_eq!(a.re, b.re, epsilon = TOL);
        }
    }
}
