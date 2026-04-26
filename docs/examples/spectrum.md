<!--
Author(s): Shibaji Chakraborty
-->

# Example 3 — Spectrum Determination (CD pair)

*Corresponds to Boteler (2012) Application 3, §4.3.*

When a time-domain signal is known to be periodic — for example a 60 Hz
power-system current distorted by a GIC-induced DC offset — the frequency
domain contains only discrete harmonics, not a continuous spectrum. The
**CD forward transform** is appropriate here.

The physical check: the CD forward of \(A\cos(2\pi f_1 t)\) produces spikes
of amplitude \(A/2\) at \(\pm f_1\). Adding them via Euler's formula recovers
the full amplitude \(A\).

!!! note "Bin alignment"
    For exact spike magnitudes, the signal frequency must fall exactly on a DFT
    bin: \(f_1 = k \cdot \Delta f\) for integer \(k\). If \(f_1\) is off-bin,
    spectral leakage reduces the spike magnitude. In this example we use
    \(f_s = 10\,\text{kHz}\), \(N = 1024\), so \(\Delta f \approx 9.77\,\text{Hz}\).

---

## Result

Waveform (top) and its CD-forward discrete spectrum (bottom). Annotations
confirm that each spike equals \(A/2\):

![Spectrum determination](../assets/figures/spectrum_example.png)

---

## Code

=== "Python"

    ```python
    import numpy as np
    from universalfft import fft_cd, freqs

    fs = 10_000.0               # 10 kHz sampling
    dt = 1.0 / fs
    N  = 1024

    t  = np.arange(N) * dt
    A1 = 1.0;  f1 = 60.0       # fundamental
    A3 = 0.3;  f3 = 180.0      # 3rd harmonic (GIC distortion)

    x = A1 * np.cos(2*np.pi*f1*t) + A3 * np.cos(2*np.pi*f3*t)

    X  = fft_cd(x)              # CD forward
    f  = freqs(N, dt)

    i1 = np.argmin(np.abs(f - f1))
    i3 = np.argmin(np.abs(f - f3))

    print(f"|X[{f1:.0f} Hz]| = {abs(X[i1]):.4f}  (expect {A1/2:.4f} = A1/2)")
    print(f"|X[{f3:.0f} Hz]| = {abs(X[i3]):.4f}  (expect {A3/2:.4f} = A3/2)")
    ```

=== "C"

    ```c
    #include "universalfft.h"
    #include <math.h>
    #include <stdio.h>

    /* N=1024, dt=1e-4 */
    for (size_t n = 0; n < N; n++) {
        x_re[n] = cos(2*M_PI*60.0*n*dt) + 0.3*cos(2*M_PI*180.0*n*dt);
        x_im[n] = 0.0;
    }
    ufft_cd_forward(x_re, x_im, X_re, X_im, N);   /* CD forward */

    /* bin index = round(f1 * N * dt) */
    size_t i1 = (size_t)(60.0 * N * dt + 0.5);
    double amp = sqrt(X_re[i1]*X_re[i1] + X_im[i1]*X_im[i1]);
    printf("|X[60 Hz]| = %.4f  (expect 0.5000)\n", amp);
    ```

=== "C++"

    ```cpp
    #include "universalfft.hpp"
    #include <cmath>
    #include <cstdio>
    using namespace ufft;

    int N = 1024; double dt = 1e-4;
    cvec x(N);
    for (int n = 0; n < N; ++n)
        x[n] = {std::cos(2*M_PI*60.0*n*dt) + 0.3*std::cos(2*M_PI*180.0*n*dt), 0.0};

    cvec X = fft_cd(x);            // CD forward
    auto f = freqs(N, dt);

    // find bin closest to 60 Hz
    int i1 = std::min_element(f.begin(), f.end(), [](double a, double b){
        return std::abs(a - 60.0) < std::abs(b - 60.0);
    }) - f.begin();

    printf("|X[60 Hz]| = %.4f  (expect 0.5000)\n", std::abs(X[i1]));
    ```

=== "Fortran"

    ```fortran
    use universalfft_mod
    use iso_fortran_env, only: dp => real64
    implicit none
    integer,  parameter :: N  = 1024
    real(dp), parameter :: dt = 1.0e-4_dp
    real(dp), parameter :: PI = 3.14159265358979323846_dp
    real(dp) :: x_r(0:N-1), x_i(0:N-1), X_r(0:N-1), X_i(0:N-1)
    real(dp) :: f(0:N-1), amp
    integer  :: n, k, i1, rc

    do n = 0, N-1
        x_r(n) = cos(2.0_dp*PI*60.0_dp*real(n,dp)*dt) &
                + 0.3_dp*cos(2.0_dp*PI*180.0_dp*real(n,dp)*dt)
        x_i(n) = 0.0_dp
    end do

    rc = ufft_cd_forward(x_r, x_i, X_r, X_i, N)     ! CD forward
    call ufft_freqs(f, N, dt)

    i1 = 0
    do k = 1, N-1
        if (abs(f(k) - 60.0_dp) < abs(f(i1) - 60.0_dp)) i1 = k
    end do

    amp = sqrt(X_r(i1)**2 + X_i(i1)**2)
    print '(A,F8.4,A)', '|X[60 Hz]| = ', amp, '  (expect 0.5000)'
    ```

=== "Julia"

    ```julia
    using UniversalFFT

    fs = 10_000.0; dt = 1/fs; N = 1024
    t  = (0:N-1) .* dt
    x  = cos.(2π*60.0*t) + 0.3*cos.(2π*180.0*t)

    X = fft_cd(complex.(x))           # CD forward
    f = freqs(N, dt)

    i1 = argmin(abs.(f .- 60.0))
    i3 = argmin(abs.(f .- 180.0))

    println("|X[60 Hz]|  = $(round(abs(X[i1]), digits=4))  (expect 0.5000)")
    println("|X[180 Hz]| = $(round(abs(X[i3]), digits=4))  (expect 0.1500)")
    ```

=== "Rust"

    ```rust
    use universalfft::*;
    use std::f64::consts::PI;

    let n: usize = 1024;
    let dt = 1e-4f64;
    let x: Vec<C64> = (0..n)
        .map(|k| {
            let t = k as f64 * dt;
            C64::new((2.0*PI*60.0*t).cos() + 0.3*(2.0*PI*180.0*t).cos(), 0.0)
        })
        .collect();

    let X = fft_cd(&x);               // CD forward
    let f = freqs(n, dt);

    let i1 = f.iter().enumerate()
        .min_by(|(_, a), (_, b)| ((*a-60.0).abs()).partial_cmp(&((*b-60.0).abs())).unwrap())
        .map(|(i, _)| i).unwrap();

    println!("|X[60 Hz]| = {:.4}  (expect 0.5000)", X[i1].norm());
    ```

=== "MATLAB"

    ```matlab
    addpath('matlab/')
    fs = 10000.0;  dt = 1/fs;  N = 1024;
    t  = (0:N-1)' * dt;
    x  = cos(2*pi*60*t) + 0.3*cos(2*pi*180*t);

    X  = ufft_cd_forward(x);          % CD forward
    f  = ufft_freqs(N, dt);

    [~, i1] = min(abs(f - 60));
    [~, i3] = min(abs(f - 180));
    fprintf('|X[60 Hz]|  = %.4f  (expect 0.5000)\n', abs(X(i1)));
    fprintf('|X[180 Hz]| = %.4f  (expect 0.1500)\n', abs(X(i3)));
    ```

=== "R"

    ```r
    source("r/universalfft.R")
    fs <- 10000;  dt <- 1/fs;  N <- 1024L
    t  <- (0:(N-1)) * dt
    x  <- cos(2*pi*60*t) + 0.3*cos(2*pi*180*t)

    X  <- ufft_cd_forward(x)          # CD forward
    f  <- ufft_freqs(N, dt)

    i1 <- which.min(abs(f - 60))
    i3 <- which.min(abs(f - 180))
    cat(sprintf("|X[60 Hz]|  = %.4f  (expect 0.5000)\n", Mod(X[i1])))
    cat(sprintf("|X[180 Hz]| = %.4f  (expect 0.1500)\n", Mod(X[i3])))
    ```

=== "JavaScript"

    ```js
    import { fftCDReal, freqs } from "./universalfft.js";

    const N = 1024, dt = 1e-4;
    const x = new Float64Array(N).map((_, k) =>
        Math.cos(2*Math.PI*60*k*dt) + 0.3*Math.cos(2*Math.PI*180*k*dt));

    const X = fftCDReal(x);           // CD forward
    const f = freqs(N, dt);

    let i1 = 0;
    for (let k = 1; k < N; k++)
        if (Math.abs(f[k] - 60) < Math.abs(f[i1] - 60)) i1 = k;

    const mag = Math.sqrt(X.re[i1]**2 + X.im[i1]**2);
    console.log(`|X[60 Hz]| = ${mag.toFixed(4)}  (expect 0.5000)`);
    ```

=== "Octave"

    ```matlab
    source('octave/universalfft.m')
    fs = 10000;  dt = 1/fs;  N = 1024;
    t  = (0:N-1)' * dt;
    x  = cos(2*pi*60*t) + 0.3*cos(2*pi*180*t);

    X  = ufft_cd_forward(x);          % CD forward
    f  = ufft_freqs(N, dt);

    [~, i1] = min(abs(f - 60));
    [~, i3] = min(abs(f - 180));
    fprintf('|X[60 Hz]|  = %.4f  (expect 0.5000)\n', abs(X(i1)));
    fprintf('|X[180 Hz]| = %.4f  (expect 0.1500)\n', abs(X(i3)));
    ```

=== "CUDA / HIP"

    ```c
    #include "universalfft.cuh"
    #include <math.h>
    #include <stdio.h>

    /* N=1024, dt=1e-4 — host arrays; device memory managed internally */
    for (int n = 0; n < N; n++) {
        x_re[n] = cos(2*M_PI*60.0*n*dt) + 0.3*cos(2*M_PI*180.0*n*dt);
        x_im[n] = 0.0;
    }
    ufft_cd_forward_host(x_re, x_im, X_re, X_im, N);   /* CD forward */

    int i1 = (int)(60.0 * N * dt + 0.5);
    double amp = sqrt(X_re[i1]*X_re[i1] + X_im[i1]*X_im[i1]);
    printf("|X[60 Hz]| = %.4f  (expect 0.5000)\n", amp);
    ```

=== "IDL / GDL"

    ```idl
    ; GDL (free): gdl idl/universalfft.pro
    ; IDL:        idl -e "@idl/universalfft.pro"
    @universalfft.pro

    N = 1024L & dt = 1D-4
    k = DINDGEN(N)
    x = COS(2D*!DPI*60D*k*dt) + 0.3D*COS(2D*!DPI*180D*k*dt)

    X = ufft_cd_forward(x)             ; CD forward
    f = ufft_freqs(N, dt)

    i1 = (MIN(ABS(f - 60D), /SUBSCRIPT_MAX, sub))[0]   ; nearest bin to 60 Hz
    amp = ABS(X[i1])
    PRINT, '|X[60 Hz]| =', amp, '  (expect 0.5000)'
    ```
