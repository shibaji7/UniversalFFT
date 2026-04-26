# Fortran — universalfft.f90

Source: [`fortran/universalfft.f90`](https://github.com/shibaji7/UniversalFFT/blob/main/fortran/universalfft.f90)

## Backend

Self-contained **Cooley-Tukey radix-2 DIT FFT** with bit-reversal permutation.
Both forward and inverse produce raw sums (no 1/N factor internally).

## Boteler Mapping

| Subroutine | Expression |
|------------|------------|
| `ufft_cc_forward(x, N, dt, X)` | `rawFFT(x) * dt` |
| `ufft_cc_inverse(X, N, dt, x)` | `rawIFFT(X) * df`  (df = 1/(N·dt)) |
| `ufft_cd_forward(x, N, X)` | `rawFFT(x) / N` |
| `ufft_cd_inverse(X, N, x)` | `rawIFFT(X)` |

## Build

```bash
cd fortran/
make        # produces ufft_demo_f90
make test   # runs demo and writes ../tests/data/fortran_*.csv
```

Requires `gfortran` with Fortran 2008 support.

## Quick Start

```fortran
use universalfft_mod
integer,  parameter :: N = 128
real(dp), parameter :: dt = 1e-3_dp
complex(dp) :: x(N), X_cc(N), x_rec(N)
! ... fill x ...
call ufft_cc_forward(x, N, dt, X_cc)
call ufft_cc_inverse(X_cc, N, dt, x_rec)
```

## API

### `ufft_cc_forward(x, N, dt, X)`
CC forward — Boteler Eq. 21a.

### `ufft_cc_inverse(X, N, dt, x)`
CC inverse — Boteler Eq. 21b.

### `ufft_cd_forward(x, N, X)`
CD forward — Boteler Eq. 22a.

### `ufft_cd_inverse(X, N, x)`
CD inverse — Boteler Eq. 22b.

### `ufft_freqs(N, dt, f)`
Frequency bin array in Hz, FFT-output order.

### `ufft_filter(x, H, N, dt, y)`
Filter via CC pair.
