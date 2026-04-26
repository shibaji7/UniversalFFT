# Julia — UniversalFFT.jl

Source: [`julia/src/UniversalFFT.jl`](https://github.com/shibaji7/UniversalFFT/blob/main/julia/src/UniversalFFT.jl)

## Backend

Uses **FFTW.jl**, which follows the same convention as NumPy:

| Direction | FFTW.jl output |
|-----------|----------------|
| `fft(x)`  | raw forward sum — no scaling |
| `ifft(X)` | raw sum ÷ N   |

## Boteler Mapping

| Function       | Expression |
|---------------|------------|
| `fft_cc(x, dt)` | `FFTW.fft(x) * dt` |
| `ifft_cc(X, dt)` | `FFTW.ifft(X) * N * df` = `FFTW.ifft(X) / dt` |
| `fft_cd(x)` | `FFTW.fft(x) / N` |
| `ifft_cd(X)` | `FFTW.ifft(X) * N` |

## Installation

```julia
# From the julia/ directory:
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

## Quick Start

```julia
using UniversalFFT

N, dt, f1, A = 128, 1e-3, 60.0, 2.5
x  = [A * cos(2π * f1 * n * dt) for n in 0:N-1]
cx = complex.(x)

X    = fft_cc(cx, dt)
xrec = ifft_cc(X, dt)
```

## Running Tests

```bash
julia --project=. test/runtests.jl
```

## API

### `fft_cc(x, dt)`
CC forward — Boteler Eq. 21a. Returns `Vector{ComplexF64}`.

### `ifft_cc(X, dt)`
CC inverse — Boteler Eq. 21b. Returns `Vector{ComplexF64}`.

### `fft_cd(x)`
CD forward — Boteler Eq. 22a. Returns `Vector{ComplexF64}`.

### `ifft_cd(X)`
CD inverse — Boteler Eq. 22b. Returns `Vector{ComplexF64}`.

### `freqs(N, dt)`
Frequency bin array in Hz, FFT-output order. Returns `Vector{Float64}`.

### `low_pass_response(f, fc)`
Brick-wall low-pass response. Returns `Vector{ComplexF64}`.

### `fft_filter(x, H, dt)`
Filter via CC pair. Returns `Vector{ComplexF64}`.
