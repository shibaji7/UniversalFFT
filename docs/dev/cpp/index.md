# C++ — universalfft.hpp

Source: [`cpp/universalfft.hpp`](https://github.com/shibaji7/UniversalFFT/blob/main/cpp/universalfft.hpp)

## Design

Header-only, C++17, no external dependencies.
All functions live in the `ufft` namespace.
Uses `std::complex<double>` (`ufft::cdouble`) throughout.

## Backend

Self-contained **Cooley-Tukey radix-2 DIT FFT** (`fft_inplace`).
Raw sums in both directions; scaling applied at the wrapper level.

## Boteler Mapping

| Function | Expression |
|----------|------------|
| `fft_cc(x, dt)` | `rawFFT(x) * dt` |
| `ifft_cc(X, dt)` | `rawIFFT(X) * df`  (df = 1/(N·dt)) |
| `fft_cd(x)` | `rawFFT(x) / N` |
| `ifft_cd(X)` | `rawIFFT(X)` |

## Build

```bash
cd cpp/
make        # produces ufft_demo_cpp
make test   # runs demo and writes ../tests/data/cpp_*.csv
```

## Quick Start

```cpp
#include "universalfft.hpp"
#include <vector>
#include <complex>

using namespace ufft;

int N = 128; double dt = 1e-3;
cvec x(N);
// ... fill x with signal ...
cvec X    = fft_cc(x, dt);
cvec xrec = ifft_cc(X, dt);
```

## API

### `fft_cc(x, dt)` / `fft_cc_real(x, dt)`
CC forward — Boteler Eq. 21a.

### `ifft_cc(X, dt)`
CC inverse — Boteler Eq. 21b.

### `fft_cd(x)` / `fft_cd_real(x)`
CD forward — Boteler Eq. 22a.

### `ifft_cd(X)`
CD inverse — Boteler Eq. 22b.

### `freqs(N, dt)`
Frequency bin array in Hz, FFT-output order.

### `fft_filter(x, H, dt)`
Filter via CC pair.
