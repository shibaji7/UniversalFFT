# Rust — universalfft crate

Source: [`rust/src/lib.rs`](https://github.com/shibaji7/UniversalFFT/blob/main/rust/src/lib.rs)

## Backend

Uses the **`rustfft`** crate (version 6).
Both `plan_fft_forward` and `plan_fft_inverse` produce raw sums — neither divides by N.

## Boteler Mapping

| Function | Expression |
|----------|------------|
| `fft_cc(x, dt)` | `rustfft_forward(x) * dt` |
| `ifft_cc(X, dt)` | `rustfft_inverse(X) * df`  (df = 1/(N·dt)) |
| `fft_cd(x)` | `rustfft_forward(x) / N` |
| `ifft_cd(X)` | `rustfft_inverse(X)` |

## Build & Test

```bash
cd rust/
cargo build --release
cargo test
cargo run -- ../tests/data    # writes rust_*.csv
```

## Quick Start

```rust
use universalfft::*;

let n = 128usize;
let dt = 1e-3f64;
let x: Vec<C64> = /* ... */;

let X    = fft_cc(&x, dt);
let xrec = ifft_cc(&X, dt);
```

## API

### `fft_cc(x: &[C64], dt: f64) -> Vec<C64>`
CC forward — Boteler Eq. 21a.

### `ifft_cc(X: &[C64], dt: f64) -> Vec<C64>`
CC inverse — Boteler Eq. 21b.

### `fft_cd(x: &[C64]) -> Vec<C64>`
CD forward — Boteler Eq. 22a.

### `ifft_cd(X: &[C64]) -> Vec<C64>`
CD inverse — Boteler Eq. 22b.

### `freqs(n: usize, dt: f64) -> Vec<f64>`
Frequency bin array in Hz, FFT-output order.

### `low_pass_response(f: &[f64], fc: f64) -> Vec<C64>`
Brick-wall low-pass frequency response.

### `fft_filter(x: &[f64], H: &[C64], dt: f64) -> Vec<C64>`
Filter via CC pair.

### `sinc_integral_check(h: &[C64], dt: f64, tol: f64) -> bool`
Validates that `∫h(t)dt ≈ 1` (impulse response normalisation check).
