# JavaScript — universalfft.js

Source: [`js/universalfft.js`](https://github.com/shibaji7/UniversalFFT/blob/main/js/universalfft.js)

## Design

ES module (ESM), no external dependencies, works in Node.js ≥ 14 and modern browsers.
All arrays use `Float64Array` for numerical precision.

## Backend

Self-contained **Cooley-Tukey radix-2 DIT FFT** (`_fftInplace`).
Raw sums in both directions; scaling applied at the wrapper level.

## Boteler Mapping

| Function | Expression |
|----------|------------|
| `fftCC(re, im, dt)` | `rawFFT * dt` |
| `ifftCC(X, dt)` | `rawIFFT * df`  (df = 1/(N·dt)) |
| `fftCD(re, im)` | `rawFFT / N` |
| `ifftCD(X)` | `rawIFFT` |

## Installation & Usage

```bash
# No install needed — pure ES module, zero dependencies
node demo.js [output_dir]

# Tests (Node ≥ 18 built-in test runner)
node --test test.js
```

## Quick Start

```js
import { fftCCReal, ifftCC, freqs } from "./universalfft.js";

const N = 128, dt = 1e-3;
const x = new Float64Array(N).map((_, i) => Math.cos(2*Math.PI*60*i*dt));

const X    = fftCCReal(x, dt);
const xrec = ifftCC(X, dt);
```

## API

### `fftCC(re, im, dt)` / `fftCCReal(x, dt)`
CC forward — Boteler Eq. 21a.

### `ifftCC(X, dt)`
CC inverse — Boteler Eq. 21b. `X` is `{ re: Float64Array, im: Float64Array }`.

### `fftCD(re, im)` / `fftCDReal(x)`
CD forward — Boteler Eq. 22a.

### `ifftCD(X)`
CD inverse — Boteler Eq. 22b.

### `freqs(N, dt)`
Frequency bin array in Hz, FFT-output order. Returns `Float64Array`.

### `lowPassResponse(f, fc)`
Brick-wall low-pass response. Returns `{ re, im }`.

### `fftFilter(x, H, dt)`
Filter via CC pair. Returns `{ re, im }`.
