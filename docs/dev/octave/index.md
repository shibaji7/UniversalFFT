# Octave — universalfft.m

Source: [`octave/universalfft.m`](https://github.com/shibaji7/UniversalFFT/blob/main/octave/universalfft.m)

## Compatibility

Compatible with both **GNU Octave** and **MATLAB** (same FFT conventions).

## Backend

Uses Octave/MATLAB's built-in `fft` / `ifft`:

| Built-in | Output |
|----------|--------|
| `fft(x)` | raw forward sum — no scaling |
| `ifft(X)` | raw sum ÷ N |

## Boteler Mapping

| Function | Expression |
|----------|------------|
| `ufft_cc_forward(x, dt)` | `fft(x) * dt` |
| `ufft_cc_inverse(X, dt)` | `ifft(X) * N * df` = `ifft(X) / dt` |
| `ufft_cd_forward(x)` | `fft(x) / N` |
| `ufft_cd_inverse(X)` | `ifft(X) * N` |

## Running the Demo

```bash
cd octave/
octave --no-gui ufft_demo_octave.m [output_dir]
```

## Quick Start

```matlab
source('universalfft.m')   % Octave
% addpath('octave/')        % MATLAB alternative

N = 128; dt = 1e-3;
x = 2.5 * cos(2*pi*60*(0:N-1)'*dt);

X    = ufft_cc_forward(x, dt);
xrec = ufft_cc_inverse(X, dt);
```

## API

### `ufft_cc_forward(x, dt)`
CC forward — Boteler Eq. 21a.

### `ufft_cc_inverse(X, dt)`
CC inverse — Boteler Eq. 21b.

### `ufft_cd_forward(x)`
CD forward — Boteler Eq. 22a.

### `ufft_cd_inverse(X)`
CD inverse — Boteler Eq. 22b.

### `ufft_freqs(N, dt)`
Frequency bin array in Hz, FFT-output order.

### `ufft_low_pass(f, fc)`
Brick-wall low-pass response.

### `ufft_filter(x, H, dt)`
Filter via CC pair.
