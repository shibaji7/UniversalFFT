# IDL / GDL — universalfft.pro

Source: [`idl/universalfft.pro`](https://github.com/shibaji7/UniversalFFT/blob/main/idl/universalfft.pro)

## GDL — Free Open-Source IDL Runtime

**IDL** (Interactive Data Language) is a proprietary product by NV5 Geospatial.
The `.pro` stubs are 100% compatible with **GDL (GNU Data Language)**, the free alternative:

```bash
# Install GDL
sudo apt install gnudatalanguage      # Debian / Ubuntu
brew install gnudatalanguage          # macOS (Homebrew)

# Run the demo
gdl idl/universalfft.pro
```

## IDL / GDL Normalisation — Important!

IDL's `FFT(x, -1)` already **includes 1/N** in the **forward** direction —
the opposite convention from NumPy, MATLAB, Julia, and every other library here.

| IDL call | Output |
|----------|--------|
| `FFT(x, -1)` | `(1/N) * rawForwardSum` |
| `FFT(X, +1)` | raw inverse sum (no 1/N) |

## Boteler Mapping

| Function | Expression |
|----------|------------|
| `ufft_cc_forward(x, dt)` | `FFT(x,-1) * N * dt`  (undo IDL's 1/N, then × dt) |
| `ufft_cc_inverse(X, dt)` | `FFT(X,+1) * df`  (df = 1/(N·dt)) |
| `ufft_cd_forward(x)` | `FFT(x,-1)`  (IDL's 1/N is exactly the CD definition) |
| `ufft_cd_inverse(X)` | `FFT(X,+1)` |

## Quick Start

```idl
; Source definitions
@universalfft.pro

N = 128L & dt = 1D-3 & f1 = 60D & A = 2.5D
t = DINDGEN(N) * dt
x = A * COS(2D * !DPI * f1 * t)

X_cc  = ufft_cc_forward(x, dt)
x_rec = ufft_cc_inverse(X_cc, dt)

PRINT, MAX(ABS(REAL_PART(x_rec) - x))   ; expect < 1e-9
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

### `universalfft_demo [, outdir=outdir]`
Run the full demo (cosine CC round-trip + impulse response).
