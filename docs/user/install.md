<!--
Author(s): Shibaji Chakraborty
-->

# Installation

UniversalFFT has implementations in **11 languages**. Install only the ones you need.

---

## Python

### Requirements

- Python ≥ 3.9
- NumPy ≥ 1.23

### User install (pip)

```bash
pip install -e "path/to/UniversalFFT/python/"
```

### Developer install (with test dependencies)

```bash
git clone https://github.com/shibaji7/UniversalFFT.git
cd UniversalFFT
pip install -e "python/.[dev]"
```

### Verify

```python
import universalfft
print(universalfft.__version__)   # 1.0.0
```

---

## C

### Requirements

- C99-compliant compiler (GCC ≥ 9, Clang ≥ 11, MSVC 2019+)
- `libm` (math library — standard on all platforms)
- No other external dependencies; the FFT is self-contained Cooley–Tukey

### Build

```bash
cd UniversalFFT/c
make all          # builds libuniversalfft.so + ufft_demo
make test         # runs the demo and writes CSV outputs to tests/data/
```

### Link against the shared library

```bash
gcc -O2 mycode.c -L/path/to/UniversalFFT/c -luniversalfft -lm -o myprogram
```

Or compile the source directly:

```bash
gcc -O2 mycode.c /path/to/UniversalFFT/c/universalfft.c -lm -o myprogram
```

---

## MATLAB

### Requirements

- MATLAB R2019b or later (uses complex arrays and standard `fft`/`ifft`)

### Setup

Add the `matlab/` directory to your MATLAB path once per session:

```matlab
addpath('/path/to/UniversalFFT/matlab')
```

Or permanently via **Home → Set Path → Add Folder**.

### Verify

```matlab
x  = cos(2*pi*60*(0:127)*1e-3);
Xc = ufft_cc_forward(x, 1e-3);
disp(size(Xc))   % [128 1]
```

---

## R

### Requirements

- R ≥ 4.1
- No additional packages required; uses only `stats::fft`

### Setup

Source the script at the top of your analysis:

```r
source("/path/to/UniversalFFT/r/universalfft.R")
```

### Verify

```r
x  <- cos(2*pi*60*(0:127)*1e-3)
Xc <- ufft_cc_forward(x, 1e-3)
cat(length(Xc), "\n")   # 128
```

---

## Running the tests

```bash
# Python unit tests (38 tests, 100 % coverage)
cd UniversalFFT
make test

# Generate reference vectors
make reference-vectors

# C demo
make c-test

# Cross-language validation (after running C/MATLAB/R demos)
make validate
```

!!! note "MATLAB and R demos"
    Run `ufft_demo` in MATLAB and `Rscript r/ufft_demo.R` to produce their
    respective CSV files in `tests/data/`, then call `make validate`.

---

## Julia

### Requirements

- Julia ≥ 1.6
- [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) (auto-installed via `Pkg.instantiate`)

### Setup

```bash
cd UniversalFFT/julia
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

### Verify

```julia
using UniversalFFT
x  = [2.5*cos(2π*60*n*1e-3) for n in 0:127]
Xc = fft_cc(complex.(x), 1e-3)
println(length(Xc))   # 128
```

### Tests

```bash
julia --project=. test/runtests.jl
```

---

## Fortran

### Requirements

- `gfortran` ≥ 9 (Fortran 2008 support)

### Build

```bash
cd UniversalFFT/fortran
make all          # builds ufft_demo_f90
make test         # runs demo, writes fortran_*.csv to tests/data/
```

### Use in your project

```fortran
use universalfft_mod
! compile with: gfortran -std=f2008 mycode.f90 universalfft.f90 -o myprogram
```

---

## C++

### Requirements

- C++17-capable compiler (GCC ≥ 9, Clang ≥ 7, MSVC 2019+)
- Header-only — no dependencies

### Build

```bash
cd UniversalFFT/cpp
make all          # builds ufft_demo_cpp
make test         # runs demo, writes cpp_*.csv to tests/data/
```

### Include in your project

```cpp
#include "universalfft.hpp"    // single header, copy to your project
```

---

## Rust

### Requirements

- Rust 2021 edition (rustup toolchain ≥ 1.65)
- `rustfft = "6"` (declared in `Cargo.toml` — fetched automatically)

### Build & Test

```bash
cd UniversalFFT/rust
cargo build --release
cargo test
cargo run -- ../tests/data    # writes rust_*.csv
```

---

## CUDA / HIP (GPU)

### Requirements

- **CUDA:**  NVIDIA GPU, CUDA Toolkit ≥ 11, `nvcc`, `libcufft`
- **HIP:**   AMD GPU, ROCm ≥ 5, `hipcc`, `libhipfft`

### Build

```bash
cd UniversalFFT/cuda

# NVIDIA CUDA
make

# AMD HIP
make HIP=1
```

### Test

```bash
make test        # or: make HIP=1 test
```

---

## JavaScript

### Requirements

- Node.js ≥ 14 (ES module support)
- No external npm packages

### Run

```bash
cd UniversalFFT/js
node demo.js [../tests/data]

# Tests (Node ≥ 18)
node --test test.js
```

---

## Octave

### Requirements

- GNU Octave ≥ 6 (or MATLAB R2019b+)

### Run

```bash
cd UniversalFFT/octave
octave --no-gui ufft_demo_octave.m [../tests/data]
```

---

## IDL / GDL

### Requirements

- **GDL (free):** `sudo apt install gnudatalanguage`
- **IDL (commercial):** any version with `FFT` support

### Run

```bash
gdl idl/universalfft.pro          # GDL
idl -e "@idl/universalfft.pro"    # IDL
```

!!! info "IDL normalisation note"
    IDL's `FFT(x,-1)` already includes the 1/N factor — the opposite convention
    from every other library. The `universalfft.pro` wrappers correct for this
    transparently so results match the Boteler conventions.
