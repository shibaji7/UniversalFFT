# UniversalFFT: A Multi-Language Software Library for Convention-Consistent Fourier Transforms in Geoscience

**Shibaji Chakraborty**

*Affiliation TBD*

---

## Code Metadata

| Field | Value |
|---|---|
| Current code version | 1.0.0 |
| Permanent link to code / repository | https://github.com/shibaji7/UniversalFFT |
| Legal Software License | MIT |
| Computing platforms / Operating Systems | Linux, macOS, Windows |
| Programming languages | Python, C, C++, Fortran, Julia, Rust, MATLAB, R, JavaScript, GNU Octave, CUDA/HIP |
| Dependencies | Python ≥ 3.9, NumPy ≥ 1.23; all other backends are self-contained or use host-language FFT primitives |
| Link to developer documentation | https://universalfft.readthedocs.io |
| Software code repository | https://github.com/shibaji7/UniversalFFT |

---

## Abstract

Every numerical FFT library makes independent, silent choices about normalisation scaling. NumPy divides the inverse transform by N; R's `fft(inverse=TRUE)` does not; IDL's `FFT()` divides the forward transform by N. For a filter applied as a forward–inverse pair these factors cancel, and practitioners never notice the discrepancy. For a *single* transform — recovering a filter's causal impulse response or computing the discrete spectrum of a periodic waveform — they do not cancel, and the result is silently mis-scaled by a factor of N·Δt or N·Δf. Boteler (2012) formalised the physically correct convention for both cases in two transform pairs, CC (continuous–continuous) and CD (continuous–discrete), but provided no implementation. UniversalFFT closes that gap: it delivers a uniform six-function API across 11 scientific computing environments, corrects each backend's native normalisation internally, and validates cross-language numerical agreement to better than 10⁻⁹. All code is open-source under the MIT licence with 100% test coverage.

---

## 1. Introduction

The discrete Fourier transform (DFT) occupies a central role in geoscience signal processing — magnetotelluric analysis, geomagnetic induction modelling, and geomagnetically induced current (GIC) calculations all depend on it. Yet the DFT is not a single algorithm: it is a family of conventions, differentiated by which direction carries a 1/N factor (or Δt, or Δf), and which sign appears in the complex exponential. Any specific library implementation commits silently to one choice from this family.

Table 1 illustrates the diversity among widely used libraries. For a forward transform followed immediately by an inverse, the normalisation factors multiply to 1 regardless of which convention is used, so roundtrip reconstruction succeeds under any consistent convention. The problem emerges for single transforms. The inverse DFT of a brick-wall low-pass transfer function should recover the sinc impulse response, whose integral over time equals 1 — a consequence of the filter passing DC with unit gain. With NumPy's convention this integral is approximately 1/Δt; with R's it is N/Δt. Neither is 1. The shape of the impulse response is correct in every case, but the physical amplitude is wrong, and the error is proportional to the number of samples — the kind of bug that is easy to overlook and hard to attribute to a normalisation mismatch.

**Table 1.** Native normalisation conventions of common FFT libraries.

| Library | Forward scale | Inverse scale |
|---|:---:|:---:|
| NumPy `fft` / `ifft` | 1 | 1/N |
| MATLAB `fft` / `ifft` | 1 | 1/N |
| R `stats::fft(inverse=TRUE)` | 1 | 1 (raw sum) |
| FFTW `FFTW_FORWARD / BACKWARD` | 1 | 1 |
| IDL `FFT(x, -1)` / `FFT(X, +1)` | 1/N | 1 |

Boteler (2012) addressed this problem rigorously. Starting from the Fourier integral written in frequency *f* (not angular frequency ω, which would introduce a 2π factor in Parseval's theorem), he derived two distinct DFT approximations appropriate for geoscience use:

- **CC pair (Eqs. 21a/21b)** — both the time and frequency domains are treated as samples of continuous functions. The forward transform scales by Δt, so that the DFT sum approximates the Fourier integral F(f) = ∫ x(t) e^{−i2πft} dt. The inverse scales by Δf = 1/(NΔt). This pair is the correct choice when recovering an impulse response from a transfer function, because the sinc integral ∫ h(t) dt = 1 holds exactly under CC convention.

- **CD pair (Eqs. 22a/22b)** — the time domain is continuous but the frequency domain contains only discrete harmonics (e.g. the fundamental and overtones of a power-system waveform). The forward transform divides by N, yielding dimensionless Fourier-series coefficients: a cosine of amplitude A produces spikes of amplitude A/2 at ±f₁, which Euler's formula recovers as A. The inverse is a raw summation with no scaling.

Boteler (2012) is a *theory* paper. It stops at the equations and explicitly notes that the appropriate software implementation depends on the language and library in use. In the fourteen years since, practitioners working in different environments — Python, MATLAB, Fortran, Julia, R, and others — have been left to derive and verify the per-language corrections independently. UniversalFFT provides the missing implementation layer.

---

## 2. Software Description

### 2.1 Mathematical Foundation

UniversalFFT exposes two transform pairs matching Boteler (2012) Equations 21 and 22. Let x[n] denote N time-domain samples at spacing Δt seconds, X[k] denote their spectrum, and Δf = 1/(NΔt) denote the frequency resolution.

**CC forward (Boteler Eq. 21a)**

$$X_{\text{CC}}[k] = \sum_{n=0}^{N-1} x[n]\, e^{-i2\pi kn/N}\, \Delta t$$

**CC inverse (Boteler Eq. 21b)**

$$x_{\text{CC}}[n] = \sum_{k=0}^{N-1} X[k]\, e^{+i2\pi kn/N}\, \Delta f$$

**CD forward (Boteler Eq. 22a)**

$$X_{\text{CD}}[k] = \frac{1}{N}\sum_{n=0}^{N-1} x[n]\, e^{-i2\pi kn/N}$$

**CD inverse (Boteler Eq. 22b)**

$$x_{\text{CD}}[n] = \sum_{k=0}^{N-1} X[k]\, e^{+i2\pi kn/N}$$

The roundtrip identity Δt · Δf = 1/N guarantees that applying the CC forward followed by the CC inverse (or CD forward followed by CD inverse) reconstructs the original signal exactly. For filter applications — where a forward transform, frequency-domain multiplication, and inverse transform are chained — either convention yields the identical result because the scaling factors cancel.

**Convention selection guide:**

| Application | Correct pair | Reason |
|---|:---:|---|
| Filter (FFT → multiply → IFFT) | Either | Δt·Δf = 1/N cancels across the pair |
| Impulse response from transfer function | CC inverse | Approximates Fourier integral; sinc integral = 1 |
| Discrete spectrum of periodic waveform | CD forward | Fourier-series coefficients; spike amplitude = A/2 |

### 2.2 API Design

UniversalFFT presents six public functions in every supported language:

| Function | Role |
|---|---|
| `fft_cc(x, dt)` | CC forward transform |
| `ifft_cc(X, dt)` | CC inverse transform |
| `fft_cd(x)` | CD forward transform |
| `ifft_cd(X)` | CD inverse transform |
| `freqs(N, dt)` | FFT frequency bin array (Hz) |
| `fft_filter(x, H, dt)` | Convenience: forward → multiply → inverse |

Function names and argument order are identical across all 11 languages. A script written in Julia can be translated line-by-line to Python (or Fortran, or Rust) by replacing only the import statement; the transform calls remain unchanged.

The Python implementation in `python/universalfft/core.py` wraps NumPy:

```python
import numpy as np
import universalfft as ufft

N, dt = 1024, 1e-3                       # 1024 samples, 1 ms spacing
t = np.arange(N) * dt
x = np.cos(2 * np.pi * 10 * t)           # 10 Hz cosine

# Discrete spectrum (CD forward)
f  = ufft.freqs(N, dt)
Xc = ufft.fft_cd(x)                      # spike at ±10 Hz with amplitude 0.5

# Impulse response of a lowpass filter (CC inverse)
H = ufft.low_pass_response(f, f_cut=30.0)
h = ufft.ifft_cc(H, dt).real             # sinc-like; h.sum() * dt ≈ 1.0
```

### 2.3 Backend Normalisation Corrections

Each language backend is responsible for one transformation: mapping its host library's native FFT output to Boteler-compliant values. Table 2 summarises the corrections.

**Table 2.** Per-backend normalisation corrections applied by UniversalFFT.

| Backend | Native forward | Correction to CC | Correction to CD |
|---|:---:|:---:|:---:|
| Python (NumPy) | raw | × Δt | / N |
| MATLAB / GNU Octave | raw | × Δt | / N |
| R | raw | × Δt | / N |
| C / C++ / Fortran / Rust / JavaScript / CUDA | raw (self-contained) | × Δt | / N |
| Julia (FFTW.jl) | raw | × Δt | / N |
| R `fft(inverse=TRUE)` | raw (no 1/N) | × Δf directly | identity |
| IDL `FFT(x, -1)` | 1/N (IDL convention) | × N·Δt | identity |

The self-contained backends (C, C++, Fortran, Rust, JavaScript, CUDA/HIP) implement a radix-2 Cooley–Tukey FFT in-language, with no external library dependency. This ensures that the library is portable to environments where FFTW or equivalent is unavailable.

### 2.4 Cross-Language Validation Framework

Numerical consistency is verified by a purpose-built cross-language validation framework in `tests/cross_language/`. The workflow is:

1. `generate_reference.py` (Python) produces reference `.npz` files containing inputs and expected outputs for four canonical test cases: CC roundtrip, CD roundtrip, CD cosine spectrum, and CC impulse response integral.
2. A demo program in each language reads the same input, computes the transforms, and writes results as CSV.
3. `validate_all.py` loads every CSV and checks that each value agrees with the Python reference to within 10⁻⁹ in absolute error.

All 11 backends currently pass this check. Figure 1 shows the agreement visually for the CC roundtrip test across all languages.

*[Figure 1: `cross_language_demo.png` — cross-language numerical agreement for the CC roundtrip test. Each bar represents the maximum absolute deviation from the Python reference; all are below 10⁻⁹.]*

---

## 3. Illustrative Examples

### 3.1 Impulse Response of a Lowpass Filter

A brick-wall lowpass filter with cutoff at f_c = 30 Hz is defined on N = 1024 samples at Δt = 1 ms. Its transfer function H[k] = 1 for |f_k| ≤ f_c and 0 otherwise.

```python
import numpy as np
import universalfft as ufft

N, dt, f_cut = 1024, 1e-3, 30.0
f = ufft.freqs(N, dt)
H = ufft.low_pass_response(f, f_cut)
h = ufft.ifft_cc(H, dt).real            # CC inverse: correct physical scaling

integral = h.sum() * dt                 # ≈ 1.0000000000  (Parseval/DC gain check)
```

The CC inverse produces a causal sinc-like impulse response whose discrete integral ∫ h(n)Δt equals 1.0 to floating-point precision. Using `ifft_cd` instead would yield a value of N·Δf = 1/(Δt) ≈ 1000, indistinguishable in shape but off by three orders of magnitude in amplitude — the silent error that motivated this library.

Figure 2 shows the impulse response and its cumulative integral converging to 1.

*[Figure 2: `impulse_response.png` — CC-convention impulse response of the 30 Hz lowpass filter. The cumulative integral (dashed) converges to 1.0.]*

### 3.2 Discrete Spectrum of a Cosine

A pure cosine of amplitude A = 2 at frequency f₁ = 10 Hz is analysed with the CD forward transform.

```python
t  = np.arange(N) * dt
x  = 2.0 * np.cos(2 * np.pi * 10 * t)
Xc = ufft.fft_cd(x)

# Spike at the positive-frequency bin k=10:
print(abs(Xc[10]))   # → 1.0  (= A/2, Fourier-series coefficient)
```

The CD forward delivers exactly A/2 = 1.0 at ±f₁. Using `fft_cc` instead would yield Δt = 0.001 at those bins — dimensionally a spectrum of a continuous function, not a dimensionless Fourier coefficient.

---

## 4. Impact

UniversalFFT is designed for the geoscience modelling community, where Fourier transforms serve as the primary tool for working with magnetometer time series, transfer functions derived from magnetotelluric surveys, and electromagnetic induction calculations for GIC assessment in power grids and submarine cables.

The library addresses a concrete reproducibility gap. When one group uses Python and another uses MATLAB or Julia, their single-transform results may differ by a factor of N·Δt or N·Δf even when both implement the same physical model — because each language's raw FFT has a different native normalisation. UniversalFFT makes those results numerically identical to 10⁻⁹, verified automatically by the cross-language validation suite.

The library is also useful outside geoscience for any workflow that interprets single DFT outputs physically: spectral energy estimates, matched-filter amplitudes, or system-identification impulse responses. The explicit CC/CD distinction provides a natural vocabulary for documenting *which* transform was intended, reducing the risk of silent convention drift across software versions or collaborators.

The software is available under the MIT licence, carries a permanent DOI on Zenodo (placeholder: 10.5281/zenodo.XXXXXXX), and is documented at https://universalfft.readthedocs.io.

---

## 5. Conclusions

Boteler (2012) proved that two distinct FFT normalisation conventions are appropriate for geoscience applications, each for physically distinct reasons, and derived the correct scaling factors from first principles. The fourteen years since that publication have left practitioners with the mathematical answer but no production software to apply it. UniversalFFT provides that software: a six-function API, implemented across 11 languages, that corrects each host library's normalisation internally and validates cross-language numerical agreement to better than 10⁻⁹. The library requires no modification to use the correct convention — it is enforced by construction.

---

## Acknowledgements

The author thanks David H. Boteler (Natural Resources Canada) for the foundational mathematical conventions implemented here and for encouraging the development of a dedicated software companion to Boteler (2012).

---

## References

- Boteler, D.H. (2012). On Choosing Fourier Transforms for Practical Geoscience Applications. *International Journal of Geosciences*, **3**, 952–959. doi:10.4236/ijg.2012.325096

- Bracewell, R.N. (1978). *The Fourier Transform and Its Applications*. McGraw-Hill.

- Brigham, E.O. (1974). *The Fast Fourier Transform*. Prentice-Hall.

- Harris, C.R., et al. (2020). Array programming with NumPy. *Nature*, **585**, 357–362. doi:10.1038/s41586-020-2649-2

- Price, A.T. (1962). The theory of magnetotelluric methods when the source field is considered. *Journal of Geophysical Research*, **67**, 1907–1918.

- Virtanen, P., et al. (2020). SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. *Nature Methods*, **17**, 261–272. doi:10.1038/s41592-019-0686-2

- Ward, S.H. & Hohmann, G.W. (1988). Electromagnetic Theory for Geophysical Applications. In M.N. Nabighian (Ed.), *Electromagnetic Methods in Applied Geophysics*, Vol. 1. Society of Exploration Geophysicists.
