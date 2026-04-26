# CUDA / HIP — universalfft.cuh / .cu

Sources:
[`cuda/universalfft.cuh`](https://github.com/shibaji7/UniversalFFT/blob/main/cuda/universalfft.cuh),
[`cuda/universalfft.cu`](https://github.com/shibaji7/UniversalFFT/blob/main/cuda/universalfft.cu)

## Backend

Wraps **cuFFT** (NVIDIA) or **hipFFT** (AMD ROCm) selected at compile time.
Both libraries are raw-sum in both directions (no 1/N).
A custom GPU kernel `ufft_scale_kernel` applies the Boteler scaling factor in parallel on the device.

## Compile-Time Backend Selection

```c
// CUDA (default)
nvcc -O2 -o ufft_demo_cu  demo.cu universalfft.cu -lcufft

// HIP (AMD ROCm)
hipcc -DUFFT_USE_HIP -O2 -o ufft_demo_hip demo.cu universalfft.cu -lhipfft
```

## Boteler Mapping

| Function | Scale applied on GPU |
|----------|---------------------|
| `ufft_cc_forward_host(…, dt)` | `dt` |
| `ufft_cc_inverse_host(…, dt)` | `df = 1/(N·dt)` |
| `ufft_cd_forward_host(…)` | `1/N` |
| `ufft_cd_inverse_host(…)` | `1.0` (no scaling) |

## Host API

All functions allocate/free device memory internally (convenience API).

```c
int ufft_cc_forward_host(
    const double* in_re, const double* in_im,
    double* out_re, double* out_im,
    int N, double dt);

int ufft_cc_inverse_host(
    const double* in_re, const double* in_im,
    double* out_re, double* out_im,
    int N, double dt);

int ufft_cd_forward_host(
    const double* in_re, const double* in_im,
    double* out_re, double* out_im,
    int N);

int ufft_cd_inverse_host(
    const double* in_re, const double* in_im,
    double* out_re, double* out_im,
    int N);

void ufft_freqs_cpu(double* f, int N, double dt);

int ufft_filter_host(
    const double* in_re, const double* in_im,
    const double* H_re,  const double* H_im,
    double* out_re, double* out_im,
    int N, double dt);
```

All functions return `0` on success, `-1` on error.

## GPU Kernel

```c
__global__ void ufft_scale_kernel(UFFT_COMPLEX* d, double scale, int N);
```

Launched as `ceil(N/256)` blocks × 256 threads. Multiplies each element by `scale`.
