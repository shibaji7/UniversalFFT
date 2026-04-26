/**
 * universalfft.cuh — Boteler (2012) compliant FFT/IFFT for CUDA and HIP.
 *
 * Supports both NVIDIA CUDA (cuFFT) and AMD HIP (hipFFT) through
 * compile-time macros.  Build with:
 *
 *   CUDA:  nvcc -O2 -o ufft_demo_cu demo.cu universalfft_cuda.cu -lcufft
 *   HIP:   hipcc -O2 -o ufft_demo_hip demo.cu universalfft_cuda.cu -lhipfft
 *
 * cuFFT / hipFFT convention (both are raw summations — no 1/N):
 *   FORWARD: Σ x[n] e^{-i2πkn/N}
 *   INVERSE: Σ X[k] e^{+i2πkn/N}
 *
 * Boteler mapping:
 *   CC forward  = cufft_forward * dt          (on device)
 *   CC inverse  = cufft_inverse * df          (Δf = 1/(N·Δt))
 *   CD forward  = cufft_forward / N
 *   CD inverse  = cufft_inverse               (raw sum already correct)
 */
#pragma once

/* ── Backend selection ────────────────────────────────────────────────────── */
#ifdef UFFT_USE_HIP
  #include <hipfft/hipfft.h>
  #define UFFT_HANDLE       hipfftHandle
  #define UFFT_PLAN_1D      hipfftPlan1d
  #define UFFT_EXEC_Z2Z     hipfftExecZ2Z
  #define UFFT_DESTROY      hipfftDestroy
  #define UFFT_Z2Z          HIPFFT_Z2Z
  #define UFFT_FORWARD      HIPFFT_FORWARD
  #define UFFT_INVERSE      HIPFFT_BACKWARD
  #define UFFT_COMPLEX      hipfftDoubleComplex
  #define CUDA_MALLOC       hipMalloc
  #define CUDA_FREE         hipFree
  #define CUDA_MEMCPY_H2D   hipMemcpy, hipMemcpyHostToDevice
  #define CUDA_MEMCPY_D2H   hipMemcpy, hipMemcpyDeviceToHost
  #define CUDA_LAUNCH_ARGS  <<< , >>>
#else
  #include <cufft.h>
  #include <cuda_runtime.h>
  #define UFFT_HANDLE       cufftHandle
  #define UFFT_PLAN_1D      cufftPlan1d
  #define UFFT_EXEC_Z2Z     cufftExecZ2Z
  #define UFFT_DESTROY      cufftDestroy
  #define UFFT_Z2Z          CUFFT_Z2Z
  #define UFFT_FORWARD      CUFFT_FORWARD
  #define UFFT_INVERSE      CUFFT_INVERSE
  #define UFFT_COMPLEX      cufftDoubleComplex
  #define CUDA_MALLOC       cudaMalloc
  #define CUDA_FREE         cudaFree
#endif

#include <cstddef>

/* ── GPU kernels for post-FFT scaling ─────────────────────────────────────── */

/**
 * Scale every element of a complex device array by a real scalar.
 * Launched with ceil(N/256) blocks of 256 threads.
 */
__global__ void ufft_scale_kernel(UFFT_COMPLEX* d, double scale, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        d[idx].x *= scale;
        d[idx].y *= scale;
    }
}

/* ── High-level wrappers (manage device memory internally) ────────────────── */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * CC forward on host arrays → host arrays.
 * Allocates/frees device memory each call (convenience API).
 * For performance-critical code, use the device-resident API below.
 */
int ufft_cc_forward_host(
    const double* in_re, const double* in_im,
    double* out_re,      double* out_im,
    int N, double dt);

int ufft_cc_inverse_host(
    const double* in_re, const double* in_im,
    double* out_re,      double* out_im,
    int N, double dt);

int ufft_cd_forward_host(
    const double* in_re, const double* in_im,
    double* out_re,      double* out_im,
    int N);

int ufft_cd_inverse_host(
    const double* in_re, const double* in_im,
    double* out_re,      double* out_im,
    int N);

/** Frequency bin array (Hz, FFT-output order) — CPU side only. */
void ufft_freqs_cpu(double* f, int N, double dt);

/** Filter via CC pair (host arrays). */
int ufft_filter_host(
    const double* in_re, const double* in_im,
    const double* H_re,  const double* H_im,
    double* out_re,      double* out_im,
    int N, double dt);

#ifdef __cplusplus
}
#endif
