/**
 * universalfft.cu — Boteler (2012) compliant FFT/IFFT implementation
 *                   for CUDA (cuFFT) and HIP (hipFFT).
 *
 * Build:
 *   CUDA:  nvcc -O2 -o ufft_demo_cu  demo.cu universalfft.cu -lcufft
 *   HIP:   hipcc -O2 -o ufft_demo_hip demo.cu universalfft.cu -lhipfft
 */

#include "universalfft.cuh"
#include <cstdlib>
#include <cstring>
#include <cmath>

/* ── Internal helper: alloc device, copy H→D, exec FFT, copy D→H, free ─── */

static int exec_fft_host(
        const double* in_re, const double* in_im,
        double*       out_re, double*      out_im,
        int N, int direction, double scale)
{
    size_t bytes = (size_t)N * sizeof(UFFT_COMPLEX);

    UFFT_COMPLEX* d_buf = nullptr;
    if (cudaMalloc((void**)&d_buf, bytes) != cudaSuccess) return -1;

    /* Pack interleaved host buffer */
    UFFT_COMPLEX* h_buf = (UFFT_COMPLEX*)malloc(bytes);
    if (!h_buf) { cudaFree(d_buf); return -1; }
    for (int i = 0; i < N; ++i) { h_buf[i].x = in_re[i]; h_buf[i].y = in_im[i]; }

    cudaMemcpy(d_buf, h_buf, bytes, cudaMemcpyHostToDevice);

    UFFT_HANDLE plan;
    UFFT_PLAN_1D(&plan, N, UFFT_Z2Z, 1);
    UFFT_EXEC_Z2Z(plan, d_buf, d_buf, direction);
    UFFT_DESTROY(plan);

    /* Scale kernel */
    if (scale != 1.0) {
        int threads = 256;
        int blocks  = (N + threads - 1) / threads;
        ufft_scale_kernel<<<blocks, threads>>>(d_buf, scale, N);
    }

    cudaMemcpy(h_buf, d_buf, bytes, cudaMemcpyDeviceToHost);
    cudaFree(d_buf);

    for (int i = 0; i < N; ++i) { out_re[i] = h_buf[i].x; out_im[i] = h_buf[i].y; }
    free(h_buf);
    return 0;
}

/* ── Public API ──────────────────────────────────────────────────────────── */

extern "C" {

int ufft_cc_forward_host(
        const double* in_re, const double* in_im,
        double* out_re,      double* out_im,
        int N, double dt)
{
    /* CC forward: X[k] = Σ x[n] e^{-i2πkn/N} · Δt  (Boteler Eq. 21a) */
    return exec_fft_host(in_re, in_im, out_re, out_im, N, UFFT_FORWARD, dt);
}

int ufft_cc_inverse_host(
        const double* in_re, const double* in_im,
        double* out_re,      double* out_im,
        int N, double dt)
{
    /* CC inverse: x[n] = Σ X[k] e^{+i2πkn/N} · Δf,  Δf = 1/(N·Δt) */
    double df = 1.0 / ((double)N * dt);
    return exec_fft_host(in_re, in_im, out_re, out_im, N, UFFT_INVERSE, df);
}

int ufft_cd_forward_host(
        const double* in_re, const double* in_im,
        double* out_re,      double* out_im,
        int N)
{
    /* CD forward: X[k] = (1/N) Σ x[n] e^{-i2πkn/N}  (Boteler Eq. 22a) */
    return exec_fft_host(in_re, in_im, out_re, out_im, N, UFFT_FORWARD,
                         1.0 / (double)N);
}

int ufft_cd_inverse_host(
        const double* in_re, const double* in_im,
        double* out_re,      double* out_im,
        int N)
{
    /* CD inverse: x[n] = Σ X[k] e^{+i2πkn/N}  (raw sum — scale = 1) */
    return exec_fft_host(in_re, in_im, out_re, out_im, N, UFFT_INVERSE, 1.0);
}

void ufft_freqs_cpu(double* f, int N, double dt)
{
    double df = 1.0 / ((double)N * dt);
    for (int k = 0; k < N; ++k)
        f[k] = (k < N / 2) ? (double)k * df : ((double)k - (double)N) * df;
}

int ufft_filter_host(
        const double* in_re, const double* in_im,
        const double* H_re,  const double* H_im,
        double* out_re,      double* out_im,
        int N, double dt)
{
    /* Forward CC */
    double* X_re = (double*)malloc(N * sizeof(double));
    double* X_im = (double*)malloc(N * sizeof(double));
    if (!X_re || !X_im) { free(X_re); free(X_im); return -1; }

    if (ufft_cc_forward_host(in_re, in_im, X_re, X_im, N, dt) != 0) {
        free(X_re); free(X_im); return -1;
    }

    /* Multiply spectrum by transfer function */
    double* XH_re = (double*)malloc(N * sizeof(double));
    double* XH_im = (double*)malloc(N * sizeof(double));
    if (!XH_re || !XH_im) { free(X_re); free(X_im); free(XH_re); free(XH_im); return -1; }

    for (int k = 0; k < N; ++k) {
        XH_re[k] = X_re[k] * H_re[k] - X_im[k] * H_im[k];
        XH_im[k] = X_re[k] * H_im[k] + X_im[k] * H_re[k];
    }
    free(X_re); free(X_im);

    /* Inverse CC */
    int ret = ufft_cc_inverse_host(XH_re, XH_im, out_re, out_im, N, dt);
    free(XH_re); free(XH_im);
    return ret;
}

} /* extern "C" */
