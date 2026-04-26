/**
 * demo.cu — UniversalFFT CUDA/HIP demo.
 *
 * Validates CC round-trip and writes CSV files to ../tests/data/.
 *
 * Build:
 *   CUDA:  nvcc -O2 -o ufft_demo_cu  demo.cu universalfft.cu -lcufft
 *   HIP:   hipcc -O2 -o ufft_demo_hip demo.cu universalfft.cu -lhipfft
 */

#include "universalfft.cuh"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

static void write_csv(const char* path, const double* re, const double* im, int N)
{
    FILE* fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return; }
    for (int i = 0; i < N; ++i)
        fprintf(fp, "%+.17e%+.17ej\n", re[i], im[i]);
    fclose(fp);
}

int main(int argc, char* argv[])
{
    const char* outdir = (argc > 1) ? argv[1] : "../tests/data";

    /* ── cosine: N=128, dt=1e-3, f1=60 Hz, A=2.5 ─────────────────────── */
    const int N = 128;
    const double dt = 1e-3, f1 = 60.0, A = 2.5;

    double* x_re  = (double*)calloc(N, sizeof(double));
    double* x_im  = (double*)calloc(N, sizeof(double));
    double* X_re  = (double*)calloc(N, sizeof(double));
    double* X_im  = (double*)calloc(N, sizeof(double));
    double* xr_re = (double*)calloc(N, sizeof(double));
    double* xr_im = (double*)calloc(N, sizeof(double));

    for (int i = 0; i < N; ++i)
        x_re[i] = A * cos(2.0 * M_PI * f1 * i * dt);

    ufft_cc_forward_host(x_re, x_im, X_re,  X_im,  N, dt);
    ufft_cc_inverse_host(X_re, X_im, xr_re, xr_im, N, dt);

    /* Round-trip check */
    double max_err = 0.0;
    for (int i = 0; i < N; ++i) max_err = fmax(max_err, fabs(xr_re[i] - x_re[i]));
    printf("[CUDA] CC round-trip max error: %.3e  (expect < 1e-9)\n", max_err);

    /* Write CSV */
    char path[512];
    snprintf(path, sizeof(path), "%s/cuda_X_cc_cosine.csv", outdir);
    write_csv(path, X_re, X_im, N);

    snprintf(path, sizeof(path), "%s/cuda_x_cc_rec_cosine.csv", outdir);
    write_csv(path, xr_re, xr_im, N);
    printf("[CUDA] cosine CC forward/inverse written.\n");

    /* ── impulse response: N=256, dt=60 s ─────────────────────────────── */
    const int N2 = 256;
    const double dt2 = 60.0, fc = 1.0 / 3600.0;

    double* f2    = (double*)calloc(N2, sizeof(double));
    double* H_re  = (double*)calloc(N2, sizeof(double));
    double* H_im  = (double*)calloc(N2, sizeof(double));
    double* h_re  = (double*)calloc(N2, sizeof(double));
    double* h_im  = (double*)calloc(N2, sizeof(double));

    ufft_freqs_cpu(f2, N2, dt2);
    for (int k = 0; k < N2; ++k)
        H_re[k] = (fabs(f2[k]) <= fc) ? 1.0 : 0.0;

    ufft_cc_inverse_host(H_re, H_im, h_re, h_im, N2, dt2);

    double integral = 0.0;
    for (int i = 0; i < N2; ++i) integral += h_re[i] * dt2;
    printf("[CUDA] impulse response integral = %.8f  (expect ≈ 1.0)\n", integral);

    snprintf(path, sizeof(path), "%s/cuda_h_impulse.csv", outdir);
    write_csv(path, h_re, h_im, N2);

    free(x_re); free(x_im); free(X_re); free(X_im); free(xr_re); free(xr_im);
    free(f2); free(H_re); free(H_im); free(h_re); free(h_im);

    printf("[CUDA] Demo complete.\n");
    return 0;
}
