/**
 * demo.cpp — C++ cross-language validation demo.
 * Compile:  make demo
 * Run:      ./ufft_demo_cpp [output_dir]
 */
#include "universalfft.hpp"
#include <cstdio>
#include <cstring>
#include <cmath>

static void write_csv(const char* path, const ufft::cvec& z) {
    FILE* fp = fopen(path, "w");
    if (!fp) { perror(path); return; }
    for (const auto& v : z) {
        double r = v.real(), im = v.imag();
        if (im >= 0) fprintf(fp, "%.17g+%.17gj\n", r, im);
        else         fprintf(fp, "%.17g%.17gj\n",  r, im);
    }
    fclose(fp);
}

int main(int argc, char* argv[]) {
    const char* outdir = (argc > 1) ? argv[1] : "../tests/data";

    // ── cosine: N=128, dt=1e-3, f1=60 Hz, A=2.5
    const std::size_t N  = 128;
    const double dt = 1e-3, f1 = 60.0, A = 2.5;
    ufft::rvec x(N);
    for (std::size_t n = 0; n < N; ++n)
        x[n] = A * std::cos(2.0 * ufft::PI * f1 * n * dt);

    auto X_cc  = ufft::fft_cc(x, dt);
    auto x_rec = ufft::ifft_cc(X_cc, dt);

    char path[512];
    snprintf(path, sizeof(path), "%s/cpp_X_cc_cosine.csv", outdir);
    write_csv(path, X_cc);
    printf("[C++] Wrote %s\n", path);

    snprintf(path, sizeof(path), "%s/cpp_x_cc_rec_cosine.csv", outdir);
    write_csv(path, x_rec);
    printf("[C++] Wrote %s\n", path);

    // ── impulse response: N=256, dt=60 s
    const std::size_t N2 = 256;
    const double dt2 = 60.0, fc = 1.0 / 3600.0;
    auto f2 = ufft::freqs(N2, dt2);
    auto H  = ufft::low_pass_response(f2, fc);
    auto h  = ufft::ifft_cc(H, dt2);

    snprintf(path, sizeof(path), "%s/cpp_h_impulse.csv", outdir);
    write_csv(path, h);
    printf("[C++] Wrote %s\n", path);

    double integral = 0.0;
    for (const auto& v : h) integral += v.real();
    integral *= dt2;
    printf("[C++] Impulse integral = %.8f  (expect ≈1.0)\n", integral);

    puts("[C++] Demo complete.");
    return 0;
}
