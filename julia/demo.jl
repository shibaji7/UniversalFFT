#!/usr/bin/env julia
"""
Julia cross-language validation demo.

Writes CSV files to ../tests/data/ matching the Python reference vectors.

Usage:
    julia --project=. demo.jl [output_dir]
"""
using UniversalFFT

outdir = length(ARGS) > 0 ? ARGS[1] : joinpath(@__DIR__, "../tests/data")
mkpath(outdir)

function write_complex_csv(path, z)
    open(path, "w") do io
        for v in z
            r, im = real(v), imag(v)
            println(io, im >= 0 ? "$(r)+$(im)j" : "$(r)$(im)j")
        end
    end
end

# ── 1. Cosine signal: N=128, dt=1e-3, f1=60 Hz, A=2.5 ──────────────────────
N, dt, f1, A = 128, 1e-3, 60.0, 2.5
t = (0:N-1) .* dt
x = A .* cos.(2π .* f1 .* t)

X_cc  = fft_cc(x, dt)
x_rec = ifft_cc(X_cc, dt)

write_complex_csv(joinpath(outdir, "julia_X_cc_cosine.csv"),    X_cc)
write_complex_csv(joinpath(outdir, "julia_x_cc_rec_cosine.csv"), x_rec)
println("[Julia] cosine CC forward/inverse written.")

# ── 2. LP filter: N=256, dt=60 s ────────────────────────────────────────────
N2, dt2, fc = 256, 60.0, 1/3600
using Random; rng = MersenneTwister(0)
x2 = randn(rng, N2)
f2 = freqs(N2, dt2)
H2 = complex.(low_pass_response(f2, fc))

y_filtered = fft_filter(x2, H2, dt2)
write_complex_csv(joinpath(outdir, "julia_y_filtered.csv"), y_filtered)
println("[Julia] filter output written.")

# ── 3. Impulse response ──────────────────────────────────────────────────────
H3 = complex.(low_pass_response(f2, fc))
h  = ifft_cc(H3, dt2)
write_complex_csv(joinpath(outdir, "julia_h_impulse.csv"), h)
integral = sum(real.(h)) * dt2
@printf("[Julia] impulse response written. Integral = %.8f (expect ≈1.0)\n", integral)
