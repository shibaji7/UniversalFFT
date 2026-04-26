"""
    UniversalFFT

Boteler (2012) compliant FFT/IFFT wrappers for Julia.

Julia/FFTW.jl convention (what we wrap)
----------------------------------------
    fft(x)   → X_raw[k] = Σ x[n] e^{-i2π(k-1)(n-1)/N}   (no scaling, 1-indexed)
    ifft(X)  → x_raw[n] = (1/N) Σ X[k] e^{+i2π(k-1)(n-1)/N}

Boteler mapping
---------------
    CC forward  = fft(x)   * dt
    CC inverse  = ifft(X)  * N * df  =  ifft(X) / dt    (since N·df = 1/dt)
    CD forward  = fft(x)   / N
    CD inverse  = ifft(X)  * N                           (undo the 1/N)
"""
module UniversalFFT

using FFTW

export fft_cc, ifft_cc, fft_cd, ifft_cd, fft_filter, freqs,
       low_pass_response, sinc_integral_check, next_power_of_two

# ── CC forward  (Boteler Eq. 21a) ───────────────────────────────────────────
"""
    fft_cc(x, dt) -> Vector{ComplexF64}

Continuous–Continuous forward DFT.

    X[k] = Σ x[n] e^{-i2πkn/N} Δt

# Arguments
- `x`  : time-domain samples (real or complex), length N
- `dt` : sampling interval in seconds
"""
function fft_cc(x::AbstractVector, dt::Real)
    return FFTW.fft(complex.(x)) .* dt
end

# ── CC inverse  (Boteler Eq. 21b) ───────────────────────────────────────────
"""
    ifft_cc(X, dt) -> Vector{ComplexF64}

Continuous–Continuous inverse DFT.

    x[n] = Σ X[k] e^{+i2πkn/N} Δf,    Δf = 1/(N·Δt)

FFTW.ifft applies 1/N, so the full scale factor is N·Δf = 1/Δt.
"""
function ifft_cc(X::AbstractVector, dt::Real)
    N  = length(X)
    df = 1.0 / (N * dt)
    return FFTW.ifft(complex.(X)) .* (N * df)   # == ifft(X) / dt
end

# ── CD forward  (Boteler Eq. 22a) ───────────────────────────────────────────
"""
    fft_cd(x, dt=nothing) -> Vector{ComplexF64}

Continuous–Discrete forward DFT.

    X[k] = (1/N) Σ x[n] e^{-i2πkn/N}

`dt` is unused; accepted for API symmetry.
"""
function fft_cd(x::AbstractVector, dt=nothing)
    N = length(x)
    return FFTW.fft(complex.(x)) ./ N
end

# ── CD inverse  (Boteler Eq. 22b) ───────────────────────────────────────────
"""
    ifft_cd(X, dt=nothing) -> Vector{ComplexF64}

Continuous–Discrete inverse DFT (raw summation).

    x[n] = Σ X[k] e^{+i2πkn/N}

FFTW.ifft divides by N, so we undo that factor.
"""
function ifft_cd(X::AbstractVector, dt=nothing)
    N = length(X)
    return FFTW.ifft(complex.(X)) .* N
end

# ── Frequency bin array ───────────────────────────────────────────────────────
"""
    freqs(N, dt) -> Vector{Float64}

FFT frequency bin array in Hz, in FFT-output order (matching numpy.fft.fftfreq).

    f[k] =  k/(N·dt)       for k = 0 … N/2-1
    f[k] = (k-N)/(N·dt)    for k = N/2 … N-1
"""
function freqs(N::Integer, dt::Real)
    df = 1.0 / (N * dt)
    f  = collect(0:N-1) .* df
    half = N ÷ 2
    f[half+1:end] .= (collect(half:N-1) .- N) .* df
    return f
end

# ── Convenience filter ────────────────────────────────────────────────────────
"""
    fft_filter(x, H, dt) -> Vector{ComplexF64}

Apply transfer function `H` (at FFT frequency bins) to `x` via the CC pair.
Either CC or CD pair gives identical results for a filter pair (Boteler §4.1).
"""
function fft_filter(x::AbstractVector, H::AbstractVector, dt::Real)
    length(x) == length(H) || throw(
        DimensionMismatch("x and H must have the same length")
    )
    return ifft_cc(fft_cc(x, dt) .* complex.(H), dt)
end

# ── Utility helpers ───────────────────────────────────────────────────────────
"""
    low_pass_response(f, fc) -> Vector{Float64}

Brick-wall low-pass transfer function: 1.0 where |f| ≤ fc, 0.0 elsewhere.
"""
low_pass_response(f::AbstractVector, fc::Real) = Float64.(abs.(f) .<= fc)

"""
    sinc_integral_check(h, dt; tol=1e-4) -> Bool

Verify that the impulse response integrates to 1.0 (Boteler §4.2).
"""
function sinc_integral_check(h::AbstractVector, dt::Real; tol::Real=1e-4)
    return abs(sum(real.(h)) * dt - 1.0) < tol
end

"""
    next_power_of_two(n) -> Int

Smallest power of two ≥ n.
"""
next_power_of_two(n::Integer) = n <= 1 ? 1 : 2^(ceil(Int, log2(n)))

end # module
