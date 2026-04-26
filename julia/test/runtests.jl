using Test
using UniversalFFT

const TOL = 1e-9

@testset "UniversalFFT.jl" begin

    # shared fixtures
    N, dt, f1, A = 100, 1e-3, 60.0, 3.0
    t = (0:N-1) .* dt
    x = A .* cos.(2π .* f1 .* t)

    @testset "CC roundtrip" begin
        X    = fft_cc(x, dt)
        xrec = ifft_cc(X, dt)
        @test maximum(abs.(real.(xrec) .- x)) < TOL
    end

    @testset "CD roundtrip" begin
        X    = fft_cd(x)
        xrec = ifft_cd(X)
        @test maximum(abs.(real.(xrec) .- x)) < TOL
    end

    @testset "CC forward = fft(x)*dt" begin
        import FFTW
        expected = FFTW.fft(complex.(x)) .* dt
        @test maximum(abs.(fft_cc(x, dt) .- expected)) < TOL
    end

    @testset "CD spectrum — spike amplitude A/2" begin
        X = fft_cd(x)
        f = freqs(N, dt)
        i1 = argmin(abs.(f .- f1))
        @test abs(abs(X[i1]) - A/2) < 1e-3
    end

    @testset "Impulse response integrates to 1 (CC inverse)" begin
        N2, dt2 = 2048, 60.0
        f2  = freqs(N2, dt2)
        H   = complex.(low_pass_response(f2, 1/3600))
        h   = ifft_cc(H, dt2)
        @test sinc_integral_check(h, dt2; tol=1e-4)
    end

    @testset "Filter equivalence: CC == CD" begin
        rng   = MersenneTwister(42)
        xr    = randn(rng, 2048)
        dt2   = 60.0
        f2    = freqs(2048, dt2)
        H     = complex.(low_pass_response(f2, 1/3600))
        y_cc  = ifft_cc(fft_cc(xr, dt2) .* H, dt2)
        y_cd  = ifft_cd(fft_cd(xr) .* H)
        @test maximum(abs.(y_cc .- y_cd)) < TOL
    end

    @testset "Parseval (CC)" begin
        df = 1.0 / (N * dt)
        X  = fft_cc(x, dt)
        E_t = sum(x .^ 2) * dt
        E_f = sum(abs.(X) .^ 2) * df
        @test abs(E_t - E_f) / E_t < 1e-10
    end

    @testset "freqs DC and Nyquist" begin
        f = freqs(256, 1e-3)
        @test f[1] ≈ 0.0
        @test abs(f[129]) ≈ 500.0
    end

    @testset "next_power_of_two" begin
        @test next_power_of_two(1)   == 1
        @test next_power_of_two(3)   == 4
        @test next_power_of_two(128) == 128
        @test next_power_of_two(129) == 256
    end

end
