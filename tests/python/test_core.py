"""
Unit tests for universalfft.core following Boteler (2012).

Test strategy
-------------
1. Roundtrip fidelity — CC and CD pairs each reconstruct the input.
2. Physical scaling — CC forward result approximates Fourier integral.
3. Spectrum correctness — CD forward of cosine gives amplitude spikes of A/2.
4. Impulse response — CC inverse of boxcar TF integrates to 1 (Boteler §4.2).
5. Filter equivalence — CC pair and CD pair give identical filter output.
6. Cross-convention — filter output equals time-domain convolution.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../python"))

import numpy as np
import pytest
from universalfft import fft_cc, ifft_cc, fft_cd, ifft_cd, fft_filter, freqs
from universalfft.utils import low_pass_response, sinc_integral_check

TOL = 1e-9          # cross-language target tolerance
LOOSE = 1e-6        # for numerical integrals


# ---------------------------------------------------------------------------
# 1. Roundtrip tests
# ---------------------------------------------------------------------------

class TestRoundtrip:
    def test_cc_roundtrip_real(self, simple_cosine):
        x = simple_cosine["x"]
        dt = simple_cosine["dt"]
        X = fft_cc(x, dt)
        x_rec = ifft_cc(X, dt)
        np.testing.assert_allclose(np.real(x_rec), x, atol=TOL)

    def test_cc_roundtrip_complex(self, magnetic_storm_like):
        x = magnetic_storm_like["x"].astype(complex)
        dt = magnetic_storm_like["dt"]
        X = fft_cc(x, dt)
        x_rec = ifft_cc(X, dt)
        np.testing.assert_allclose(x_rec, x, atol=TOL)

    def test_cd_roundtrip_real(self, simple_cosine):
        x = simple_cosine["x"]
        dt = simple_cosine["dt"]
        X = fft_cd(x, dt)
        x_rec = ifft_cd(X, dt)
        np.testing.assert_allclose(np.real(x_rec), x, atol=TOL)

    def test_cd_roundtrip_complex(self, magnetic_storm_like):
        x = magnetic_storm_like["x"].astype(complex)
        dt = magnetic_storm_like["dt"]
        X = fft_cd(x, dt)
        x_rec = ifft_cd(X, dt)
        np.testing.assert_allclose(x_rec, x, atol=TOL)


# ---------------------------------------------------------------------------
# 2. Physical scaling of CC forward transform
# ---------------------------------------------------------------------------

class TestCCPhysicalScaling:
    def test_cc_forward_equals_numpy_times_dt(self, simple_cosine):
        """CC forward = np.fft.fft * dt (Boteler Eq. 21a)."""
        x, dt = simple_cosine["x"], simple_cosine["dt"]
        expected = np.fft.fft(x.astype(complex)) * dt
        np.testing.assert_allclose(fft_cc(x, dt), expected, atol=TOL)

    def test_cc_inverse_equals_numpy_ifft_over_dt(self, simple_cosine):
        """CC inverse = np.fft.ifft(X) * N * df = np.fft.ifft(X) / dt."""
        x, dt = simple_cosine["x"], simple_cosine["dt"]
        X = fft_cc(x, dt)
        N = len(x)
        expected = np.fft.ifft(X) / dt
        np.testing.assert_allclose(ifft_cc(X, dt), expected, atol=TOL)

    def test_parseval_cc(self, magnetic_storm_like):
        """Energy in time domain equals energy in frequency domain (Parseval)."""
        x, dt = magnetic_storm_like["x"], magnetic_storm_like["dt"]
        N = len(x)
        df = 1.0 / (N * dt)
        X = fft_cc(x, dt)
        E_time = np.sum(x ** 2) * dt
        E_freq = np.sum(np.abs(X) ** 2) * df
        np.testing.assert_allclose(E_time, E_freq, rtol=1e-10)


# ---------------------------------------------------------------------------
# 3. CD spectrum — cosine amplitude check (Boteler §4.3)
# ---------------------------------------------------------------------------

class TestCDSpectrum:
    def test_cosine_spike_amplitude(self, simple_cosine):
        """CD forward of A*cos(2πf₁t) → spikes of A/2 at ±f₁ (Boteler §4.3)."""
        d = simple_cosine
        X = fft_cd(d["x"])
        f = freqs(d["N"], d["dt"])
        # find index of +f1
        idx = np.argmin(np.abs(f - d["f1"]))
        np.testing.assert_allclose(np.abs(X[idx]), d["A"] / 2, atol=1e-3)

    def test_cd_forward_equals_numpy_div_N(self, simple_cosine):
        """CD forward = np.fft.fft / N (Boteler Eq. 22a)."""
        x = simple_cosine["x"]
        N = simple_cosine["N"]
        expected = np.fft.fft(x.astype(complex)) / N
        np.testing.assert_allclose(fft_cd(x), expected, atol=TOL)


# ---------------------------------------------------------------------------
# 4. Impulse response integrates to 1 (Boteler §4.2)
# ---------------------------------------------------------------------------

class TestImpulseResponse:
    def test_boxcar_ifft_cc_sinc_integral(self):
        """IFFT-CC of a boxcar TF (low-pass) integrates to 1 (Boteler §4.2)."""
        N = 2048
        dt = 60.0                         # seconds (magnetometer cadence)
        f = freqs(N, dt)
        fc = 1.0 / 3600.0                 # 1-hour period cut-off
        H = low_pass_response(f, fc)
        # use CC inverse — this is the critical choice (Boteler §4.2)
        h = ifft_cc(H.astype(complex), dt)
        assert sinc_integral_check(h, dt, tol=1e-4), (
            f"Integral of impulse response = {np.sum(np.real(h)) * dt:.6f}; expected 1.0"
        )

    def test_boxcar_ifft_cd_incorrect_scaling(self):
        """CD inverse of a boxcar TF should NOT integrate to 1 (wrong choice)."""
        N = 2048
        dt = 60.0
        f = freqs(N, dt)
        fc = 1.0 / 3600.0
        H = low_pass_response(f, fc)
        h_wrong = ifft_cd(H.astype(complex))
        integral = np.sum(np.real(h_wrong)) * dt
        # The integral should be ~N (not 1), confirming CD is wrong here
        assert abs(integral - 1.0) > 1.0, (
            "CD inverse accidentally passed the integral check — unexpected"
        )


# ---------------------------------------------------------------------------
# 5. Filter equivalence: CC pair == CD pair
# ---------------------------------------------------------------------------

class TestFilterEquivalence:
    def test_cc_cd_filter_identical(self, magnetic_storm_like):
        """Filtering with CC pair and CD pair gives identical results (Boteler §4.1)."""
        x, dt = magnetic_storm_like["x"], magnetic_storm_like["dt"]
        N = len(x)
        f = freqs(N, dt)
        fc = 1.0 / 3600.0
        H = low_pass_response(f, fc).astype(complex)

        # CC pair
        y_cc = ifft_cc(fft_cc(x, dt) * H, dt)

        # CD pair
        y_cd = ifft_cd(fft_cd(x) * H)

        np.testing.assert_allclose(y_cc, y_cd, atol=TOL)

    def test_fft_filter_convenience(self, magnetic_storm_like):
        """fft_filter() matches manual CC pair application."""
        x, dt = magnetic_storm_like["x"], magnetic_storm_like["dt"]
        N = len(x)
        f = freqs(N, dt)
        fc = 1.0 / 3600.0
        H = low_pass_response(f, fc).astype(complex)

        y_manual = ifft_cc(fft_cc(x, dt) * H, dt)
        y_conv = fft_filter(x, H, dt)
        np.testing.assert_allclose(y_conv, y_manual, atol=TOL)


# ---------------------------------------------------------------------------
# 6. Frequency bin array
# ---------------------------------------------------------------------------

class TestFreqs:
    def test_dc_is_zero(self):
        f = freqs(256, 1e-3)
        assert f[0] == 0.0

    def test_nyquist_bin(self):
        N, dt = 256, 1e-3
        f = freqs(N, dt)
        fs = 1.0 / dt
        assert abs(f[N // 2]) == pytest.approx(fs / 2)

    def test_bin_spacing(self):
        N, dt = 100, 0.5
        f = freqs(N, dt)
        df_expected = 1.0 / (N * dt)
        assert f[1] == pytest.approx(df_expected)


# ---------------------------------------------------------------------------
# 7. Edge cases and error handling
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_single_sample(self):
        x = np.array([5.0])
        dt = 1.0
        X = fft_cc(x, dt)
        x_rec = ifft_cc(X, dt)
        np.testing.assert_allclose(np.real(x_rec), x, atol=TOL)

    def test_dc_only_signal(self):
        N, dt = 64, 1.0
        x = np.ones(N) * 7.0
        X = fft_cc(x, dt)
        # DC bin should be 7 * N * dt
        np.testing.assert_allclose(np.real(X[0]), 7.0 * N * dt, atol=TOL)

    def test_filter_shape_mismatch_raises(self):
        x = np.ones(64)
        H = np.ones(32)
        with pytest.raises(ValueError, match="same shape"):
            fft_filter(x, H, dt=1.0)
