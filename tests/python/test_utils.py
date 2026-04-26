"""Unit tests for universalfft.utils."""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../python"))

import numpy as np
import pytest
from universalfft.utils import (
    next_power_of_two,
    pad_to_length,
    low_pass_response,
    sinc_integral_check,
)


class TestNextPowerOfTwo:
    @pytest.mark.parametrize("n, expected", [
        (1, 1), (2, 2), (3, 4), (4, 4),
        (5, 8), (100, 128), (128, 128), (129, 256),
    ])
    def test_values(self, n, expected):
        assert next_power_of_two(n) == expected

    def test_zero_raises(self):
        with pytest.raises(ValueError):
            next_power_of_two(0)

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            next_power_of_two(-5)


class TestPadToLength:
    def test_zero_pad(self):
        x = np.array([1.0, 2.0, 3.0])
        xp = pad_to_length(x, 6)
        assert xp.shape == (6,)
        np.testing.assert_array_equal(xp, [1, 2, 3, 0, 0, 0])

    def test_no_pad_needed(self):
        x = np.array([1.0, 2.0])
        xp = pad_to_length(x, 2)
        np.testing.assert_array_equal(xp, x)

    def test_reflect_pad(self):
        x = np.array([1.0, 2.0, 3.0])
        xp = pad_to_length(x, 5, mode="reflect")
        assert xp.shape == (5,)

    def test_too_short_raises(self):
        with pytest.raises(ValueError):
            pad_to_length(np.ones(5), 3)

    def test_unknown_mode_raises(self):
        with pytest.raises(ValueError, match="Unknown mode"):
            pad_to_length(np.ones(3), 6, mode="wrap")

    def test_2d_array(self):
        x = np.ones((4, 2))
        xp = pad_to_length(x, 6)
        assert xp.shape == (6, 2)
        np.testing.assert_array_equal(xp[4:], 0)


class TestLowPassResponse:
    def test_passes_dc(self):
        f = np.array([0.0, 0.5, 1.0, 2.0])
        H = low_pass_response(f, fc=1.0)
        assert H[0] == 1.0

    def test_cuts_above_fc(self):
        f = np.array([0.0, 0.5, 1.0, 2.0])
        H = low_pass_response(f, fc=1.0)
        assert H[3] == 0.0

    def test_keeps_fc(self):
        f = np.array([-1.0, 0.0, 1.0])
        H = low_pass_response(f, fc=1.0)
        np.testing.assert_array_equal(H, [1.0, 1.0, 1.0])

    def test_negative_freqs(self):
        f = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
        H = low_pass_response(f, fc=1.0)
        np.testing.assert_array_equal(H, [0.0, 1.0, 1.0, 1.0, 0.0])


class TestSincIntegralCheck:
    def test_sinc_from_boxcar(self):
        from universalfft import ifft_cc, freqs
        N, dt = 1024, 60.0
        f = freqs(N, dt)
        fc = 1.0 / 3600.0
        H = low_pass_response(f, fc).astype(complex)
        h = ifft_cc(H, dt)
        assert sinc_integral_check(h, dt, tol=1e-4)

    def test_fails_for_wrong_integral(self):
        h = np.zeros(100)
        dt = 1.0
        assert not sinc_integral_check(h, dt)

    def test_tight_tolerance(self):
        from universalfft import ifft_cc, freqs
        N, dt = 2048, 1.0
        f = freqs(N, dt)
        H = low_pass_response(f, fc=0.25).astype(complex)
        h = ifft_cc(H, dt)
        assert sinc_integral_check(h, dt, tol=1e-3)
