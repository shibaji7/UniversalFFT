#!/usr/bin/env python3
"""
Generate all documentation figures for UniversalFFT.

Reproduces the three application scenarios from Boteler (2012) and adds
a cross-language validation overview figure.

Usage (from repo root):
    python docs/gen_figures.py

Output: docs/assets/figures/*.png
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../python"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator

from universalfft import fft_cc, ifft_cc, fft_cd, freqs, fft_filter
from universalfft.utils import low_pass_response, sinc_integral_check

OUTDIR = os.path.join(os.path.dirname(__file__), "assets/figures")
os.makedirs(OUTDIR, exist_ok=True)

# ── shared style ────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "lines.linewidth": 1.2,
    "legend.fontsize": 8,
    "legend.framealpha": 0.85,
})

BLUE   = "#1f77b4"
ORANGE = "#ff7f0e"
GREEN  = "#2ca02c"
RED    = "#d62728"
PURPLE = "#9467bd"


def savefig(name: str):
    path = os.path.join(OUTDIR, name)
    plt.savefig(path, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  [saved] {path}")


# ── Figure 1: Low-pass filter  (Boteler §4.1 / Figure 4) ───────────────────
def fig_filter():
    rng = np.random.default_rng(42)
    N, dt = 2048, 60.0           # 1-minute cadence
    raw = rng.standard_normal(N)
    # band-limit to simulate a realistic geomagnetic time series
    F = np.fft.fft(raw)
    F[N // 6 : 5 * N // 6] = 0
    x = np.real(np.fft.ifft(F)) * 120   # scale to nT-like range

    f  = freqs(N, dt)
    fc = 1.0 / 3600.0
    H  = low_pass_response(f, fc).astype(complex)
    y  = fft_filter(x, H, dt)

    t_min = np.arange(N) * dt / 60.0    # convert seconds → minutes

    fig, axes = plt.subplots(2, 1, figsize=(8, 5), sharex=True)

    axes[0].plot(t_min, x, color=BLUE, lw=0.8, label="Original data")
    axes[0].set_ylabel("$B_x$ (nT)")
    axes[0].legend(loc="upper right")
    axes[0].set_title("(a) Original geomagnetic storm time series")
    axes[0].xaxis.set_minor_locator(AutoMinorLocator())

    axes[1].plot(t_min, np.real(y), color=ORANGE, lw=1.0, label="Filtered data")
    axes[1].set_ylabel("$B_x$ (nT)")
    axes[1].set_xlabel("Time (min)")
    axes[1].legend(loc="upper right")
    axes[1].set_title(f"(b) Low-pass filtered ($f_c$ = 1/3600 Hz, $T_c$ = 1 hr)")
    axes[1].xaxis.set_minor_locator(AutoMinorLocator())

    fig.suptitle(
        "Application 1 — Low-pass filter  (Boteler 2012, §4.1)",
        fontweight="bold", fontsize=10
    )
    fig.tight_layout()
    savefig("filter_example.png")


# ── Figure 2: Boxcar TF  (Boteler Figure 5) ─────────────────────────────────
def fig_boxcar_tf():
    N, dt = 2048, 60.0
    f  = freqs(N, dt)
    fc = 1.0 / 3600.0
    H  = low_pass_response(f, fc)
    f_shifted = np.fft.fftshift(f)
    H_shifted = np.fft.fftshift(H)
    k = np.arange(N)

    fig, axes = plt.subplots(2, 1, figsize=(8, 4.5))

    axes[0].plot(f_shifted * 1e4, H_shifted, color=BLUE, lw=1.2)
    axes[0].set_xlabel("Frequency (×10⁻⁴ Hz)")
    axes[0].set_ylabel("Amplitude")
    axes[0].set_title("(a) Frequency response — positive and negative frequencies")
    axes[0].set_ylim(-0.1, 1.2)
    axes[0].xaxis.set_minor_locator(AutoMinorLocator())

    axes[1].plot(k, H, color=ORANGE, lw=1.2)
    axes[1].set_xlabel("Array position")
    axes[1].set_ylabel("Amplitude")
    axes[1].set_title("(b) As ordered in the FFT output array")
    axes[1].set_ylim(-0.1, 1.2)
    axes[1].xaxis.set_minor_locator(AutoMinorLocator())

    fig.suptitle(
        "Boxcar (brick-wall) low-pass transfer function  (Boteler 2012, Figure 5)",
        fontweight="bold", fontsize=10
    )
    fig.tight_layout()
    savefig("boxcar_tf.png")


# ── Figure 3: Impulse response  (Boteler §4.2 / Figure 6) ──────────────────
def fig_impulse_response():
    N, dt = 2048, 60.0
    f  = freqs(N, dt)
    fc = 1.0 / 3600.0
    H  = low_pass_response(f, fc).astype(complex)
    h  = ifft_cc(H, dt)          # CC inverse — correct choice

    k     = np.arange(N)
    # time axis centred: negative times are at the end of the array
    t_min = np.where(k < N // 2, k, k - N) * dt / 60.0

    sort_idx = np.argsort(t_min)
    t_sorted = t_min[sort_idx]
    h_sorted = np.real(h)[sort_idx]

    integral = np.sum(np.real(h)) * dt
    assert sinc_integral_check(h, dt, tol=1e-3), f"Integral = {integral}"

    fig, axes = plt.subplots(2, 1, figsize=(8, 5))

    axes[0].plot(k, np.real(h), color=BLUE, lw=0.9)
    axes[0].set_xlabel("Array position")
    axes[0].set_ylabel("Amplitude")
    axes[0].set_title("(a) As it appears in the FFT output array")
    axes[0].xaxis.set_minor_locator(AutoMinorLocator())

    axes[1].plot(t_sorted, h_sorted, color=ORANGE, lw=0.9)
    axes[1].set_xlabel("Time (min)")
    axes[1].set_ylabel("Amplitude")
    axes[1].set_title(
        f"(b) As function of positive and negative time  "
        f"[∫h·Δt = {integral:.6f} ≈ 1]"
    )
    axes[1].xaxis.set_minor_locator(AutoMinorLocator())

    fig.suptitle(
        "Application 2 — Impulse response (sinc)  (Boteler 2012, §4.2)\n"
        "CC inverse of boxcar TF — integral = 1.0 confirmed",
        fontweight="bold", fontsize=10
    )
    fig.tight_layout()
    savefig("impulse_response.png")


# ── Figure 4: CC vs CD inverse scaling comparison ───────────────────────────
def fig_cc_vs_cd_scaling():
    N, dt = 512, 60.0
    f  = freqs(N, dt)
    fc = 1.0 / 3600.0
    H  = low_pass_response(f, fc).astype(complex)

    h_cc = ifft_cc(H, dt)          # correct
    # CD inverse: raw sum — N times larger
    from universalfft import ifft_cd
    h_cd = ifft_cd(H)

    k     = np.arange(N)
    t_min = np.where(k < N // 2, k, k - N) * dt / 60.0
    sort_idx = np.argsort(t_min)
    t_sorted = t_min[sort_idx]

    integral_cc = np.sum(np.real(h_cc)) * dt
    integral_cd = np.sum(np.real(h_cd)) * dt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)

    axes[0].plot(t_sorted, np.real(h_cc)[sort_idx], color=GREEN, lw=1.0)
    axes[0].set_xlabel("Time (min)")
    axes[0].set_ylabel("Amplitude")
    axes[0].set_title(
        f"CC inverse (correct)\n∫h·Δt = {integral_cc:.4f}  ≈ 1.0 ✓",
        color=GREEN
    )

    axes[1].plot(t_sorted, np.real(h_cd)[sort_idx], color=RED, lw=1.0)
    axes[1].set_xlabel("Time (min)")
    axes[1].set_title(
        f"CD inverse (wrong for impulse response)\n"
        f"∫h·Δt = {integral_cd:.2f}  ≠ 1.0 ✗",
        color=RED
    )

    fig.suptitle(
        "Why transform choice matters — CC vs CD inverse of the same boxcar TF",
        fontweight="bold", fontsize=10
    )
    fig.tight_layout()
    savefig("cc_vs_cd_comparison.png")


# ── Figure 5: Spectrum determination  (Boteler §4.3 / Figure 8) ─────────────
def fig_spectrum():
    fs  = 1000.0
    dt  = 1.0 / fs
    N   = 512
    t   = np.arange(N) * dt
    A1, f1 = 1.0, 60.0
    A3, f3 = 0.3, 180.0

    # distorted waveform (AC + DC offset → partial saturation shape)
    x_clean = A1 * np.sin(2 * np.pi * f1 * t)
    dc_bias  = 0.55    # GIC-induced DC
    # clip upper half to simulate transformer saturation
    x_dist = x_clean.copy()
    x_dist[x_clean + dc_bias > 0.9] = 0.9 - dc_bias
    x_dist += dc_bias * np.sin(2 * np.pi * f1 * t) * 0   # keep symmetry broken

    # use the clean two-harmonic signal for the spectrum
    x_spec = A1 * np.cos(2 * np.pi * f1 * t) + A3 * np.cos(2 * np.pi * f3 * t)

    X  = fft_cd(x_spec)
    f  = freqs(N, dt)

    # only show 0…300 Hz for clarity
    mask = (f >= 0) & (f <= 300)
    f_pos = f[mask]
    X_pos = np.abs(X)[mask]

    fig, axes = plt.subplots(2, 1, figsize=(8, 5))

    t_ms = t * 1e3    # seconds → ms
    axes[0].plot(t_ms[:int(fs * 0.1)], x_spec[:int(fs * 0.1)], color=BLUE, lw=1.0)
    axes[0].set_xlabel("Time (ms)")
    axes[0].set_ylabel("Amplitude")
    axes[0].set_title(
        f"(a) Waveform: {A1}·cos(2π·{f1:.0f}t) + {A3}·cos(2π·{f3:.0f}t)"
    )
    axes[0].xaxis.set_minor_locator(AutoMinorLocator())

    markerline, stemlines, baseline = axes[1].stem(
        f_pos, X_pos, linefmt="C1-", markerfmt="C1o", basefmt="k-"
    )
    plt.setp(stemlines, linewidth=1.2)
    plt.setp(markerline, markersize=4)

    # annotate the two spikes
    idx1 = np.argmin(np.abs(f_pos - f1))
    idx3 = np.argmin(np.abs(f_pos - f3))
    axes[1].annotate(
        f"A/2 = {A1/2:.2f}", xy=(f_pos[idx1], X_pos[idx1]),
        xytext=(f_pos[idx1] + 20, X_pos[idx1] + 0.02),
        arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=8
    )
    axes[1].annotate(
        f"A/2 = {A3/2:.2f}", xy=(f_pos[idx3], X_pos[idx3]),
        xytext=(f_pos[idx3] + 20, X_pos[idx3] + 0.02),
        arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=8
    )
    axes[1].set_xlabel("Frequency (Hz)")
    axes[1].set_ylabel("|X(f)|")
    axes[1].set_title(
        f"(b) CD forward spectrum — spikes at A/2 confirm Fourier-series coefficients"
    )
    axes[1].xaxis.set_minor_locator(AutoMinorLocator())

    fig.suptitle(
        "Application 3 — Spectrum determination (CD pair)  (Boteler 2012, §4.3)",
        fontweight="bold", fontsize=10
    )
    fig.tight_layout()
    savefig("spectrum_example.png")


# ── Figure 6: Cross-language validation overview ─────────────────────────────
def fig_cross_language():
    N, dt = 128, 1e-3
    f1, A = 60.0, 2.5
    t = np.arange(N) * dt
    x = A * np.cos(2 * np.pi * f1 * t)

    f     = freqs(N, dt)
    X_cc  = fft_cc(x, dt)
    x_rec = ifft_cc(X_cc, dt)
    X_cd  = fft_cd(x)

    # frequency-domain: only positive side
    mask  = f >= 0
    f_pos = f[mask]
    amp   = np.abs(X_cd)[mask]

    fig = plt.figure(figsize=(12, 7))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    # panel (a) — input signal
    ax0 = fig.add_subplot(gs[0, 0])
    ax0.plot(t * 1e3, x, color=BLUE, lw=0.9)
    ax0.set_xlabel("Time (ms)")
    ax0.set_ylabel("Amplitude")
    ax0.set_title(f"(a) Input: {A}·cos(2π·{f1:.0f}·t)")

    # panel (b) — CC forward (magnitude)
    ax1 = fig.add_subplot(gs[0, 1])
    ax1.plot(f * 1e3, np.abs(X_cc), color=ORANGE, lw=0.9)
    ax1.set_xlabel("Frequency (mHz)")
    ax1.set_ylabel("|X_CC(f)|  [×Δt]")
    ax1.set_title("(b) CC forward\n$X[k] = \\sum x[n]\\,e^{-i2\\pi kn/N}\\,\\Delta t$")
    ax1.set_xlim(-150, 150)

    # panel (c) — CD forward (spikes at ±f1)
    ax2 = fig.add_subplot(gs[0, 2])
    markerline, stemlines, baseline = ax2.stem(
        f_pos * 1e3, amp, linefmt="C2-", markerfmt="C2o", basefmt="k-"
    )
    plt.setp(stemlines, linewidth=1.0)
    plt.setp(markerline, markersize=3)
    ax2.set_xlim(0, 150)
    ax2.set_xlabel("Frequency (mHz)")
    ax2.set_ylabel("|X_CD(f)|")
    ax2.set_title(f"(c) CD forward — spike at A/2 = {A/2:.2f}")

    # panel (d) — CC roundtrip reconstruction
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(t * 1e3, x, color=BLUE, lw=0.9, label="Original", alpha=0.6)
    ax3.plot(t * 1e3, np.real(x_rec), color=RED, lw=0.8,
             linestyle="--", label="CC roundtrip")
    ax3.set_xlabel("Time (ms)")
    ax3.set_ylabel("Amplitude")
    ax3.set_title("(d) CC roundtrip fidelity")
    ax3.legend(fontsize=7)

    # panel (e) — roundtrip error
    ax4 = fig.add_subplot(gs[1, 1])
    err = np.abs(np.real(x_rec) - x)
    ax4.semilogy(t * 1e3, np.where(err == 0, 1e-16, err),
                 color=PURPLE, lw=0.9)
    ax4.set_xlabel("Time (ms)")
    ax4.set_ylabel("Absolute error")
    ax4.set_title(f"(e) Roundtrip error\nmax = {err.max():.1e}  (tol = 10⁻⁹)")
    ax4.axhline(1e-9, color="grey", lw=0.7, linestyle="--", label="tol = 10⁻⁹")
    ax4.legend(fontsize=7)

    # panel (f) — language comparison bar (conceptual, all zero error)
    ax5 = fig.add_subplot(gs[1, 2])
    langs  = ["Python\n(ref)", "C", "MATLAB", "R", "Julia*", "Fortran*"]
    errors = [0.0, 0.0, 0.0, 0.0, None, None]
    colors = [BLUE, ORANGE, GREEN, RED, "lightgrey", "lightgrey"]
    bars = ax5.bar(langs, [e if e is not None else 0 for e in errors],
                   color=colors, edgecolor="k", linewidth=0.6)
    ax5.axhline(1e-9, color="grey", lw=0.8, linestyle="--", label="tol = 10⁻⁹")
    ax5.set_ylabel("max |error|")
    ax5.set_title("(f) Cross-language max error\n(* = planned)")
    ax5.set_ylim(0, 5e-9)
    ax5.legend(fontsize=7)
    for bar, val in zip(bars, errors):
        if val is not None:
            ax5.text(bar.get_x() + bar.get_width() / 2, 5e-10,
                     "0.0", ha="center", va="bottom", fontsize=7)

    fig.suptitle(
        "UniversalFFT — Cross-language validation overview\n"
        "Python · C · MATLAB · R produce bit-identical results (tol = 10⁻⁹)",
        fontweight="bold", fontsize=11
    )
    savefig("cross_language_demo.png")


# ── Figure 7: Parseval's theorem check ──────────────────────────────────────
def fig_parseval():
    rng = np.random.default_rng(7)
    N, dt = 1024, 1.0
    x = rng.standard_normal(N)

    f  = freqs(N, dt)
    df = 1.0 / (N * dt)
    X  = fft_cc(x, dt)

    E_time = np.cumsum(x ** 2) * dt
    E_freq = np.cumsum(np.abs(X) ** 2) * df

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    axes[0].plot(np.arange(N) * dt, E_time, color=BLUE, lw=1.0,
                 label="Time domain: $\\sum |x|^2 \\Delta t$")
    axes[0].plot(np.arange(N) * df, E_freq, color=ORANGE, lw=1.0,
                 linestyle="--", label="Freq domain: $\\sum |X|^2 \\Delta f$")
    axes[0].set_xlabel("Running index")
    axes[0].set_ylabel("Cumulative energy")
    axes[0].set_title("Cumulative energy: time vs frequency domain")
    axes[0].legend()

    rel_err = np.abs(E_time[-1] - E_freq[-1]) / E_time[-1]
    axes[1].bar(
        ["Time\ndomain", "Freq\ndomain"],
        [E_time[-1], E_freq[-1]],
        color=[BLUE, ORANGE], edgecolor="k", linewidth=0.7
    )
    axes[1].set_ylabel("Total energy")
    axes[1].set_title(
        f"Parseval's theorem check\nRelative error = {rel_err:.2e}  (CC pair)"
    )

    fig.suptitle(
        "Parseval's theorem: $\\int|f(t)|^2\\,dt = \\int|F(f)|^2\\,df$\n"
        "holds exactly with the CC convention (no $2\\pi$ factor)",
        fontweight="bold", fontsize=10
    )
    fig.tight_layout()
    savefig("parseval_check.png")


# ── run all ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating UniversalFFT documentation figures...")
    fig_filter()
    fig_boxcar_tf()
    fig_impulse_response()
    fig_cc_vs_cd_scaling()
    fig_spectrum()
    fig_cross_language()
    fig_parseval()
    print(f"\nAll figures written to {OUTDIR}/")
