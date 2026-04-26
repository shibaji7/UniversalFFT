#!/usr/bin/env python3
"""
Cross-language validation: load reference vectors and compare language outputs.

Run after each language produces its own output CSVs.
Usage:
    python tests/cross_language/validate_all.py
"""

import sys
import os
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(__file__), "../data")
TOL = 1e-9


def load_npz(name):
    return np.load(os.path.join(DATA_DIR, name))


def load_csv(path):
    """Load a complex-per-line CSV written as 're+imj' or 're-imj'.

    Strips internal spaces so Fortran's 'ES24.17' width-padded output
    (e.g. '-4.1E-03+ 0.0E+00j') is also handled correctly.
    """
    values = []
    with open(path) as fh:
        for line in fh:
            line = line.strip().replace(" ", "")
            if not line:
                continue
            try:
                values.append(complex(line))
            except ValueError:
                # Fallback: split on last sign before 'j'
                s = line[:-1]  # strip trailing 'j'
                for idx in range(len(s) - 1, 0, -1):
                    if s[idx] in ("+", "-"):
                        try:
                            values.append(complex(float(s[:idx]),
                                                  float(s[idx:])))
                            break
                        except ValueError:
                            pass
    return np.array(values, dtype=complex)


def validate(lang: str, filename: str,
             ref_re_key: str, ref_im_key: str, npz_name: str) -> bool:
    csv_path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(csv_path):
        print(f"  [SKIP] {lang}: {filename} not found — run that language's demo first")
        return True  # don't fail CI on missing optional language output

    ref = load_npz(npz_name)
    result   = load_csv(csv_path)
    expected = ref[ref_re_key] + 1j * ref[ref_im_key]

    max_err = np.max(np.abs(result - expected))
    ok = max_err <= TOL
    print(f"  [{'PASS' if ok else 'FAIL'}] {lang}: max_abs_error = {max_err:.2e}  (tol={TOL:.0e})")
    return ok


def main():
    print("=" * 65)
    print("UniversalFFT cross-language validation")
    print(f"Tolerance: {TOL:.0e}   Data dir: {DATA_DIR}")
    print("=" * 65)

    all_pass = True

    # ── (1) CC forward of cosine signal ──────────────────────────────────────
    print("\n[1] CC forward transform of cosine signal")
    for lang, fname in [
        ("C",          "c_X_cc_cosine.csv"),
        ("MATLAB",     "matlab_X_cc_cosine.csv"),
        ("R",          "r_X_cc_cosine.csv"),
        ("Julia",      "julia_X_cc_cosine.csv"),
        ("Fortran",    "fortran_X_cc_cosine.csv"),
        ("C++",        "cpp_X_cc_cosine.csv"),
        ("Rust",       "rust_X_cc_cosine.csv"),
        ("CUDA",       "cuda_X_cc_cosine.csv"),
        ("JavaScript", "js_X_cc_cosine.csv"),
        ("Octave",     "octave_X_cc_cosine.csv"),
    ]:
        ok = validate(lang, fname, "X_cc_real", "X_cc_imag", "cosine_vectors.npz")
        all_pass = all_pass and ok

    # ── (2) CC round-trip (forward then inverse) ──────────────────────────────
    print("\n[2] CC round-trip (forward → inverse) on cosine signal")
    for lang, fname in [
        ("C",          "c_x_cc_rec_cosine.csv"),
        ("MATLAB",     "matlab_x_cc_rec_cosine.csv"),
        ("R",          "r_x_cc_rec_cosine.csv"),
        ("Julia",      "julia_x_cc_rec_cosine.csv"),
        ("Fortran",    "fortran_x_cc_rec_cosine.csv"),
        ("C++",        "cpp_x_cc_rec_cosine.csv"),
        ("Rust",       "rust_x_cc_rec_cosine.csv"),
        ("CUDA",       "cuda_x_cc_rec_cosine.csv"),
        ("JavaScript", "js_x_cc_rec_cosine.csv"),
        ("Octave",     "octave_x_cc_rec_cosine.csv"),
    ]:
        ok = validate(lang, fname, "x_cc_rec_real", "x_cc_rec_imag", "cosine_vectors.npz")
        all_pass = all_pass and ok

    # ── (3) Impulse response ──────────────────────────────────────────────────
    print("\n[3] Impulse response (CC inverse of boxcar TF)")
    for lang, fname in [
        ("C",          "c_h_impulse.csv"),
        ("MATLAB",     "matlab_h_impulse.csv"),
        ("R",          "r_h_impulse.csv"),
        ("Julia",      "julia_h_impulse.csv"),
        ("Fortran",    "fortran_h_impulse.csv"),
        ("C++",        "cpp_h_impulse.csv"),
        ("Rust",       "rust_h_impulse.csv"),
        ("CUDA",       "cuda_h_impulse.csv"),
        ("JavaScript", "js_h_impulse.csv"),
        ("Octave",     "octave_h_impulse.csv"),
    ]:
        ok = validate(lang, fname, "h_real", "h_imag", "impulse_response_vectors.npz")
        all_pass = all_pass and ok

    print("\n" + "=" * 65)
    print("OVERALL:", "PASS" if all_pass else "FAIL (see above)")
    print("=" * 65)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
