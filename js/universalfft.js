/**
 * universalfft.js — Boteler (2012) compliant FFT/IFFT for JavaScript (ES module).
 *
 * Uses a self-contained Cooley-Tukey radix-2 DIT FFT with bit-reversal.
 * Both forward and inverse are raw summations (no 1/N factor internally).
 *
 * Boteler mapping:
 *   CC forward  = rawFFT(x)  * dt
 *   CC inverse  = rawIFFT(X) * df          (df = 1/(N·dt))
 *   CD forward  = rawFFT(x)  / N
 *   CD inverse  = rawIFFT(X)               (raw sum — already correct)
 */

"use strict";

/* ── Cooley-Tukey radix-2 DIT in-place ──────────────────────────────────── */

/**
 * In-place FFT on paired Float64Arrays re[] and im[].
 * N must be a power of 2. inverse=true gives raw summation (no 1/N).
 */
function _fftInplace(re, im, inverse) {
  const N = re.length;
  /* Bit-reversal permutation */
  let j = 0;
  for (let i = 1; i < N; ++i) {
    let bit = N >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      [re[i], re[j]] = [re[j], re[i]];
      [im[i], im[j]] = [im[j], im[i]];
    }
  }
  /* Butterfly stages */
  const sign = inverse ? 1.0 : -1.0;
  for (let len = 2; len <= N; len <<= 1) {
    const ang = sign * 2.0 * Math.PI / len;
    const wRe = Math.cos(ang);
    const wIm = Math.sin(ang);
    for (let i = 0; i < N; i += len) {
      let uRe = 1.0, uIm = 0.0;
      for (let k = 0; k < len / 2; ++k) {
        const eRe = re[i + k];
        const eIm = im[i + k];
        const oRe = re[i + k + len / 2] * uRe - im[i + k + len / 2] * uIm;
        const oIm = re[i + k + len / 2] * uIm + im[i + k + len / 2] * uRe;
        re[i + k]           = eRe + oRe;
        im[i + k]           = eIm + oIm;
        re[i + k + len / 2] = eRe - oRe;
        im[i + k + len / 2] = eIm - oIm;
        const tmpRe = uRe * wRe - uIm * wIm;
        uIm = uRe * wIm + uIm * wRe;
        uRe = tmpRe;
      }
    }
  }
}

function _runFFT(reIn, imIn, inverse) {
  const re = Float64Array.from(reIn);
  const im = Float64Array.from(imIn);
  _fftInplace(re, im, inverse);
  return { re, im };
}

/* ── CC forward  (Boteler Eq. 21a) ─────────────────────────────────────── */

/**
 * CC forward transform.
 * @param {number[]|Float64Array} re  Real part of input
 * @param {number[]|Float64Array} im  Imaginary part of input (may be zeros)
 * @param {number} dt  Sample interval (s)
 * @returns {{ re: Float64Array, im: Float64Array }}
 */
export function fftCC(re, im, dt) {
  const r = _runFFT(re, im, false);
  for (let i = 0; i < r.re.length; ++i) { r.re[i] *= dt; r.im[i] *= dt; }
  return r;
}

/** Real-input overload — imaginary part assumed zero. */
export function fftCCReal(x, dt) {
  return fftCC(x, new Float64Array(x.length), dt);
}

/* ── CC inverse  (Boteler Eq. 21b) ─────────────────────────────────────── */

/**
 * CC inverse transform.
 * @param {{ re, im }} X  Frequency-domain spectrum
 * @param {number} dt  Original sample interval (s)
 * @returns {{ re: Float64Array, im: Float64Array }}
 */
export function ifftCC(X, dt) {
  const N  = X.re.length;
  const df = 1.0 / (N * dt);
  const r  = _runFFT(X.re, X.im, true);
  for (let i = 0; i < N; ++i) { r.re[i] *= df; r.im[i] *= df; }
  return r;
}

/* ── CD forward  (Boteler Eq. 22a) ─────────────────────────────────────── */

export function fftCD(re, im) {
  const N = re.length;
  const r = _runFFT(re, im, false);
  for (let i = 0; i < N; ++i) { r.re[i] /= N; r.im[i] /= N; }
  return r;
}

export function fftCDReal(x) {
  return fftCD(x, new Float64Array(x.length));
}

/* ── CD inverse  (Boteler Eq. 22b) ─────────────────────────────────────── */

export function ifftCD(X) {
  return _runFFT(X.re, X.im, true);
}

/* ── Frequency bin array ─────────────────────────────────────────────────── */

/**
 * FFT frequency bins in Hz, FFT-output order (same as numpy.fft.fftfreq).
 * @param {number} N   Number of samples
 * @param {number} dt  Sample interval (s)
 * @returns {Float64Array}
 */
export function freqs(N, dt) {
  const df = 1.0 / (N * dt);
  const f  = new Float64Array(N);
  for (let k = 0; k < N; ++k)
    f[k] = k < N / 2 ? k * df : (k - N) * df;
  return f;
}

/* ── Low-pass brick-wall response ────────────────────────────────────────── */

export function lowPassResponse(f, fc) {
  const re = new Float64Array(f.length);
  const im = new Float64Array(f.length);
  for (let k = 0; k < f.length; ++k)
    re[k] = Math.abs(f[k]) <= fc ? 1.0 : 0.0;
  return { re, im };
}

/* ── Filter via CC pair ──────────────────────────────────────────────────── */

export function fftFilter(x, H, dt) {
  const X  = fftCCReal(x, dt);
  const XH = { re: new Float64Array(x.length), im: new Float64Array(x.length) };
  for (let k = 0; k < x.length; ++k) {
    XH.re[k] = X.re[k] * H.re[k] - X.im[k] * H.im[k];
    XH.im[k] = X.re[k] * H.im[k] + X.im[k] * H.re[k];
  }
  return ifftCC(XH, dt);
}
