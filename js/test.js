/**
 * test.js — Unit tests for universalfft.js using Node built-in test runner.
 * Run: node --test test.js  (Node ≥ 18)
 */

import { test } from "node:test";
import assert from "node:assert/strict";
import { fftCC, ifftCC, fftCD, ifftCD, fftCCReal, freqs, lowPassResponse } from "./universalfft.js";

const N   = 128;
const DT  = 1e-3;
const F1  = 60.0;
const A   = 2.5;
const TOL = 1e-9;

function cosine() {
  const x = new Float64Array(N);
  for (let i = 0; i < N; ++i) x[i] = A * Math.cos(2 * Math.PI * F1 * i * DT);
  return x;
}

function maxErr(a, b) {
  let e = 0;
  for (let i = 0; i < a.length; ++i) e = Math.max(e, Math.abs(a[i] - b[i]));
  return e;
}

test("CC round-trip", () => {
  const x    = cosine();
  const X    = fftCCReal(x, DT);
  const xrec = ifftCC(X, DT);
  assert.ok(maxErr(x, xrec.re) < TOL, `CC round-trip error too large: ${maxErr(x, xrec.re)}`);
});

test("CD round-trip", () => {
  const x    = cosine();
  const im   = new Float64Array(N);
  const X    = fftCD(x, im);
  const xrec = ifftCD(X);
  assert.ok(maxErr(x, xrec.re) < TOL, `CD round-trip error too large`);
});

test("CD spectrum spike at bin-aligned f1", () => {
  // Use f1=62.5 Hz so k = f1*N*DT = 62.5*128*1e-3 = 8 (exact bin, no leakage).
  const f1Exact = 62.5;
  const x  = new Float64Array(N).map((_, i) => A * Math.cos(2 * Math.PI * f1Exact * i * DT));
  const im = new Float64Array(N);
  const X  = fftCD(x, im);
  const f  = freqs(N, DT);
  let best = 0, bestDist = Infinity;
  for (let k = 0; k < N; ++k) {
    const d = Math.abs(f[k] - f1Exact);
    if (d < bestDist) { bestDist = d; best = k; }
  }
  const mag = Math.sqrt(X.re[best] ** 2 + X.im[best] ** 2);
  assert.ok(Math.abs(mag - A / 2) < 1e-9, `CD spike magnitude ${mag} ≠ A/2=${A/2}`);
});

test("impulse response integrates to 1", () => {
  const N2 = 2048, dt2 = 60.0, fc = 1 / 3600;
  const f  = freqs(N2, dt2);
  const H  = lowPassResponse(f, fc);
  const h  = ifftCC(H, dt2);
  let integral = 0;
  for (let i = 0; i < N2; ++i) integral += h.re[i] * dt2;
  assert.ok(Math.abs(integral - 1.0) < 1e-4, `Integral = ${integral}`);
});
