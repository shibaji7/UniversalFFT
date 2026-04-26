/**
 * demo.js — UniversalFFT JavaScript demo.
 *
 * Run with Node ≥ 14 (ESM):
 *   node --experimental-vm-modules demo.js [output_dir]
 * Or with Node ≥ 18 (native ESM):
 *   node demo.js [output_dir]
 *
 * Writes CSV files to ../tests/data/ for cross-language validation.
 */

import { fftCC, ifftCC, freqs, lowPassResponse } from "./universalfft.js";
import { writeFileSync, mkdirSync } from "node:fs";

const outdir = process.argv[2] ?? "../tests/data";
mkdirSync(outdir, { recursive: true });

function writeCSV(path, re, im) {
  const lines = [];
  for (let i = 0; i < re.length; ++i) {
    const sign = im[i] >= 0 ? "+" : "";
    lines.push(`${re[i].toExponential(17)}${sign}${im[i].toExponential(17)}j`);
  }
  writeFileSync(path, lines.join("\n") + "\n");
}

/* ── cosine: N=128, dt=1e-3, f1=60 Hz, A=2.5 ─────────────────────────── */
const N = 128, dt = 1e-3, f1 = 60.0, A = 2.5;
const x = new Float64Array(N).map((_, i) => A * Math.cos(2 * Math.PI * f1 * i * dt));

const Xcc  = fftCC(x, new Float64Array(N), dt);
const xrec = ifftCC(Xcc, dt);

writeCSV(`${outdir}/js_X_cc_cosine.csv`,     Xcc.re,  Xcc.im);
writeCSV(`${outdir}/js_x_cc_rec_cosine.csv`, xrec.re, xrec.im);
console.log("[JS] cosine CC forward/inverse written.");

/* ── impulse response: N=256, dt=60 s ──────────────────────────────────── */
const N2 = 256, dt2 = 60.0, fc = 1 / 3600;
const f2 = freqs(N2, dt2);
const H  = lowPassResponse(f2, fc);
const h  = ifftCC(H, dt2);

writeCSV(`${outdir}/js_h_impulse.csv`, h.re, h.im);

let integral = 0.0;
for (let i = 0; i < N2; ++i) integral += h.re[i] * dt2;
console.log(`[JS] impulse response integral = ${integral.toFixed(8)}  (expect ≈ 1.0)`);

console.log("[JS] Demo complete.");
