//! Rust cross-language validation demo.
//! Writes CSV files to ../tests/data/ matching the Python reference vectors.
//!
//! Run: cargo run -- [output_dir]

use universalfft::*;
use std::fs::File;
use std::io::{BufWriter, Write as IoWrite};
use std::path::Path;

fn write_csv(path: &str, z: &[C64]) {
    let f = File::create(path).expect("cannot create file");
    let mut w = BufWriter::new(f);
    for v in z {
        if v.im >= 0.0 {
            writeln!(w, "{:+.17e}+{:.17e}j", v.re, v.im).unwrap();
        } else {
            writeln!(w, "{:+.17e}{:.17e}j", v.re, v.im).unwrap();
        }
    }
}

fn main() {
    let outdir = std::env::args().nth(1).unwrap_or_else(|| "../tests/data".to_string());
    std::fs::create_dir_all(&outdir).unwrap();

    // ── cosine: N=128, dt=1e-3, f1=60 Hz, A=2.5 ───────────────────────────
    let (n, dt, f1, a) = (128usize, 1e-3f64, 60.0f64, 2.5f64);
    let x: Vec<f64> = (0..n)
        .map(|i| a * (2.0 * std::f64::consts::PI * f1 * i as f64 * dt).cos())
        .collect();
    let cx: Vec<C64> = x.iter().map(|&v| C64::new(v, 0.0)).collect();

    let X_cc  = fft_cc(&cx, dt);
    let x_rec = ifft_cc(&X_cc, dt);

    write_csv(&format!("{}/rust_X_cc_cosine.csv",    outdir), &X_cc);
    write_csv(&format!("{}/rust_x_cc_rec_cosine.csv", outdir), &x_rec);
    println!("[Rust] cosine CC forward/inverse written.");

    // ── impulse response: N=256, dt=60 s ────────────────────────────────────
    let (n2, dt2) = (256usize, 60.0f64);
    let f2  = freqs(n2, dt2);
    let H   = low_pass_response(&f2, 1.0 / 3600.0);
    let h   = ifft_cc(&H, dt2);

    write_csv(&format!("{}/rust_h_impulse.csv", outdir), &h);
    let integral: f64 = h.iter().map(|v| v.re).sum::<f64>() * dt2;
    println!("[Rust] impulse response written. Integral = {:.8}  (expect ≈1.0)", integral);

    println!("[Rust] Demo complete.");
}
