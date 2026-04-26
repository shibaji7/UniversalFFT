#!/usr/bin/env Rscript
# ufft_demo.R — R cross-language validation demo.
#
# Writes CSV files to ../tests/data/ for validate_all.py.
# Usage:
#   Rscript ufft_demo.R [output_dir]

source(file.path(dirname(sys.frame(1)$ofile %||% "."), "universalfft.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

args    <- commandArgs(trailingOnly = TRUE)
outdir  <- if (length(args) > 0) args[1] else "../tests/data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write_complex_csv <- function(path, z) {
  lines <- vapply(z, function(v) {
    r  <- Re(v)
    im <- Im(v)
    if (im >= 0) sprintf("%.17g+%.17gj", r, im)
    else         sprintf("%.17g%.17gj",  r, im)
  }, character(1))
  writeLines(lines, path)
}

# -------------------------------------------------------------------
# 1. Cosine signal
# -------------------------------------------------------------------
N  <- 128L;  dt <- 1e-3;  f1 <- 60.0;  A <- 2.5
t  <- (0:(N-1)) * dt
x  <- A * cos(2 * pi * f1 * t)

X_cc  <- ufft_cc_forward(x, dt)
x_rec <- ufft_cc_inverse(X_cc, dt)

write_complex_csv(file.path(outdir, "r_X_cc_cosine.csv"),    X_cc)
write_complex_csv(file.path(outdir, "r_x_cc_rec_cosine.csv"), x_rec)
cat("[R] cosine CC forward/inverse written.\n")

# -------------------------------------------------------------------
# 2. LP filter (N=256, dt=60 s)
# -------------------------------------------------------------------
N2  <- 256L;  dt2 <- 60.0;  fc <- 1/3600
set.seed(42)
x2  <- rnorm(N2)
f2  <- ufft_freqs(N2, dt2)
H2  <- ufft_lowpass(f2, fc) + 0i

y_filtered <- ufft_filter(x2, H2, dt2)
write_complex_csv(file.path(outdir, "r_y_filtered.csv"), y_filtered)
cat("[R] filter output written.\n")

# -------------------------------------------------------------------
# 3. Impulse response
# -------------------------------------------------------------------
H3 <- ufft_lowpass(f2, fc) + 0i
h  <- ufft_cc_inverse(H3, dt2)
write_complex_csv(file.path(outdir, "r_h_impulse.csv"), h)
integral <- sum(Re(h)) * dt2
cat(sprintf("[R] impulse response written. Integral = %.8f (expect ~1.0)\n", integral))
