#' UniversalFFT — Boteler (2012) conventions for R
#'
#' Wraps stats::fft() to implement the same two transform pairs as the
#' Python and MATLAB packages.
#'
#' R convention (stats::fft):
#'   fft(x)          → X_raw[k] = Σ x[n] exp(-i2π(k-1)(n-1)/N)  (no scaling)
#'   fft(X, inverse=TRUE) → Σ X[k] exp(+i2π(k-1)(n-1)/N)         (no 1/N)
#'
#' So:
#'   CC forward  = fft(x) * dt
#'   CC inverse  = fft(X, inverse=TRUE) * df   where df = 1/(N*dt)
#'   CD forward  = fft(x) / N
#'   CD inverse  = fft(X, inverse=TRUE) / N * N  = fft(X, inverse=TRUE)
#'
#' NOTE: R's fft(inverse=TRUE) does NOT divide by N (unlike numpy.ifft).
#' This actually simplifies the CD inverse: it is exactly fft(X, inverse=TRUE).


# ---- CC forward (Boteler Eq. 21a) ----------------------------------------

#' @param x  Numeric or complex vector of N time-domain samples.
#' @param dt Sampling interval in seconds.
#' @return   Complex vector of length N in FFT-output order.
ufft_cc_forward <- function(x, dt) {
  stats::fft(x) * dt
}


# ---- CC inverse (Boteler Eq. 21b) ----------------------------------------

#' @param X  Complex frequency-domain vector (FFT-output order), length N.
#' @param dt Sampling interval of the original time series (seconds).
#' @return   Complex vector of N reconstructed time-domain samples.
ufft_cc_inverse <- function(X, dt) {
  N  <- length(X)
  df <- 1.0 / (N * dt)
  # R's fft(inverse=TRUE) returns the raw summation (no 1/N division),
  # so we simply multiply by df.
  stats::fft(X, inverse = TRUE) * df
}


# ---- CD forward (Boteler Eq. 22a) ----------------------------------------

#' @param x  Numeric or complex vector of N time-domain samples.
#' @param dt Unused; kept for API symmetry.
#' @return   Complex vector of length N.
ufft_cd_forward <- function(x, dt = NULL) {
  N <- length(x)
  stats::fft(x) / N
}


# ---- CD inverse (Boteler Eq. 22b) ----------------------------------------

#' Raw summation x[n] = Σ X[k] e^{+i2πkn/N}.
#' R's fft(inverse=TRUE) already returns the raw summation without 1/N,
#' so this is a direct call.
#'
#' @param X  Complex frequency-domain vector (CD-scaled), length N.
#' @param dt Unused; kept for API symmetry.
#' @return   Complex vector of N time-domain samples.
ufft_cd_inverse <- function(X, dt = NULL) {
  stats::fft(X, inverse = TRUE)
}


# ---- Frequency bin array -------------------------------------------------

#' @param N  Number of samples.
#' @param dt Sampling interval in seconds.
#' @return   Real vector of length N containing frequency bins in Hz,
#'           in FFT-output order (same as numpy.fft.fftfreq).
ufft_freqs <- function(N, dt) {
  df <- 1.0 / (N * dt)
  k  <- 0:(N - 1)
  f  <- k * df
  f[k >= N %/% 2] <- (k[k >= N %/% 2] - N) * df
  f
}


# ---- Filter via CC pair ---------------------------------------------------

#' Apply a frequency-domain filter to a time series.
#'
#' @param x  Input time-domain vector (N samples).
#' @param H  Transfer function at FFT frequency bins (length N, complex).
#' @param dt Sampling interval in seconds.
#' @return   Complex filtered time-domain vector.
ufft_filter <- function(x, H, dt) {
  X <- ufft_cc_forward(x, dt)
  Y <- X * H
  ufft_cc_inverse(Y, dt)
}


# ---- Low-pass brick-wall transfer function --------------------------------

#' @param f  Frequency bin vector (from ufft_freqs).
#' @param fc Cut-off frequency in Hz.
#' @return   Real vector: 1.0 where |f| <= fc, 0.0 elsewhere.
ufft_lowpass <- function(f, fc) {
  as.numeric(abs(f) <= fc)
}
