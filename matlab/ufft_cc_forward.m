function X = ufft_cc_forward(x, dt)
% UFFT_CC_FORWARD  Continuous–Continuous forward DFT (Boteler 2012, Eq. 21a).
%
%   X = UFFT_CC_FORWARD(x, dt)
%
%   X(k) = sum_{n=0}^{N-1} x(n) * exp(-i*2*pi*k*n/N) * dt
%
%   MATLAB's built-in fft computes: X_raw(k) = sum x(n)*exp(-i*2*pi*(k-1)*(n-1)/N)
%   (1-indexed, no scaling), so the Boteler CC forward is:
%
%       X = fft(x) * dt
%
%   Parameters
%   ----------
%   x  : column or row vector of N real or complex time-domain samples.
%   dt : sampling interval in seconds.
%
%   Returns
%   -------
%   X  : complex vector of length N in FFT-output order:
%        [DC, +f1, +f2, ..., +f_Nyquist, -f_Nyquist+df, ..., -f1]

X = fft(x(:)) * dt;
end
