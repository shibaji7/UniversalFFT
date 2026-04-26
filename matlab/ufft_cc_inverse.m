function x = ufft_cc_inverse(X, dt)
% UFFT_CC_INVERSE  Continuous–Continuous inverse DFT (Boteler 2012, Eq. 21b).
%
%   x = UFFT_CC_INVERSE(X, dt)
%
%   x(n) = sum_{k=0}^{N-1} X(k) * exp(+i*2*pi*k*n/N) * df
%          where df = 1 / (N * dt)
%
%   MATLAB's ifft divides by N:
%       ifft(X) = (1/N) * sum X(k) * exp(+i*2*pi*(k-1)*(n-1)/N)
%
%   So the Boteler CC inverse = ifft(X) * N * df = ifft(X) / dt
%
%   Parameters
%   ----------
%   X  : complex frequency-domain vector (FFT-output order), length N.
%   dt : sampling interval of the original time series (seconds).
%
%   Returns
%   -------
%   x  : complex vector of N reconstructed time-domain samples.

N = numel(X);
df = 1.0 / (N * dt);
x = ifft(X(:)) * N * df;   % == ifft(X) / dt
end
