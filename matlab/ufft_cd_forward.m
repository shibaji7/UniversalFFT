function X = ufft_cd_forward(x, ~)
% UFFT_CD_FORWARD  Continuous–Discrete forward DFT (Boteler 2012, Eq. 22a).
%
%   X = UFFT_CD_FORWARD(x)   or   X = UFFT_CD_FORWARD(x, dt)
%
%   X(k) = (1/N) * sum_{n=0}^{N-1} x(n) * exp(-i*2*pi*k*n/N)
%
%   dt is accepted but unused (kept for API symmetry with ufft_cc_forward).
%
%   Use for spectrum determination of periodic signals: a cosine of
%   amplitude A yields spikes of A/2 at ±f1 (Boteler §4.3).

N = numel(x);
X = fft(x(:)) / N;
end
