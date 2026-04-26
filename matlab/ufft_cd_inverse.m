function x = ufft_cd_inverse(X, ~)
% UFFT_CD_INVERSE  Continuous–Discrete inverse DFT (Boteler 2012, Eq. 22b).
%
%   x = UFFT_CD_INVERSE(X)   or   x = UFFT_CD_INVERSE(X, dt)
%
%   x(n) = sum_{k=0}^{N-1} X(k) * exp(+i*2*pi*k*n/N)   (raw summation)
%
%   MATLAB's ifft divides by N, so we undo that: x = N * ifft(X).

N = numel(X);
x = N * ifft(X(:));
end
