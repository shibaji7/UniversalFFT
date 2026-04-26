function f = ufft_freqs(N, dt)
% UFFT_FREQS  FFT frequency bin array in Hz (FFT output order).
%
%   f = UFFT_FREQS(N, dt)
%
%   Matches numpy.fft.fftfreq(N, d=dt):
%     f(k) =  k/(N*dt)       for k = 0 .. N/2-1
%     f(k) = (k-N)/(N*dt)    for k = N/2 .. N-1
%
%   Use fftshift(f) to get [-f_Nyquist, ..., f_Nyquist] ordering.

df = 1.0 / (N * dt);
k  = (0 : N-1)';
f  = k * df;
half = N / 2;
f(k >= half) = (k(k >= half) - N) * df;
end
