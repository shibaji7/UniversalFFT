%% universalfft.m — Boteler (2012) compliant FFT/IFFT for Octave (and MATLAB).
%%
%% Octave's fft/ifft follow the same convention as MATLAB:
%%   fft(x)    = raw sum (no 1/N)
%%   ifft(X)   = raw_sum / N
%%
%% Boteler mapping:
%%   CC forward  = fft(x) * dt
%%   CC inverse  = ifft(X) * N * df  = ifft(X) / dt
%%   CD forward  = fft(x) / N
%%   CD inverse  = ifft(X) * N

function X = ufft_cc_forward(x, dt)
  %% X[k] = sum_n x[n] e^{-i2pifkn/N} * dt  (Boteler Eq. 21a)
  X = fft(x(:)) * dt;
end

function x = ufft_cc_inverse(X, dt)
  %% x[n] = sum_k X[k] e^{+i2pifkn/N} * df,  df=1/(N*dt)  (Boteler Eq. 21b)
  N  = numel(X);
  df = 1 / (N * dt);
  x  = ifft(X(:)) * N * df;   %% == ifft(X) / dt
end

function X = ufft_cd_forward(x)
  %% X[k] = (1/N) sum_n x[n] e^{-i2pifkn/N}  (Boteler Eq. 22a)
  N = numel(x);
  X = fft(x(:)) / N;
end

function x = ufft_cd_inverse(X)
  %% x[n] = sum_k X[k] e^{+i2pifkn/N}  (Boteler Eq. 22b)
  N = numel(X);
  x = ifft(X(:)) * N;   %% undo Octave's 1/N
end

function f = ufft_freqs(N, dt)
  %% Frequency bins in Hz, FFT-output order (matches numpy.fft.fftfreq).
  df = 1 / (N * dt);
  k  = (0 : N-1)';
  f  = zeros(N, 1);
  f(k < N/2)  = k(k < N/2) .* df;
  f(k >= N/2) = (k(k >= N/2) - N) .* df;
end

function H = ufft_low_pass(f, fc)
  %% Brick-wall low-pass frequency response.
  H = complex(double(abs(f) <= fc));
end

function y = ufft_filter(x, H, dt)
  %% Filter via CC pair.
  X = ufft_cc_forward(x, dt);
  y = ufft_cc_inverse(X .* H(:), dt);
end
