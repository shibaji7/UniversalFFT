function y = ufft_filter(x, H, dt)
% UFFT_FILTER  Apply a linear filter H(f) to time series x (CC pair).
%
%   y = UFFT_FILTER(x, H, dt)
%
%   Performs: y = IFFT_CC( FFT_CC(x, dt) .* H, dt )
%
%   For a filter pair (forward + inverse) either CC or CD convention gives
%   identical results (Boteler 2012, §4.1), so this is purely CC.
%
%   Parameters
%   ----------
%   x  : input time-domain vector (N samples).
%   H  : transfer function vector at FFT frequency bins (N values).
%        Build using ufft_freqs and a mask function.
%   dt : sampling interval in seconds.
%
%   Returns
%   -------
%   y  : filtered time-domain vector (complex; take real() if input is real).

X = ufft_cc_forward(x, dt);
Y = X .* H(:);
y = ufft_cc_inverse(Y, dt);
end
