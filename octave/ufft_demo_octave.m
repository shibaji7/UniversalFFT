%% ufft_demo_octave.m — UniversalFFT Octave demo.
%%
%% Run: octave --no-gui ufft_demo_octave.m [output_dir]
%%
%% Writes CSV files for cross-language validation.

%% Source the function definitions
source('universalfft.m');

if nargin >= 1
  outdir = argv(){1};
else
  outdir = '../tests/data';
end
[~, ~, ~] = mkdir(outdir);

function write_csv(path, z)
  fid = fopen(path, 'w');
  for i = 1:numel(z)
    re_part = real(z(i));
    im_part = imag(z(i));
    if im_part >= 0
      fprintf(fid, '%+.17e+%.17ej\n', re_part, im_part);
    else
      fprintf(fid, '%+.17e%.17ej\n', re_part, im_part);
    end
  end
  fclose(fid);
end

%% ── cosine: N=128, dt=1e-3, f1=60 Hz, A=2.5 ───────────────────────────
N  = 128; dt = 1e-3; f1 = 60; A = 2.5;
t  = (0 : N-1)' * dt;
x  = A * cos(2*pi*f1*t);

X_cc  = ufft_cc_forward(x, dt);
x_rec = ufft_cc_inverse(X_cc, dt);

write_csv(fullfile(outdir, 'octave_X_cc_cosine.csv'),     X_cc);
write_csv(fullfile(outdir, 'octave_x_cc_rec_cosine.csv'), x_rec);
fprintf('[Octave] cosine CC forward/inverse written.\n');

%% ── impulse response: N=256, dt=60 s ──────────────────────────────────
N2  = 256; dt2 = 60; fc = 1/3600;
f2  = ufft_freqs(N2, dt2);
H   = ufft_low_pass(f2, fc);
h   = ufft_cc_inverse(H, dt2);

write_csv(fullfile(outdir, 'octave_h_impulse.csv'), h);

integral = sum(real(h)) * dt2;
fprintf('[Octave] impulse response integral = %.8f  (expect ≈ 1.0)\n', integral);

fprintf('[Octave] Demo complete.\n');
