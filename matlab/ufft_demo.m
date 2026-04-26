%% ufft_demo.m — MATLAB cross-language validation demo
%
% Produces CSV files matching the Python reference vectors in tests/data/.
% Run from the matlab/ directory:
%   >> ufft_demo
%
% Output CSVs are written to ../tests/data/.

outdir = '../tests/data';
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% 1. Cosine signal: N=128, dt=1e-3 s, f1=60 Hz, A=2.5
N  = 128;
dt = 1e-3;
f1 = 60.0;
A  = 2.5;
t  = (0:N-1)' * dt;
x  = A * cos(2*pi*f1*t);

X_cc   = ufft_cc_forward(x, dt);
x_rec  = ufft_cc_inverse(X_cc, dt);

write_complex_csv(fullfile(outdir, 'matlab_X_cc_cosine.csv'),    X_cc);
write_complex_csv(fullfile(outdir, 'matlab_x_cc_rec_cosine.csv'), x_rec);
fprintf('[MATLAB] cosine CC forward/inverse written.\n');

%% 2. LP filter on random signal: N=256, dt=60 s
N2   = 256;
dt2  = 60.0;
fc   = 1/3600;
rng(0);                         % seed matches Python rng default_rng(0) approximately
x2   = randn(N2, 1);
f2   = ufft_freqs(N2, dt2);
H2   = double(abs(f2) <= fc);

y_filtered = ufft_filter(x2, H2, dt2);
write_complex_csv(fullfile(outdir, 'matlab_y_filtered.csv'), y_filtered);
fprintf('[MATLAB] filter output written.\n');

%% 3. Impulse response: CC inverse of boxcar TF
H3 = double(abs(f2) <= fc);
h  = ufft_cc_inverse(complex(H3), dt2);
write_complex_csv(fullfile(outdir, 'matlab_h_impulse.csv'), h);
integral = sum(real(h)) * dt2;
fprintf('[MATLAB] impulse response written. Integral = %.8f (expect ~1.0)\n', integral);

%% helper
function write_complex_csv(fpath, z)
    fid = fopen(fpath, 'w');
    for k = 1:numel(z)
        r = real(z(k));
        im = imag(z(k));
        if im >= 0
            fprintf(fid, '%.17g+%.17gj\n', r, im);
        else
            fprintf(fid, '%.17g%.17gj\n', r, im);
        end
    end
    fclose(fid);
end
