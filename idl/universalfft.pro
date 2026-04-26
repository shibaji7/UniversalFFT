;;; universalfft.pro — Boteler (2012) compliant FFT/IFFT for IDL / GDL.
;;;
;;; GDL (GNU Data Language) is the free open-source IDL-compatible runtime:
;;;   sudo apt install gnudatalanguage       ; Debian/Ubuntu
;;;   brew install gnudatalanguage           ; macOS
;;;   Run: gdl universalfft.pro
;;;
;;; IMPORTANT — IDL/GDL normalisation is the OPPOSITE of every other library:
;;;   IDL FFT(x, -1)  = (1/N) * rawForwardSum   (already divided by N!)
;;;   IDL FFT(X, +1)  = rawInverseSum            (no 1/N)
;;;
;;; Boteler mapping:
;;;   CC forward  = IDL_FFT(x,-1) * N * dt   (undo IDL's 1/N, then scale by dt)
;;;   CC inverse  = IDL_FFT(X,+1) * df        (df = 1/(N*dt))
;;;   CD forward  = IDL_FFT(x,-1)             (IDL's 1/N is exactly what CD needs)
;;;   CD inverse  = IDL_FFT(X,+1)             (raw sum — correct for CD)

;; ── CC forward  (Boteler Eq. 21a) ─────────────────────────────────────────
FUNCTION ufft_cc_forward, x, dt
  N = N_ELEMENTS(x)
  RETURN, FFT(x, -1) * (N * dt)
END

;; ── CC inverse  (Boteler Eq. 21b) ─────────────────────────────────────────
FUNCTION ufft_cc_inverse, X, dt
  N  = N_ELEMENTS(X)
  df = 1.0D / (N * dt)
  RETURN, FFT(X, +1) * df
END

;; ── CD forward  (Boteler Eq. 22a) ─────────────────────────────────────────
FUNCTION ufft_cd_forward, x
  ;; IDL FFT(x,-1) already divides by N — exactly the CD forward convention.
  RETURN, FFT(x, -1)
END

;; ── CD inverse  (Boteler Eq. 22b) ─────────────────────────────────────────
FUNCTION ufft_cd_inverse, X
  RETURN, FFT(X, +1)
END

;; ── Frequency bin array (Hz, FFT-output order) ────────────────────────────
FUNCTION ufft_freqs, N, dt
  df = 1.0D / (N * dt)
  k  = DINDGEN(N)
  f  = DBLARR(N)
  idx_pos = WHERE(k LT N/2, cnt_pos)
  idx_neg = WHERE(k GE N/2, cnt_neg)
  IF cnt_pos GT 0 THEN f[idx_pos] = k[idx_pos] * df
  IF cnt_neg GT 0 THEN f[idx_neg] = (k[idx_neg] - N) * df
  RETURN, f
END

;; ── Low-pass brick-wall response ──────────────────────────────────────────
FUNCTION ufft_low_pass, f, fc
  RETURN, COMPLEX(ABS(f) LE fc, 0)
END

;; ── Filter via CC pair ─────────────────────────────────────────────────────
FUNCTION ufft_filter, x, H, dt
  X = ufft_cc_forward(x, dt)
  RETURN, ufft_cc_inverse(X * H, dt)
END

;; ── Demo / test (runs when this file is executed directly) ─────────────────
PRO universalfft_demo, outdir=outdir
  IF N_ELEMENTS(outdir) EQ 0 THEN outdir = '../tests/data'

  ;; cosine: N=128, dt=1e-3, f1=60 Hz, A=2.5
  N = 128L & dt = 1D-3 & f1 = 60D & A = 2.5D
  t = DINDGEN(N) * dt
  x = A * COS(2D * !DPI * f1 * t)

  X_cc  = ufft_cc_forward(x, dt)
  x_rec = ufft_cc_inverse(X_cc, dt)

  max_err = MAX(ABS(REAL_PART(x_rec) - x))
  PRINT, '[IDL/GDL] CC round-trip max error: ', max_err, '  (expect < 1e-9)'

  ;; impulse response: N=256, dt=60 s
  N2 = 256L & dt2 = 60D & fc = 1D/3600D
  f2 = ufft_freqs(N2, dt2)
  H  = ufft_low_pass(f2, fc)
  h  = ufft_cc_inverse(H, dt2)

  integral = TOTAL(REAL_PART(h)) * dt2
  PRINT, '[IDL/GDL] impulse response integral = ', integral, '  (expect ≈ 1.0)'

  PRINT, '[IDL/GDL] Demo complete.'
END

universalfft_demo
END
