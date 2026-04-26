!> UniversalFFT — Boteler (2012) compliant FFT/IFFT in Fortran 90.
!>
!> Self-contained Cooley-Tukey radix-2 DIT FFT — no external dependencies.
!> N must be a power of two.
!>
!> Convention (same as C implementation):
!>   CC forward : X(k) = Σ x(n) e^{-i2πkn/N} Δt
!>   CC inverse : x(n) = Σ X(k) e^{+i2πkn/N} Δf   (Δf = 1/(N·Δt))
!>   CD forward : X(k) = (1/N) Σ x(n) e^{-i2πkn/N}
!>   CD inverse : x(n) = Σ X(k) e^{+i2πkn/N}       (raw sum)
module universalfft_mod
    use iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: ufft_inplace, ufft_cc_forward, ufft_cc_inverse, &
              ufft_cd_forward, ufft_cd_inverse, ufft_freqs,   &
              ufft_filter, ufft_low_pass

    real(dp), parameter :: PI = 3.14159265358979323846_dp

contains

    ! ── Bit-reversal permutation ─────────────────────────────────────────────
    subroutine bit_reverse(xr, xi, N)
        integer,    intent(in)    :: N
        real(dp),   intent(inout) :: xr(0:N-1), xi(0:N-1)
        integer :: i, j, bit
        real(dp) :: tmp
        j = 0
        do i = 1, N-1
            bit = N / 2
            do while (iand(j, bit) /= 0)
                j = ieor(j, bit)
                bit = bit / 2
            end do
            j = ieor(j, bit)
            if (i < j) then
                tmp = xr(i); xr(i) = xr(j); xr(j) = tmp
                tmp = xi(i); xi(i) = xi(j); xi(j) = tmp
            end if
        end do
    end subroutine

    ! ── Core in-place Cooley-Tukey FFT ───────────────────────────────────────
    !> inverse=0 → forward (e^{-i2πkn/N}), inverse=1 → backward (e^{+i2πkn/N})
    !> Returns 0 on success, -1 if N is not a power of two.
    function ufft_inplace(xr, xi, N, inverse) result(rc)
        integer,    intent(in)    :: N, inverse
        real(dp),   intent(inout) :: xr(0:N-1), xi(0:N-1)
        integer :: rc, len, i, j
        real(dp) :: ang, wr, wi, wre, wim, tre, tim, wre_new
        ! check power of two
        if (N <= 0 .or. iand(N, N-1) /= 0) then
            rc = -1; return
        end if
        call bit_reverse(xr, xi, N)
        len = 2
        do while (len <= N)
            ang = merge(-1.0_dp, 1.0_dp, inverse == 0) * 2.0_dp * PI / real(len, dp)
            wr = cos(ang); wi = sin(ang)
            i = 0
            do while (i < N)
                wre = 1.0_dp; wim = 0.0_dp
                do j = 0, len/2 - 1
                    tre = wre * xr(i+j+len/2) - wim * xi(i+j+len/2)
                    tim = wre * xi(i+j+len/2) + wim * xr(i+j+len/2)
                    xr(i+j+len/2) = xr(i+j) - tre
                    xi(i+j+len/2) = xi(i+j) - tim
                    xr(i+j) = xr(i+j) + tre
                    xi(i+j) = xi(i+j) + tim
                    wre_new = wre*wr - wim*wi
                    wim     = wre*wi + wim*wr
                    wre     = wre_new
                end do
                i = i + len
            end do
            len = len * 2
        end do
        rc = 0
    end function

    ! ── CC forward ────────────────────────────────────────────────────────────
    function ufft_cc_forward(in_r, in_i, out_r, out_i, N, dt) result(rc)
        integer,  intent(in)  :: N
        real(dp), intent(in)  :: in_r(0:N-1), in_i(0:N-1), dt
        real(dp), intent(out) :: out_r(0:N-1), out_i(0:N-1)
        integer :: rc
        out_r = in_r; out_i = in_i
        rc = ufft_inplace(out_r, out_i, N, 0)
        if (rc == 0) then
            out_r = out_r * dt
            out_i = out_i * dt
        end if
    end function

    ! ── CC inverse ────────────────────────────────────────────────────────────
    function ufft_cc_inverse(in_r, in_i, out_r, out_i, N, dt) result(rc)
        integer,  intent(in)  :: N
        real(dp), intent(in)  :: in_r(0:N-1), in_i(0:N-1), dt
        real(dp), intent(out) :: out_r(0:N-1), out_i(0:N-1)
        integer :: rc
        real(dp) :: df
        df = 1.0_dp / (real(N, dp) * dt)
        out_r = in_r; out_i = in_i
        rc = ufft_inplace(out_r, out_i, N, 1)   ! inverse twiddles
        if (rc == 0) then
            out_r = out_r * df   ! raw sum * Δf  (ufft_inplace has no 1/N)
            out_i = out_i * df
        end if
    end function

    ! ── CD forward ────────────────────────────────────────────────────────────
    function ufft_cd_forward(in_r, in_i, out_r, out_i, N) result(rc)
        integer,  intent(in)  :: N
        real(dp), intent(in)  :: in_r(0:N-1), in_i(0:N-1)
        real(dp), intent(out) :: out_r(0:N-1), out_i(0:N-1)
        integer :: rc
        out_r = in_r; out_i = in_i
        rc = ufft_inplace(out_r, out_i, N, 0)
        if (rc == 0) then
            out_r = out_r / real(N, dp)
            out_i = out_i / real(N, dp)
        end if
    end function

    ! ── CD inverse ────────────────────────────────────────────────────────────
    function ufft_cd_inverse(in_r, in_i, out_r, out_i, N) result(rc)
        integer,  intent(in)  :: N
        real(dp), intent(in)  :: in_r(0:N-1), in_i(0:N-1)
        real(dp), intent(out) :: out_r(0:N-1), out_i(0:N-1)
        integer :: rc
        out_r = in_r; out_i = in_i
        rc = ufft_inplace(out_r, out_i, N, 1)   ! raw sum — no extra factor
    end function

    ! ── Frequency bin array ──────────────────────────────────────────────────
    subroutine ufft_freqs(f, N, dt)
        integer,  intent(in)  :: N
        real(dp), intent(in)  :: dt
        real(dp), intent(out) :: f(0:N-1)
        real(dp) :: df
        integer  :: k
        df = 1.0_dp / (real(N, dp) * dt)
        do k = 0, N-1
            if (k < N/2) then
                f(k) = real(k, dp) * df
            else
                f(k) = real(k - N, dp) * df
            end if
        end do
    end subroutine

    ! ── Low-pass brick-wall response ─────────────────────────────────────────
    subroutine ufft_low_pass(H, f, N, fc)
        integer,  intent(in)  :: N
        real(dp), intent(in)  :: f(0:N-1), fc
        real(dp), intent(out) :: H(0:N-1)
        integer :: k
        do k = 0, N-1
            H(k) = merge(1.0_dp, 0.0_dp, abs(f(k)) <= fc)
        end do
    end subroutine

    ! ── Filter via CC pair ────────────────────────────────────────────────────
    function ufft_filter(in_r, in_i, H_r, H_i, out_r, out_i, N, dt) result(rc)
        integer,  intent(in)  :: N
        real(dp), intent(in)  :: in_r(0:N-1), in_i(0:N-1)
        real(dp), intent(in)  :: H_r(0:N-1),  H_i(0:N-1)
        real(dp), intent(in)  :: dt
        real(dp), intent(out) :: out_r(0:N-1), out_i(0:N-1)
        integer :: rc, k
        real(dp) :: Xr(0:N-1), Xi(0:N-1), yr, yi
        rc = ufft_cc_forward(in_r, in_i, Xr, Xi, N, dt)
        if (rc /= 0) return
        do k = 0, N-1
            yr = Xr(k)*H_r(k) - Xi(k)*H_i(k)
            yi = Xr(k)*H_i(k) + Xi(k)*H_r(k)
            Xr(k) = yr; Xi(k) = yi
        end do
        rc = ufft_cc_inverse(Xr, Xi, out_r, out_i, N, dt)
    end function

end module universalfft_mod
