!> Fortran cross-language validation demo.
!> Writes CSV files to ../tests/data/ matching the Python reference vectors.
!>
!> Compile:  make demo
!> Run:      ./ufft_demo_f [output_dir]
program ufft_demo
    use universalfft_mod
    use iso_fortran_env, only: dp => real64
    implicit none

    integer,  parameter :: N  = 128
    real(dp), parameter :: dt = 1.0e-3_dp
    real(dp), parameter :: f1 = 60.0_dp
    real(dp), parameter :: A  = 2.5_dp
    real(dp), parameter :: PI = 3.14159265358979323846_dp

    ! Input signal
    real(dp) :: sig_re(0:N-1), sig_im(0:N-1)
    ! Spectrum
    real(dp) :: spec_re(0:N-1), spec_im(0:N-1)
    ! Round-trip reconstructed signal
    real(dp) :: rec_re(0:N-1), rec_im(0:N-1)
    integer  :: k, rc

    character(len=256) :: outdir
    character(len=512) :: fpath

    if (command_argument_count() >= 1) then
        call get_command_argument(1, outdir)
    else
        outdir = "../tests/data"
    end if

    ! ── 1. Cosine signal ─────────────────────────────────────────────────────
    do k = 0, N-1
        sig_re(k) = A * cos(2.0_dp * PI * f1 * real(k, dp) * dt)
        sig_im(k) = 0.0_dp
    end do

    rc = ufft_cc_forward(sig_re, sig_im, spec_re, spec_im, N, dt)
    write(fpath,'(A,"/fortran_X_cc_cosine.csv")') trim(outdir)
    call write_complex_csv(trim(fpath), spec_re, spec_im, N)
    write(*,'(A,A)') "[Fortran] Wrote ", trim(fpath)

    rc = ufft_cc_inverse(spec_re, spec_im, rec_re, rec_im, N, dt)
    write(fpath,'(A,"/fortran_x_cc_rec_cosine.csv")') trim(outdir)
    call write_complex_csv(trim(fpath), rec_re, rec_im, N)
    write(*,'(A,A)') "[Fortran] Wrote ", trim(fpath)

    ! ── 2. Impulse response (N=256, dt=60 s) ─────────────────────────────────
    block
        integer,  parameter :: N2  = 256
        real(dp), parameter :: dt2 = 60.0_dp
        real(dp), parameter :: fc  = 1.0_dp / 3600.0_dp
        real(dp) :: f2(0:N2-1), H2r(0:N2-1), H2i(0:N2-1)
        real(dp) :: imp_re(0:N2-1), imp_im(0:N2-1)
        real(dp) :: integral

        call ufft_freqs(f2, N2, dt2)
        call ufft_low_pass(H2r, f2, N2, fc)
        H2i = 0.0_dp

        rc = ufft_cc_inverse(H2r, H2i, imp_re, imp_im, N2, dt2)
        write(fpath,'(A,"/fortran_h_impulse.csv")') trim(outdir)
        call write_complex_csv(trim(fpath), imp_re, imp_im, N2)
        write(*,'(A,A)') "[Fortran] Wrote ", trim(fpath)

        integral = sum(imp_re) * dt2
        write(*,'(A,F12.8,A)') "[Fortran] Impulse integral = ", integral, "  (expect ≈1.0)"
    end block

    write(*,*) "[Fortran] Demo complete."

contains

    subroutine write_complex_csv(path, zr, zi, M)
        character(len=*), intent(in) :: path
        integer,          intent(in) :: M
        real(dp),         intent(in) :: zr(0:M-1), zi(0:M-1)
        integer :: u, j
        open(newunit=u, file=path, status='replace', action='write')
        do j = 0, M-1
            if (zi(j) >= 0.0_dp) then
                write(u,'(ES24.17,"+",ES24.17,"j")') zr(j), zi(j)
            else
                write(u,'(ES24.17,ES24.17,"j")')     zr(j), zi(j)
            end if
        end do
        close(u)
    end subroutine

end program ufft_demo
