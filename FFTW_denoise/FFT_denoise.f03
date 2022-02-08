! TODO : Using FFTW to denoise a noisy Data

program FFTW_denoise
    USE, INTRINSIC :: iso_c_binding     ! We can use C syntaxes for fftw now

    implicit none
    INCLUDE '/home/anirbankopty/Softwares/FFTW/fftw-install/include/fftw3.f03'
    ! The file where the definitions are for Fortran to C

    REAL*8, PARAMETER :: dt = 0.001, PI = 4*ATAN(1.0)
    INTEGER, PARAMETER :: N = 1001
    INTEGER :: i
    INTEGER :: freq(1001) = [(i, i=0,1000, 1)]
    REAL*8 :: t(N) = [(i*dt, i=0, 1000, 1)]                 ! t array holding t axes
    REAL*8 :: temp

    ! FFTW DataTypes
    TYPE(C_PTR) :: planFFT, planIFFT
    REAL(C_DOUBLE) :: f(N), f_clean(N)
    COMPLEX(C_DOUBLE_COMPLEX) :: f_hat(N)

    ! clean signal
    f_clean = SIN(2*PI*50*t) + SIN(2*PI*150*t) + SIN(2*PI*200*t)

    ! for Plotting
    OPEN(UNIT=10, FILE='f_clean.dat', STATUS='REPLACE')
    OPEN(UNIT=11, FILE='f_noisy.dat', STATUS="REPLACE")
    OPEN(UNIT=12, FILE='f_hat.dat', STATUS="REPLACE")
    OPEN(UNIT=13, FILE='f_clean_hat.dat', STATUS="REPLACE")
    OPEN(UNIT=14, FILE='f_denoised.dat', STATUS="REPLACE")


    do i=1, N
        ! Exporting f_clean
        WRITE(10,*) t(i), f_clean(i)

        ! Adding random noise
        CALL RANDOM_NUMBER(temp)
        ! Transforming the limit from [0,1) to [-1,1)
        temp = -1 + 2*temp
        f = f_clean + 3*temp        ! 3 times the noise

        ! Exporting f_noisy
        WRITE(11,*) t(i), f(i)
    end do

    ! FFTW Section
    planFFT = fftw_plan_dft_r2c_1d(N, f, f_hat, FFTW_ESTIMATE)

    CALL fftw_execute_dft_r2c(planFFT, f, f_hat)

    ! Plotting the Fourier Transformed Space
    ! real(), aimag()
    do i=1, N
        WRITE(12,*) freq(i), 1.0/REAL(N) * ( REAL(f_hat(i))**2 + AIMAG(f_hat(i))**2 )      ! Magnitude
    end do

    ! We only select 50 to 250 region 
    !(from the plot, there's a spike in 0the freq, since the random dist is not gaussian with infinite freq series)
    do i=1, N
        if ((i-1) < 50 .OR. (i-1) > 250) then
            f_hat(i) = 0.0 ! same as (0.0, 0.0)
        end if
    end do

    ! To plot
    do i=1, N
        WRITE(13,*) freq(i), 1.0/REAL(N) * ( REAL(f_hat(i))**2 + AIMAG(f_hat(i))**2 )      ! Magnitude
    end do

    ! IFFTW Section
    planIFFT = fftw_plan_dft_c2r_1d(N, f_hat, f, FFTW_ESTIMATE)

    CALL fftw_execute_dft_c2r(planIFFT, f_hat, f)

    ! The IFFTW doesn't scale down the inverse transform
    f = 1.0/REAL(N) * f

    ! Plotting back the real space
    do i=1, N
        WRITE(14,*) t(i), f(i)
    end do

    ! Termination Section
    CALL fftw_destroy_plan(planFFT)
    CALL fftw_destroy_plan(planIFFT)
end program FFTW_denoise