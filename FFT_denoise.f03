! TODO : Using FFTW to denoise a noisy Data

program FFTW_denoise
    USE, INTRINSIC :: iso_c_binding     ! We can use C syntaxes for fftw now
    INCLUDE 'fftw3.f03'                 ! The file where the definitions are for Fortran to C
    
    implicit none
    REAL, PARAMETER :: dt = 0.001, PI = 4*ATAN(1.0)
    INTEGER :: i
    REAL :: t(1001) = (/(i*dt, i=0, 1000, 1)/)                 ! t array holding t axes
    REAL :: f(1001), f_clean(1001)
    REAL :: temp

    ! FFTW DataTypes
    TYPE(C_PTR) :: plan
    COMPLEX(C_DOUBLE_COMPLEX) :: in(1001), out(1001)

    ! clean signal
    f_clean = SIN(2*PI*50*t) + SIN(2*PI*150*t) + SIN(2*PI*200*t)

    ! for Plotting
    OPEN(UNIT=10, FILE='f_clean.dat', STATUS='REPLACE')
    OPEN(UNIT=11, FILE='f_noisy.dat', STATUS="REPLACE")


    do i=1, 1001
        ! Exporting f_clean
        WRITE(10,*) f_clean(i)

        ! Adding random noise
        CALL RANDOM_NUMBER(temp)
        f = f_clean + temp

        ! Exporting f_noisy
        WRITE(11,*) f(i)
    end do

    ! FFTW Section
    plan = fftw_plan_dft_1d




end program FFTW_denoise