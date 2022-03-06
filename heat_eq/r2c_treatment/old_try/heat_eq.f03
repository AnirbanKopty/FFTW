! TODO: To solve 1D heat equation del/del t (theta) = kappa nabla^2 (theta)  (kappa = 1)
! Using Spectral Methods by FFTW

! We first go to the Fourier Space, thus nabla^2 = -k^2, then solve for t in that space and
! Finally, return back to real space to get the final solution

!? Here, we will be dealing with 2D array, but the fourier transformation will be on the space dimension only which is 1D

program heat_eq
    !!!!!!!!!! DECLARATION
    USE, INTRINSIC :: iso_c_binding     ! We can use C syntaxes for fftw now
    ! USE :: ode                        ! My own module for solving odes

    implicit none
    INCLUDE '/home/anirbankopty/Softwares/FFTW/fftw-install/include/fftw3.f03'
    ! The file where the definitions are there for Fortran to C

    REAL*8, PARAMETER :: PI = 4*ATAN(1.0)
    INTEGER, PARAMETER :: Nt = 250, Nx = 16                             ! Nt = t space discretization, Nx = x space discretization
    !!!!!! For some Nt, Nx values the solution diverges!! Like taking (100, 16) like Nx > 2^3
    REAL*8, PARAMETER :: L = 2*PI
    REAL*8, PARAMETER :: dt = 0.01, dx = REAL(L)/REAL(Nx), kappa = 1        !! Check CFL stability Condition for dt/dx : kappa*dt/dx <= 1 << Checked
    ! If CFL condition is not met, then it's found that the solution is not stable, it's fluctuating in the begnining, then reaches stability

    INTEGER :: i, j                                                  ! parameter for iteration
    REAL*8 :: k(Nx)                                                  ! k space array

    ! FFTW DataTypes
    TYPE(C_PTR) :: planFFT, planIFFT
    REAL(C_DOUBLE) :: theta(Nt, Nx) = 0.0                            ! the variable
    COMPLEX(C_DOUBLE_COMPLEX) :: theta_hat(Nt, Nx) = 0.0             ! transformed variable
    ! (for r2c transformation in multi-dimension, first dimension of the complex data is chopped roughly in half)
    ! Not an issue here, since we'll be transforming 1D

    
    !!!!!!!!!! Initial Condition
    ! The initial distribution : theta(x,0) = 100 0 0 0 ... 0 0 33
    theta(1,1) = 100.0
    theta(1,Nx) = 33.0

    ! Initial Config
    OPEN(UNIT=10, FILE='init_cond_x_space.dat', STATUS="REPLACE")
    do i=1, Nt
        WRITE (10,*) theta(i,:)
    end do

    !!? Source: https://stackoverflow.com/a/51232221/11348113
    ! k space array declaration - k = (2PI/L)*[0,1,2,...,0,...,-2,-1,-0]
    ! do i = 1, Nx/2
    !     k(i) = (i-1) * 2*PI / L
    ! end do
    ! k(Nx/2 +1) = 0.0
    ! do i = 1, Nx/2 -1
    !     k(Nx/2+1 +i) = -k(Nx/2+1 -i)
    ! end do

    !!? Source: numpy.fft.fft.freq
    if ( MODULO(Nx,2) == 0) then        ! even N case
        do i = 1, Nx/2
            k(i) = (i-1) * 2*PI/L
        end do
        k(Nx/2 +1) = - Nx/2
        do i = 1, Nx/2 -1
            k(Nx/2+1 +i) = -k(Nx/2+1 -i)
        end do
    else                                ! odd N case
        do i = 1, (Nx+1)/2
            k(i) = (i-1) * 2*PI/L
        end do
        do i = 1, (Nx+1)/2 -1
            k((Nx+1)/2 +i) = -k((Nx+1)/2 -i+1)
        end do
    end if

    ! !!! k space array declaration - k = (2PI/L)*[0,1,2,...,,...,-2,-1]
    ! k = 2*PI* fftw_freq(Nx, dx)
    ! !  printing the result
    do i = 1, Nx
        PRINT *, k(i)
    end do

    !!!!!!!!!!!! Transformation of Initial Condition
    ! TODO: Time for FFTW to transform initial condition in k space
    planFFT = fftw_plan_dft_r2c_1d(Nx, theta(1,:), theta_hat(1,:), FFTW_ESTIMATE)    
    ! Note that the dimension is reversed (this is due to the column-major order in FORTRAN)
    ! We don't need that here, since we're transforming 1D
    
    CALL fftw_execute_dft_r2c(planFFT, theta(1,:), theta_hat(1,:))

    ! Initial Condition Transformation
    OPEN(UNIT=11, FILE='init_cond_k_space.dat', STATUS="REPLACE")
    do i=1, Nt
        WRITE (11,*) theta_hat(i,:)
    end do
    
    
    !!!!!!!!!!!!! Solving ODE in k space to get time evolution
    ! TODO: Now, we use the **Euler method** to iterate through the time to get the solution in k space
    ! The equation now: d/dt theta_hat = - k^2 theta_hat for each k value in k space
    do j = 1, Nx             ! Along k space
        do i = 1, Nt-1       ! Along time dimension
            theta_hat(i+1,j) = theta_hat(i,j) + dt * (-k(j)**2 * theta_hat(i,j) * kappa)            ! yn+1 = yn + dx*F(xn,yn)
        end do
    end do

    ! After ODE solve in k space
    OPEN(UNIT=12, FILE='ode_solved_k_space.dat', STATUS="REPLACE")
    do i=1, Nt
        WRITE (12,*) theta_hat(i,:)
    end do
    

    !!!!!!!!!!!! Inverse FFTW to get to the Real Space Solution
    ! TODO: Now, it's time to get the real space solution by using the IFFTW
    planIFFT = fftw_plan_dft_c2r_1d(Nx, theta_hat(1,:), theta(1,:), FFTW_ESTIMATE)

    do i = 1, Nt
        CALL fftw_execute_dft_c2r(planIFFT, theta_hat(i,:), theta(i,:))
    end do
    
    ! Normalization
    theta = 1/REAL(Nx) * theta

    ! Got the result, now to export it
    OPEN(UNIT=13, FILE='heat_sol.dat', STATUS="REPLACE")
    do i=1, Nt
        WRITE (13,*) theta(i,:)
    end do
    

    ! Termination
    CALL fftw_destroy_plan(planFFT)
    CALL fftw_destroy_plan(planIFFT)
end program heat_eq


! Why not use the regular FTCS Method?
! Because there, the nabla^2 is found by finite difference method, which has a O(x^2) error, which is not here

! So, we would rather use RK Method here for better accuracy in time derivative (1st trial : Euler Method)