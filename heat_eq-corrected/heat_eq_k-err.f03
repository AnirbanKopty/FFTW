! TODO: To solve heat equation del/del t (theta) = kappa nabla^2 (theta)  (kappa = 1)
! Using Spectral Methods by FFTW

! We first go to the Fourier Space, thus nabla^2 = -k^2, then solve for t in that space and
! Finally, return back to real space to get the final solution

!? Here, we have to deal with 2D array, but the fourier transformation will be on the space dimension only which is 1D

program heat_eq
    USE, INTRINSIC :: iso_c_binding     ! We can use C syntaxes for fftw now
    ! USE :: ode

    implicit none
    INCLUDE '/home/anirbankopty/Softwares/FFTW/fftw-install/include/fftw3.f03'
    ! The file where the definitions are for Fortran to C

    INTEGER, PARAMETER :: T = 5, L = 5                              ! L = x space dim, T = t space dim
    REAL*8, PARAMETER :: dx = 0.01, dt = 0.01, PI = 4*ATAN(1.0), kappa = 1

    INTEGER :: k, i                                                 ! parameter for iteration

    ! ! These two are required for plotting only, but if we use Seaborn.headmap in python (since that plots matrix)
    ! REAL*8 :: x(L) = [(i*dx, i=0, 100, 1)]                        ! x space array
    ! INTEGER :: k(L) = [(i, i=0,100, 1)]                           ! k space array

    ! FFTW DataTypes
    TYPE(C_PTR) :: planFFT, planIFFT
    REAL(C_DOUBLE) :: theta(T, L) = 0.0                             ! the variable
    COMPLEX(C_DOUBLE_COMPLEX) :: theta_hat(T, L)              ! transformed variable     
    ! (for r2c transformation, first dimension of the complex data is chopped roughly in half)
    ! Not an issue here, since we'll be transforming 1D

    
    ! The initial distribution : theta(x,0) = Reservoir at 100 0 0 0 ... 0 0 33
    theta(:,1) = 100.0
    theta(:,L) = 33.0

    ! Initial Config
    OPEN(UNIT=9, FILE='start.dat', STATUS="REPLACE")
    do i=1, T
        WRITE (9,*) theta(i,:)
    end do

    
    ! TODO: Time for FFTW to transform initial condition in k space
    planFFT = fftw_plan_dft_r2c_1d(L, theta(1,:), theta_hat(1,:), FFTW_ESTIMATE)    
    ! Note that the dimension is reversed (this is due to the column-major order in FORTRAN)
    ! We don't need that here, since we're transforming 1D
    
    CALL fftw_execute_dft_r2c(planFFT, theta(1,:), theta_hat(1,:))

    ! Initial Transformation
    OPEN(UNIT=11, FILE='fftw_forward.dat', STATUS="REPLACE")
    do i=1, T
        WRITE (11,*) theta_hat(i,:)
    end do
    
    
    ! TODO: Now, we use the **Euler method** to iterate through the time to get the solution in k space
    ! The equation now: d/dt theta_hat = - k^2 theta_hat for each k value in k space
    do k = 1, L             ! Along k space
        do i = 1, T     ! Along time dimension
            theta_hat(i+1,k) = theta_hat(i,k) + dt * (-k**2 * theta_hat(i,k) * kappa)            ! yn+1 = yn + dx*F(xn,yn)
        end do
    end do
    
    ! TODO: Now, it's time to get the real space solution by using the IFFTW
    planIFFT = fftw_plan_dft_c2r_1d(L, theta_hat(1,:), theta(1,:), FFTW_ESTIMATE)

    do i = 1, T
        CALL fftw_execute_dft_c2r(planIFFT, theta_hat(i,:), theta(i,:))
    end do
    
    ! Normalization
    theta = 1/REAL(T*L) * theta

    ! Got the result, now to export it
    OPEN(UNIT=10, FILE='heat_sol.dat', STATUS="REPLACE")
    do i=1, T
        WRITE (10,*) theta(i,:)
    end do
    

    ! Termination
    CALL fftw_destroy_plan(planFFT)
    CALL fftw_destroy_plan(planIFFT)
end program heat_eq


! Why not use the regular FTCS Method?
! Because there, the nabla^2 is found by finite difference method, which has a O(x^2) error, which is not here

! So, we would rather use RK Method here for better accuracy in time derivative (1st trial : Euler Method)