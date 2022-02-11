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

    INTEGER, PARAMETER :: Nx = 8, Ny = 8                              ! Nx = t space discretization, Ny = x space discretization
    REAL*8, PARAMETER :: PI = 4*ATAN(1.0)
    REAL*8, PARAMETER :: T = 2*PI, L = 2*PI
    REAL*8, PARAMETER :: dt = REAL(T)/REAL(Nx), dx = REAL(L)/REAL(Ny), kappa = 1        ! Check CFL stability Condition for dx/dt

    INTEGER :: i, j                                                 ! parameter for iteration
    REAL*8 :: k(Ny)

    ! ! These two are required for plotting only, but if we use Seaborn.headmap in python (since that plots matrix)
    ! REAL*8 :: x(Ny) = [(i*dx, i=0, 100, 1)]                        ! x space array
    ! INTEGER :: k(Ny) = [(i, i=0,100, 1)]                           ! k space array

    ! FFTW DataTypes
    TYPE(C_PTR) :: planFFT, planIFFT
    REAL(C_DOUBLE) :: theta(Nx, Ny) = 0.0                             ! the variable
    COMPLEX(C_DOUBLE_COMPLEX) :: theta_hat(Nx, Ny)               ! transformed variable   NOTE for err when = 0
    ! (for r2c transformation, first dimension of the complex data is chopped roughly in half)
    ! Not an issue here, since we'll be transforming 1D

    
    ! The initial distribution : theta(x,0) = Reservoir at 100 0 0 0 ... 0 0 33
    theta(:,1) = 100.0
    theta(:,Ny) = 33.0

    ! Initial Config
    OPEN(UNIT=9, FILE='start.dat', STATUS="REPLACE")
    do i=1, Nx
        WRITE (9,*) theta(i,:)
    end do

    do i = 1, Ny/2
        k(i) = (i-1) / L
    end do
    k(Ny/2 +1) = 0.0
    do i = 1, Ny/2
        k(Ny/2+1 +i) = -k(Ny/2+1 -i)
    end do

    
    ! TODO: Time for FFTW to transform initial condition in k space
    planFFT = fftw_plan_dft_r2c_1d(Ny, theta(1,:), theta_hat(1,:), FFTW_ESTIMATE)    
    ! Note that the dimension is reversed (this is due to the column-major order in FORTRAN)
    ! We don't need that here, since we're transforming 1D
    
    CALL fftw_execute_dft_r2c(planFFT, theta(1,:), theta_hat(1,:))

    ! Initial Transformation
    OPEN(UNIT=11, FILE='fftw_forward.dat', STATUS="REPLACE")
    do i=1, Nx
        WRITE (11,*) theta_hat(i,:)
    end do
    
    
    ! TODO: Now, we use the **Euler method** to iterate through the time to get the solution in k space
    ! The equation now: d/dt theta_hat = - k^2 theta_hat for each k value in k space
    do j = 1, Ny             ! Along k space
        do i = 1, Nx     ! Along time dimension
            theta_hat(i+1,j) = theta_hat(i,j) + dt * (-k(j)**2 * theta_hat(i,j) * kappa)            ! yn+1 = yn + dx*F(xn,yn)
        end do
    end do
    
    ! TODO: Now, it's time to get the real space solution by using the IFFTW
    planIFFT = fftw_plan_dft_c2r_1d(Ny, theta_hat(1,:), theta(1,:), FFTW_ESTIMATE)

    do i = 1, Nx
        CALL fftw_execute_dft_c2r(planIFFT, theta_hat(i,:), theta(i,:))
    end do
    
    ! Normalization
    theta = 1/REAL(Nx*Ny) * theta

    ! Got the result, now to export it
    OPEN(UNIT=10, FILE='heat_sol.dat', STATUS="REPLACE")
    do i=1, Nx
        WRITE (10,*) theta(i,:)
    end do
    

    ! Termination
    CALL fftw_destroy_plan(planFFT)
    CALL fftw_destroy_plan(planIFFT)
end program heat_eq


! Why not use the regular FTCS Method?
! Because there, the nabla^2 is found by finite difference method, which has a O(x^2) error, which is not here

! So, we would rather use RK Method here for better accuracy in time derivative (1st trial : Euler Method)