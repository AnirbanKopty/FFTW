! Using Spectral Methods by FFTW



! We first go to the Fourier Space, thus nabla^2 = -k^2, then solve for t in that space and
! Finally, return back to real space to get the final solution

!? Here, we have to deal with 2D array

program heat_eq
    USE, INTRINSIC :: iso_c_binding     ! We can use C syntaxes for fftw now
    ! USE :: ode

    implicit none
    INCLUDE '/home/anirbankopty/Softwares/FFTW/fftw-install/include/fftw3.f03'
    ! The file where the definitions are for Fortran to C

    ! L = spatial length, Nx = time length, N = # discretization points
    REAL*8, PARAMETER :: PI = 4*ATAN(1.0)
    INTEGER, PARAMETER :: Nx = 8, Ny = 8              ! Nx = time discretization, Ny = Space discretization (Needs to be of power 2 for exact representation in binary to make no error)
    REAL*8, PARAMETER :: L=1 !!!!2*PI to test                         ! Length scale is 2PI
    REAL*8, PARAMETER :: dx = REAL(L)/REAL(Nx)
    REAL*8, PARAMETER :: dt = 0.001, kappa = 1    !kappa = diffusion coefficient        !! CFL COndition, **De-alising (aliasing error)

    INTEGER :: i, j                                                 ! parameter for iteration
    REAL*8 :: k(Ny)! = [(i*2*PI/REAL(L), i=-Ny/2, Ny/2-1, 1)]                                           ! k in k space = 2pi/L * i for i=0,N

    ! FFTW DataTypes
    TYPE(C_PTR) :: planFFT, planIFFT
    COMPLEX(C_DOUBLE_COMPLEX) :: theta(Nx,Ny) = 0, theta_hat(Nx,Ny) = 0
    
    ! The initial distribution : theta(x,0) = Reservoir at 100 0 0 0 ... 0 0 33
    theta(:,1) = 100.0
    theta(:,Ny) = 33.0
    
    ! Initial config
    OPEN(UNIT=9, FILE='start_cmplx.dat', STATUS="REPLACE")
    do i=1, Nx
        WRITE (9,*) theta(i,:)
    end do

    !? C Code for calculating frequencies
    !     //Frequencies 
    !     for (int i = 0; i< N / 2; ++i)
    !         kx[i] = i / L;
    !     kx[N / 2] = 0.00;
    !     for (int i = 0; i < ((N / 2) - 1); ++i)
    !         kx[i + 1 + N / 2] = -kx[N / 2 - i - 1];
    
    ! Calculating frequency domain
    do i = 1, Ny/2
        k(i) = (i-1) / L
    end do
    k(Ny/2 +1) = 0.0
    do i = 1, Ny/2
        k(Ny/2+1 +i) = -k(Ny/2+1 -i)
    end do

    do i = 1, Ny
        PRINT *, k(i)
    end do

    
    ! ! TODO: Time for FFTW to transform initial condition in k space
    ! planFFT = fftw_plan_dft_2d(L,Nx, theta, theta_hat, FFTW_FORWARD, FFTW_ESTIMATE)    
    ! ! Note that the dimension is reversed (this is due to the column-major order in FORTRAN)
    
    ! CALL fftw_execute_dft(planFFT, theta, theta_hat)
    ! One can make it more efficient by transforming only the t=0 row ---- let's do this
    planFFT = fftw_plan_dft_1d(Nx, theta(1,:), theta_hat(1,:), FFTW_FORWARD, FFTW_ESTIMATE)

    CALL fftw_execute_dft(planFFT, theta(1,:), theta_hat(1,:))

    ! Test
    OPEN(UNIT=11, FILE='fftw_forward_cmplx.dat', STATUS="REPLACE")
    do i=1, Nx
        WRITE (11,*) theta_hat(i,:)
    end do
    
    
    ! TODO: Now, we use the **Euler method** to iterate through the time to get the solution in k space
    ! The equation now: d/dt theta_hat = - k^2 theta_hat for each k value in k space
    do j = 1, Ny             ! Along k space
        do i = 1, Nx-1         ! Along time dimension
            theta_hat(i+1,j) = theta_hat(i,j) + dt * (-k(j)**2 * theta_hat(i,j) * kappa)            ! yn+1 = yn + dx*F(xn,yn)
            !! The issue is with k values
        end do
    end do
    
    ! TODO: Now, it's time to get the real space solution by using the IFFTW
    planIFFT = fftw_plan_dft_1d(Nx, theta_hat(1,:), theta(1,:), FFTW_BACKWARD, FFTW_ESTIMATE)
    
    do i = 1, Nx
        CALL fftw_execute_dft(planIFFT, theta_hat(i,:), theta(i,:))
    end do

    !! Normalization
    theta = 1.0/(Nx*Ny) * theta
    
    ! Got the result, now to export it
    OPEN(UNIT=10, FILE='heat_sol_cmplx.dat', STATUS="REPLACE")
    do i=1, Nx
        WRITE (10,*) theta(i,:)
    end do

    OPEN(UNIT=11, FILE='heat_sol_cmplx_ps.dat', STATUS="REPLACE")
    do i=1, Nx
        WRITE (11,*) SQRT(REAL(theta(i,:))**2 + AIMAG(theta(i,:))**2)
    end do
    

    ! Termination
    CALL fftw_destroy_plan(planFFT)
    CALL fftw_destroy_plan(planIFFT)
end program heat_eq


! Why not use the regular FTCS Method?
! Because there, the nabla^2 is found by finite difference method, which has a O(x^2) error, which is not here

! So, we would rather use RK Method here for better accuracy in time derivative (1st trial : Euler Method)

!? C Code for calculating frequencies
!     //Frequencies 
!     for (int i = 0; i< N / 2; ++i) {
!         kx[i] = i / L;
!     }

!     kx[N / 2] = 0.00;

!     for (int i = 0; i < ((N / 2) - 1); ++i) {
!         kx[i + 1 + N / 2] = -kx[N / 2 - i - 1];
!     }