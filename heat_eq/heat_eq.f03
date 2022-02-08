! TODO: To solve heat equation del/del t (theta) = kappa nabla^2 theta  (kappa = 1)
! Using Spectral Methods by FFTW

! We first go to the Fourier Space, thus nabla^2 = -k^2, then solve for t in that space and
! Finally, return back to real space to get the final solution

program heat_eq
    use ode
    implicit none
    

    
end program heat_eq


! Why not use the regular FTCS Method?
! Because there, the nabla^2 is found by finite difference method, which has a O(x^2) error, which is not here

! So, we would rather use RK Method here for better accuracy in time derivative