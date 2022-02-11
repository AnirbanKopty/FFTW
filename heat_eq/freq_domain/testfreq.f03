program testfreq
    implicit none
    INTEGER :: i
    INTEGER, PARAMETER :: N = 9
    REAL :: kx(8), L=1

    do i = 1, N/2
        kx(i) = (i-1) / L
    end do
    kx(N/2 +1) = 0.0
    do i = 1, N/2
        kx(N/2+1 +i) = -kx(N/2+1 -i)
    end do

    do i = 1, N
        PRINT *, kx(i)
    end do
end program testfreq