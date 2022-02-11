program testfreq
    implicit none
    INTEGER :: i
    INTEGER, PARAMETER :: N = 8
    REAL :: k(N), L=1

    do i = 1, N/2
        k(i) = (i-1) / L
    end do
    k(N/2 +1) = 0.0
    do i = 1, N/2 -1
        k(N/2+1 +i) = -k(N/2+1 -i)
    end do

    do i = 1, N
        PRINT *, k(i)
    end do
end program testfreq