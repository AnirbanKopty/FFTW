program testfreq
    use freq_mod
    implicit none
    INTEGER :: i
    INTEGER, PARAMETER :: N = 11
    REAL*8 :: k(N), L=1, d

    d = L/ REAL(N)

    k = fftw_freq(N, d)
    ! CALL fftw_freq_sub(k, N, d)

    !  printing the result
    do i = 1, N
        PRINT *, k(i)
    end do

end program testfreq