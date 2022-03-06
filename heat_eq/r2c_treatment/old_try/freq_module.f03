module freq_mod
    contains
        function fftw_freq(N, d) result(k)
            implicit none
            INTEGER :: i
            INTEGER, INTENT(IN) :: N
            REAL*8 :: d, L, k(N)

            L = REAL(N)*d
            ! d = L / REAL(N)

            ! Source: https://stackoverflow.com/a/51232221/11348113
            ! do i = 1, N/2
            !     k(i) = (i-1) / L
            ! end do
            ! k(N/2 +1) = 0.0
            ! do i = 1, N/2 -1
            !     k(N/2+1 +i) = -k(N/2+1 -i)
            ! end do

            
            ! Source: Documentation of numpy.fft.fftfreq
            if ( MODULO(N,2) == 0) then
                do i = 1, N/2
                    k(i) = (i-1) / L
                end do
                k(N/2 +1) = - N/2
                do i = 1, N/2 -1
                    k(N/2+1 +i) = -k(N/2+1 -i)
                end do

            else
                do i = 1, (N+1)/2
                    k(i) = (i-1) / L
                end do
                do i = 1, (N+1)/2 -1
                    k((N+1)/2 +i) = -k((N+1)/2 -i+1)
                end do

            end if


        end function fftw_freq

end module