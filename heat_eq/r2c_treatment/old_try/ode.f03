module ode
    implicit none
    INTEGER :: i        ! For iteration
    
    contains
    function euler(F, f0, x0, x, n) result(fn)
        ! f0 = f(x0), n = # slices, x = desired point of evaluation, F = Function
        implicit none
        REAL :: F, f0, x0, x, fn, dx
        INTEGER :: n

        ! step
        dx = (x - x0) / REAL(n)

        do i = 1, n
            x = x0 + dx
            fn = f0 + dx * F(x0, f0)

            x0 = x
            f0 = fn
        end do
    end function euler

    function RK2(F, f0, x0, x, n) result(fn)
        implicit none
        REAL :: F, f0, x0, x, fn, dx
        INTEGER :: n
        REAL :: k1, k2

        ! step
        dx = (x - x0) / REAL(n)

        do i = 1, n
            x = x0 + dx

            k1 = dx * F(x0, f0)
            k2 = dx * F(x0 + 0.5*dx, f0 + 0.5*k1)

            fn = f0 + k2

            x0 = x
            f0 = fn
        end do
    end function RK2

    function RK4(F, f0, x0, x, n) result(fn)
        implicit none
        REAL :: F, f0, x0, x, fn, dx
        INTEGER :: n
        REAL :: k1, k2, k3, k4

        ! step
        dx = (x - x0) / REAL(n)

        do i = 1, n
            x = x0 + dx

            k1 = dx * F(x0, f0)
            k2 = dx * F(x0 + 0.5*dx, f0 + 0.5*k1)
            k3 = dx * F(x0 + 0.5*dx, f0 + 0.5*k2)
            k4 = dx * F(x0 + dx, f0 + k3)

            fn = f0 + (1.0/3.0) * (0.5*k1 + k2 + k3 + 0.5*k4)

            x0 = x
            f0 = fn
        end do
    end function RK4

end module ode