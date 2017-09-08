! curvefit_lowess_example.f90

program example
    use curvefit_core
    use curvefit_regression
    implicit none

    ! Parameters
    integer(i32), parameter :: n = 100
    real(dp), parameter :: maxX = 1.0d0
    real(dp), parameter :: minX = 0.0d0

    ! Local Variables
    integer(i32) :: i, id
    real(dp) :: x(n), y(n), yr(n), ys(n), ys2(n), dx
    type(lowess_smoothing) :: fit

    ! Initialization
    dx = (maxX - minX) / (n - 1.0d0)
    x(1) = minX
    do i = 2, n
        x(i) = x(i-1) + dx
    end do
    y = 0.5d0 * sin(2.0d1 * x) + cos(5.0d0 * x) * exp(-0.1d0 * x)
    call random_number(yr)
    yr = y + (yr - 0.5d0)

    ! Generate the fit
    call fit%initialize(x, yr)
    ys = fit%smooth(0.2d0)
    ys2 = fit%smooth(0.8d0)

    ! Write the results to a text file
    open(newunit = id, file = "lowess.txt", action = "write", &
        status = "replace")
    do i = 1, n
        write(id, '(F14.10AF14.10AF14.10AF14.10AF14.10)') x(i), ",", y(i), &
            ",", yr(i), ",", ys(i), ",", ys2(i)
    end do
    close(id)
end program
