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
    real(dp) :: x(n), y(n), yr(n), ys(n), ys2(n), dx, cnl(5), ynl(n), res(n)
    type(lowess_smoothing) :: fit
    type(nonlinear_regression) :: solver
    procedure(reg_fcn), pointer :: fcn

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

    ! For comparison purposes, consider a nonlinear regression fit.  As we know
    ! the coefficients, they provide a very good starting guess.
    cnl = [0.5d0, 2.0d0, 20.0d0, 5.0d0, -0.1d0]
    fcn => nrfun
    call solver%initialize(x, yr, fcn, size(cnl))
    call solver%solve(cnl, res)
    do i = 1, n
        ynl(i) = fcn(x(i), cnl)
    end do

    ! Display the computed coefficients
    print '(A)', "f(x) = c0 * sin(c1 * x) + c2 * cos(c3 * x) * exp(c4 * x):"
    print '(AF12.10)', "c0: ", cnl(1)
    print '(AF13.10)', "c1: ", cnl(2)
    print '(AF12.10)', "c2: ", cnl(3)
    print '(AF12.10)', "c3: ", cnl(4)
    print '(AF13.10)', "c4: ", cnl(5)

    ! Write the results to a text file
    open(newunit = id, file = "lowess.txt", action = "write", &
        status = "replace")
    do i = 1, n
        write(id, '(F14.10AF14.10AF14.10AF14.10AF14.10AF14.10)') x(i), ",", &
            y(i), ",", yr(i), ",", ys(i), ",", ys2(i), ",", ynl(i)
    end do
    close(id)

contains
    function nrfun(xp, c) result(fn)
        real(dp), intent(in) :: xp
        real(dp), intent(in), dimension(:) :: c
        real(dp) :: fn
        fn = c(1) * sin(c(2) * xp) + c(3) * cos(c(4) * xp) * exp(c(5) * xp)
    end function

end program
