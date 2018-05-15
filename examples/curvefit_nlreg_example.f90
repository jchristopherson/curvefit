! curvefit_nlreg_example.f90

program example
    use iso_fortran_env
    use curvefit_core
    use curvefit_regression
    implicit none

    ! Local Variables
    type(nonlinear_regression) :: solver
    procedure(reg_fcn), pointer :: fcn
    real(real64) :: xp(21), yp(21), cp(4), yf(21)
    integer(int32) :: i, id

    ! Define the data to fit (21 data points)
    xp = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
        0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
        1.8d0, 1.9d0, 2.0d0]
    yp = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
        1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
        1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
        3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
        5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
        10.30252404d0]

    ! Define an initial estimate of the coefficients
    cp = 1.0d0 ! Set all coefficients to an initial guess of 1

    ! Set up the solver
    fcn => nrfun
    call solver%initialize(xp, yp, fcn, 4)

    ! Compute the coefficients
    call solver%solve(cp, yf)

    ! Display the coefficients
    print '(AF12.10)', "c0: ", cp(4)
    print '(AF12.10)', "c1: ", cp(3)
    print '(AF12.10)', "c2: ", cp(2)
    print '(AF12.10)', "c3: ", cp(1)
    print '(AF7.5)', "Max Residual: ", maxval(abs(yf))

    ! Write the results to file for plotting purposes
    open(newunit = id, file = "nl_regress.txt", action = "write", &
        status = "replace")
    do i = 1, size(xp)
        write(id, '(F14.10AF14.10AF14.10)') xp(i), ",", yp(i), ",", yf(i)
    end do
    close(id)

contains
    ! The function to fit
    function nrfun(x, c) result(f)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: c
        real(real64) :: f
        f = c(1) * x**3 + c(2) * x**2 + c(3) * x + c(4)
    end function
end program
