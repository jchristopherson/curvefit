! curvefit_nonlin_reg_example.f90

program example
    use iso_fortran_env
    use curvefit_regression
    use curvefit_core
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 100
    real(real64), parameter :: xMin = 0.0d0
    real(real64), parameter :: xMax = 1.0d0
    real(real64), parameter :: a = 2.0d0
    real(real64), parameter :: w = 30.0d0
    real(real64), parameter :: d = 0.2d0
    real(real64), parameter :: p = 0.0d0

    ! Local variables
    procedure(reg_fcn), pointer :: fun
    real(real64) :: x(npts), y(npts), c(4), yfit(npts)
    type(nonlinear_regression) :: solver
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2
    integer(int32) :: i

    ! Define the data to fit
    x = linspace(xMin, xMax, npts)
    call random_number(y)
    y = a * exp(-w * d * x) * sin(w * x - p) + (y - 0.5d0)

    ! Set up the solver.  Assume an equation of the form:
    !   A exp(-w d x) sin(w * x - p)
    fun => fcn
    c = [1.0d0, 20.0d0, 0.0d0, 0.0d0]
    call solver%initialize(x, y, fun, size(c))

    ! Solve
    call solver%solve(c)

    ! Compute the fitted form of the routine
    do i = 1, size(x)
        yfit(i) = fun(x(i), c)
    end do

    ! Print out the coefficients
    print '(AF8.4AF8.4)', "A = ", c(1), ", Actual = ", a
    print '(AF8.4AF8.4)', "w = ", c(2), ", Actual = ", w
    print '(AF8.4AF8.4)', "d = ", c(3), ", Actual = ", d
    print '(AF8.4AF8.4)', "p = ", c(4), ", Actual = ", p

    ! Plot the results
    call plt%initialize()
    call plt%set_font_size(14)

    call d1%set_name("Raw Data")
    call d1%set_line_color(CLR_BLACK)
    call d1%set_draw_line(.false.)
    call d1%set_draw_markers(.true.)
    call d1%set_marker_style(MARKER_X)
    call d1%define_data(x, y)

    call d2%set_name("Fit")
    call d2%set_line_color(CLR_BLUE)
    call d2%set_line_width(2.0)
    call d2%define_data(x, yfit)

    call plt%push(d1)
    call plt%push(d2)
    call plt%draw()

contains
    function fcn(x, c) result(f)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: c
        real(real64) :: f

        f = c(1) * exp(-c(2) * c(3) * x) * sin(c(2) * x - c(4))
    end function
end program
