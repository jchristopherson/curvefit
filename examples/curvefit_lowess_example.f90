! curvefit_lowess_example.f90

program example
    use iso_fortran_env
    use curvefit_regression
    use curvefit_core
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: n = 100
    real(real64), parameter :: maxX = 1.0d0
    real(real64), parameter :: minX = 0.0d0
    type(color), parameter :: orange = color(255, 102, 0)

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(n), y(n), yr(n), ys(n), ys2(n), dx, cnl(5), ynl(n)
    type(lowess_smoothing) :: fit
    type(nonlinear_regression) :: solver
    procedure(reg_fcn), pointer :: fcn
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2, d3, d4, d5

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
    call solver%solve(cnl)
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

    ! Plot the data
    call plt%initialize()
    call plt%set_font_size(14)

    call d1%set_name("Original Signal")
    call d1%set_line_color(CLR_BLUE)
    call d1%set_line_width(2.0)
    call d1%define_data(x, y)

    call d2%set_name("Raw Data")
    call d2%set_line_color(CLR_BLACK)
    call d2%set_draw_line(.false.)
    call d2%set_draw_markers(.true.)
    call d2%set_marker_style(MARKER_EMPTY_CIRCLE)
    call d2%set_marker_scaling(2.0)
    call d2%set_line_width(2.0)
    call d2%define_data(x, yr)

    call d3%set_name("Smoothed (f = 0.2)")
    call d3%set_line_color(CLR_GREEN)
    call d3%set_line_width(2.0)
    call d3%set_line_style(LINE_DASHED)
    call d3%define_data(x, ys)

    call d4%set_name("Smoothed (f = 0.8)")
    call d4%set_line_color(CLR_MAGENTA)
    call d4%set_line_width(2.0)
    call d4%set_line_style(LINE_DASH_DOTTED)
    call d4%define_data(x, ys2)

    call d5%set_name("Nonlinear Regression")
    call d5%set_line_color(orange)
    call d5%set_line_width(2.0)
    call d5%set_line_style(LINE_DASH_DOT_DOT)
    call d5%define_data(x, ynl)

    call plt%push(d1)
    call plt%push(d2)
    call plt%push(d3)
    call plt%push(d4)
    call plt%push(d5)

    call plt%draw()

contains
    function nrfun(xp, c) result(fn)
        real(real64), intent(in) :: xp
        real(real64), intent(in), dimension(:) :: c
        real(real64) :: fn
        fn = c(1) * sin(c(2) * xp) + c(3) * cos(c(4) * xp) * exp(c(5) * xp)
    end function

end program
