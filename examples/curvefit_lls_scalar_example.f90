! curvefit_lls_scalar_example.f90

program example
    use iso_fortran_env
    use curvefit_regression
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: n = 20
    real(real64), parameter :: xMin = 0.0d0
    real(real64), parameter :: xMax = 5.0d0
    real(real64), parameter :: slope = 10.0d0

    ! Local Variables
    real(real64) :: m, x(n), y(n), yfit(n)
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2
    type(legend), pointer :: lgnd
    character(len = 60) :: txt

    ! Initialization
    x = linspace(xMin, xMax, n)
    call random_number(y)
    y = slope * x + 2.0d0 * (y - 0.5d0)

    ! Fit the data
    m = linear_least_squares(x, y)
    yfit = m * x

    ! Plot the data
    call plt%initialize()
    call plt%set_font_size(14)

    lgnd => plt%get_legend()
    call lgnd%set_horizontal_position(LEGEND_LEFT)

    write(txt, '(AF7.3AF7.3)') "Actual Slope: ", slope, ", Fitted Slope: ", m
    call plt%set_title(trim(txt))

    call d1%set_name("Data")
    call d1%set_draw_line(.false.)
    call d1%set_draw_markers(.true.)
    call d1%set_marker_style(MARKER_X)
    call d1%set_line_color(CLR_BLACK)
    call d1%set_line_width(2.0)
    call d1%set_marker_scaling(1.5)
    call d1%define_data(x, y)

    call d2%set_name("Fit")
    call d2%set_line_color(CLR_BLUE)
    call d2%set_line_width(2.0)
    call d2%define_data(x, yfit)

    call plt%push(d1)
    call plt%push(d2)
    call plt%draw()
end program
