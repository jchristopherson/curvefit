! curvefit_moving_avg_example.f90

program example
    use iso_fortran_env
    use curvefit_regression
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: n = 1000
    integer(int32), parameter :: nAvg = 20
    real(real64), parameter :: xMin = 0.0d0
    real(real64), parameter :: xMax = 1.0d0
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2

    ! Local Variables
    real(real64) :: x(n), y(n), yAvg(n)

    ! Initialization
    x = linspace(xMin, xMax, n)
    call random_number(y)
    y = (y - 0.5d0) + &
        sin(15.0d0 * x) + &
        0.5d0 * sin(25.0d0 * x) + &
        0.1d0 * sin(75.0d0 * x)
    yAvg = y

    ! Apply the moving average to the data set
    call moving_average(yAvg, nAvg)

    ! Plot the results
    call plt%initialize()
    call plt%set_font_size(14)

    call d1%set_name("Noisy Data")
    call d1%set_line_color(CLR_BLUE)
    call d1%define_data(x, y)

    call d2%set_name("Averaged")
    call d2%set_line_color(CLR_RED)
    call d2%set_line_width(3.0)
    call d2%define_data(x, yAvg)

    call plt%push(d1)
    call plt%push(d2)
    call plt%draw()
end program
