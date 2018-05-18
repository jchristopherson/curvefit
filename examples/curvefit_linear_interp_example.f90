! curvefit_linear_interp_example.f90

program example
    use iso_fortran_env
    use curvefit_interp
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: n = 9
    integer(int32), parameter :: m = 100

    ! Local Variables
    type(linear_interp) :: interp
    real(real64) :: dx, x(n), y(n), xi(m), yi(m)
    integer(int32) :: i
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2

    ! Initialization
    x = [-4.0d0, -3.0d0, -2.0d0, -1.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]
    y = [0.0d0, 0.15d0, 1.12d0, 2.36d0, 2.36d0, 1.46d0, 0.49d0, 0.06d0, &
        0.0d0]
    xi(1) = minval(x)
    dx = (maxval(x) - minval(x)) / (m - 1.0d0)
    do i = 2, m
        xi(i) = xi(i-1) + dx
    end do

    ! Interpolate
    call interp%initialize(x, y)
    yi = interp%interpolate(xi)

    ! Plot the results
    call plt%initialize()
    call plt%set_font_size(14)

    call d1%set_name("Raw Data")
    call d1%set_line_color(CLR_BLACK)
    call d1%set_draw_line(.false.)
    call d1%set_draw_markers(.true.)
    call d1%set_marker_style(MARKER_EMPTY_CIRCLE)
    call d1%set_marker_scaling(2.0)
    call d1%set_line_width(2.0)
    call d1%define_data(x, y)

    call d2%set_name("Interpolated")
    call d2%set_line_color(CLR_BLUE)
    call d2%set_line_width(2.0)
    call d2%define_data(xi, yi)

    call plt%push(d1)
    call plt%push(d2)
    call plt%draw()

end program
