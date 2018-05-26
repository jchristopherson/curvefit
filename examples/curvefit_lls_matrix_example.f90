! curvefit_lls_matrix_example.f90

! The following program illustrates how to generate a linear fit to a
! multi-input system.  The following data was taken from an axial-torsional
! load cell calibration.
!
! Axial Load Data:       Torsional Load Data:
! 0	                     0
! 3000	                 0
! 6000	                 0
! 7500	                 0
! 9000	                 0
! 12000	                 0
! 15000	                 0
! 7500	                 0
! 0	                     0
! 0	                     0
! -3000	                 0
! -6000	                 0
! -7500	                 0
! -9000	                 0
! -12000	             0
! -15000	             0
! -7500	                 0
! 0	                     0
! 0	                     0
! 0	                     67.7908728
! 0	                     135.5817456
! 0	                     203.3726184
! 0	                     271.1634912
! 0	                     338.954364
! 0	                     203.3726184
! 0	                     0
! 0	                     0
! 0	                     -67.7908728
! 0	                     -135.5817456
! 0	                     -203.3726184
! 0	                     -271.1634912
! 0	                     -338.954364
! 0	                     -203.3726184
! 0	                     0
!
! Bridge 1 Output:       Bridge 2 Output:
! 0	                     0
! 0.38905	             0.00122
! 0.77816	             0.00259
! 0.97269             	 0.0029
! 1.16714	             0.00314
! 1.556	                 0.00338
! 1.94484	             0.00356
! 0.9726	             0.00477
! -0.00001	             -0.00001
! 0	                     0
! -0.388886	             0.00021
! -0.77775	             0.00051
! -0.97215	             0.00069
! -1.16654	             0.00088
! -1.55533	             0.0013
! -1.9441	             0.00175
! -0.97171	             0.00058
! 0.00004	             0.00003
! 0	                     0
! -0.00044	             0.27156
! -0.0013	             0.54329
! -0.0024	             0.81507
! -0.00382	             1.08682
! -0.00528	             1.35881
! -0.00257	             0.81553
! 0.00015	             0.0001
! 0	                     0
! 0.00144	             -0.27145
! 0.00306	             -0.54312
! 0.00446	             -0.81493
! 0.00567	             -1.0868
! 0.00688	             -1.35879
! 0.00451	             -0.81548
! -0.00002	             0
program example
    use iso_fortran_env
    use curvefit_regression
    use fplot_core
    implicit none

    ! Local Variables
    integer(int32), parameter :: npts = 34
    integer(int32), parameter :: nchan = 2
    real(real64) :: loads(npts, nchan), bridge(npts, nchan), &
        gain(nchan, nchan), calib(npts, nchan), errs(npts, nchan), &
        axialFullScale, torqueFullScale
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2
    class(plot_axis), pointer :: xAxis, yAxis

    ! Populate the applied loads matrix
    loads = reshape([0.0d0, 3000.0d0, 6000.0d0, 7500.0d0, 9000.0d0, 12000.0d0, &
        15000.0d0, 7500.0d0, 0.0d0, 0.0d0, -3000.0d0, -6000.0d0, -7500.0d0, &
        -9000.0d0, -12000.0d0, -15000.0d0, -7500.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 67.7908728d0, 135.5817456d0, 203.3726184d0, &
        271.1634912d0, 338.954364d0, 203.3726184d0, 0.0d0, 0.0d0, &
        -67.7908728d0, -135.5817456d0, -203.3726184d0, -271.1634912d0, &
        -338.954364d0, -203.3726184d0, 0.0d0], [npts, nchan])

    ! Populate the bridge output (measured) matrix
    bridge = reshape([0.0d0, 0.38905d0, 0.77816d0, 0.97269d0, 1.16714d0, &
        1.556d0, 1.94484d0, 0.9726d0, -0.00001d0, 0.0d0, -0.388886d0, &
        -0.77775d0, -0.97215d0, -1.16654d0, -1.55533d0, -1.9441d0, -0.97171d0, &
        0.00004d0, 0.0d0, -0.00044d0, -0.0013d0, -0.0024d0, -0.00382d0, &
        -0.00528d0, -0.00257d0, 0.00015d0, 0.0d0, 0.00144d0, 0.00306d0, &
        0.00446d0, 0.00567d0, 0.00688d0, 0.00451d0, -0.00002d0, 0.0d0, &
        0.00122d0, 0.00259d0, 0.0029d0, 0.00314d0, 0.00338d0, 0.00356d0, &
        0.00477d0, -0.00001d0, 0.0d0, 0.00021d0, 0.00051d0, 0.00069d0, &
        0.00088d0, 0.0013d0, 0.00175d0, 0.00058d0, 0.00003d0, 0.0d0, &
        0.27156d0, 0.54329d0, 0.81507d0, 1.08682d0, 1.35881d0, 0.81553d0, &
        0.0001d0, 0.0d0, -0.27145d0, -0.54312d0, -0.81493d0, -1.0868d0, &
        -1.35879d0, -0.81548d0, 0.0d0], [npts, nchan])

    ! Compute the coefficient matrix.  The transpose operations are necessary as
    ! the data was input in column-major format in the order illustrated in the
    ! comments.  Had the data been input differently, the transpose operations
    ! could have been avoided.
    gain = linear_least_squares(transpose(bridge), transpose(loads))

    ! Apply the matrix
    calib = matmul(bridge, transpose(gain))

    ! Compute errors in the fit
    errs = calib - loads

    ! Express the loads in terms of percent of full scale
    axialFullScale = maxval(loads(:,1))
    torqueFullScale = maxval(loads(:,2))
    errs(:,1) = 1.0d2 * errs(:,1) / axialFullScale
    errs(:,2) = 1.0d2 * errs(:,2) / torqueFullScale

    ! Plot the error data vs. index of the applied load
    call plt%initialize()
    call plt%set_font_size(14)

    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    call xAxis%set_title("Index")
    call yAxis%set_title("Error [% FS]")

    call d1%set_name("Axial")
    call d1%set_line_color(CLR_BLUE)
    call d1%set_line_width(2.0)
    call d1%set_draw_markers(.true.)
    call d1%set_marker_style(MARKER_X)
    call d1%define_data(errs(:,1))

    call d2%set_name("Torsional")
    call d2%set_line_color(CLR_RED)
    call d2%set_line_width(2.0)
    call d2%set_draw_markers(.true.)
    call d2%set_marker_style(MARKER_EMPTY_CIRCLE)
    call d2%define_data(errs(:,2))

    call plt%push(d1)
    call plt%push(d2)
    call plt%draw()
end program
