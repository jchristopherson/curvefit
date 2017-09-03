! curvefit_interp_example.f90

program example
    use curvefit_core
    use curvefit_interp
    implicit none

    ! Local Variables
    integer(i32), parameter :: knotpts = 9
    integer(i32), parameter :: npts = 1000
    integer(i32) :: i, id
    real(dp) :: dx, dstart, dend, x(knotpts), y(knotpts), xi(npts), y1(npts), &
        y2(npts), xmin, xmax
    type(spline_interp) :: interp

    ! Define a data set:
    x = [-4.0d0, -3.0d0, -2.0d0, -1.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]
    y = [0.0d0, 0.15d0, 1.12d0, 2.36d0, 2.36d0, 1.46d0, 0.49d0, 0.06d0, 0.0d0]

    ! Interpolate
    xmin = minval(x)
    xmax = maxval(x)
    dx = (xmax - xmin) / (npts - 1.0d0)
    xi(1) = xmin
    do i = 2, npts
        xi(i) = xi(i-1) + dx
    end do

    ! Allow for natural boundary conditions to the spline
    call interp%initialize(x, y)
    y1 = interp%interpolate(xi)

    ! Define the value of the first derivative at the end points
    dstart = 5.0d0
    dend = 0.0d0
    call interp%initialize_spline(x, y, &
        SPLINE_KNOWN_FIRST_DERIVATIVE, dstart, &
        SPLINE_KNOWN_FIRST_DERIVATIVE, dend)
    y2 = interp%interpolate(xi)

    ! Write the results to file
    open(newunit = id, file = "curvefit_interp.txt", action = "write", &
        status = "replace")
    do i = 1, max(knotpts, npts)
        if (knotpts < npts) then
            if (i <= knotpts) then
                write(id, '(F14.10AF14.10AF14.10AF14.10AF14.10)') x(i), ",", &
                    y(i), ",", xi(i), ",", y1(i), ",", y2(i)
            else
                write(id, '(AF14.10AF14.10AF14.10)') ",,", xi(i), ",", y1(i), &
                    ",", y2(i)
            end if
        else
            if (i <= npts) then
                write(id, '(F14.10AF14.10AF14.10AF14.10AF14.10)') x(i), ",", &
                    y(i), ",", xi(i), ",", y1(i), ",", y2(i)
            else
                write(id, '(F14.10AF14.10)') x(i), ",", y(i)
            end if
        end if
    end do
    close(id)
end program
