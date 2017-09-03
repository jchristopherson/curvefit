# curvefit
A library for fitting functions to sets of data.

## TO DO
- Localized Regression or Nonparametric Fitting
    - http://www.aliquote.org/cours/2012_biomed/biblio/Cleveland1979.pdf
    - https://en.wikipedia.org/wiki/Local_regression
    - https://github.com/JuliaStats/Loess.jl/blob/master/src/Loess.jl
- Linear & Nonlinear Regression
- Statistics & Confidence Intervals - Related to quality of fit
    - https://www.nist.gov/itl/sed/datapac
- Outlier Removal
- Moving Averages
- Derivatives & Integrals of Fitted Curves
- Write tests
- C API

## Example 1
This example illustrates the use of cubic spline interpolation using both natural and forced boundary conditions.  Notice, the forced boundary conditions are arbitrarily chosen to illustrate their use.
```fortran
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
        if (i <= knotpts) then
            write(id, '(F14.10AF14.10AF14.10AF14.10AF14.10)') x(i), ",", &
                y(i), ",", xi(i), ",", y1(i), ",", y2(i)
        else
            write(id, '(AF14.10AF14.10AF14.10)') ",,", xi(i), ",", y1(i), &
                ",", y2(i)
        end if
    end do
    close(id)
end program
```
The above program yields the following data.
![](images/spline_interp_example_1.png?raw=true)