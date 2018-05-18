! curvefit_test_interp.f90

! Tests the interpolation routines.
module curvefit_test_interp
    use iso_fortran_env
    use curvefit_core
    use curvefit_interp
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_linear_interp() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 9
        integer(int32), parameter :: m = 1000

        ! Local Variables
        logical :: rst
        type(linear_interp) :: interp
        real(real64) :: dx, x(n), y(n), xi(m), yi(m)
        integer(int32) :: i, id

        ! Initialization
        rst = .true.
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

        ! Write the results to a text file so we can plot them
        open(newunit = id, file = "linear_interp.txt", action = "write", &
            status = "replace")
        do i = 1, max(m, n)
            if (m < n) then
                if (i <= m) then
                    write(id, '(F14.10AF14.10AF14.10AF14.10)') x(i), ",", &
                        y(i), ",", xi(i), ",", yi(i)
                else
                    write(id, '(F14.10AF14.10)') x(i), ",", y(i)
                end if
            else
                if (i <= n) then
                    write(id, '(F14.10AF14.10AF14.10AF14.10)') x(i), ",", &
                        y(i), ",", xi(i), ",", yi(i)
                else
                    write(id, '(AF14.10AF14.10)') ",,", xi(i), ",", yi(i)
                end if
            end if
        end do
        close(id)
    end function

! ------------------------------------------------------------------------------
    function test_poly_interp() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 9
        integer(int32), parameter :: m = 1000

        ! Local Variables
        logical :: rst
        type(polynomial_interp) :: interp
        real(real64) :: dx, x(n), y(n), xi(m), yi(m)
        integer(int32) :: i, id

        ! Initialization
        rst = .true.
        x = [-4.0d0, -3.0d0, -2.0d0, -1.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]
        y = [0.0d0, 0.15d0, 1.12d0, 2.36d0, 2.36d0, 1.46d0, 0.49d0, 0.06d0, &
            0.0d0]
        xi(1) = minval(x)
        dx = (maxval(x) - minval(x)) / (m - 1.0d0)
        do i = 2, m
            xi(i) = xi(i-1) + dx
        end do

        ! Interpolate
        call interp%initialize(x, y, 2)
        yi = interp%interpolate(xi)

        ! Write the results to a text file so we can plot them
        open(newunit = id, file = "poly_interp.txt", action = "write", &
            status = "replace")
        do i = 1, max(m, n)
            if (m < n) then
                if (i <= m) then
                    write(id, '(F14.10AF14.10AF14.10AF14.10)') x(i), ",", &
                        y(i), ",", xi(i), ",", yi(i)
                else
                    write(id, '(F14.10AF14.10)') x(i), ",", y(i)
                end if
            else
                if (i <= n) then
                    write(id, '(F14.10AF14.10AF14.10AF14.10)') x(i), ",", &
                        y(i), ",", xi(i), ",", yi(i)
                else
                    write(id, '(AF14.10AF14.10)') ",,", xi(i), ",", yi(i)
                end if
            end if
        end do
        close(id)
    end function

! ------------------------------------------------------------------------------
    function test_spline_interp() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 9
        integer(int32), parameter :: m = 100

        ! Local Variables
        logical :: rst
        type(spline_interp) :: interp
        real(real64) :: dx, x(n), y(n), xi(m), yi(m)
        integer(int32) :: i, id

        ! Initialization
        rst = .true.
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

        ! Write the results to a text file so we can plot them
        open(newunit = id, file = "spline_interp.txt", action = "write", &
            status = "replace")
        do i = 1, max(m, n)
            if (m < n) then
                if (i <= m) then
                    write(id, '(F14.10AF14.10AF14.10AF14.10)') x(i), ",", &
                        y(i), ",", xi(i), ",", yi(i)
                else
                    write(id, '(F14.10AF14.10)') x(i), ",", y(i)
                end if
            else
                if (i <= n) then
                    write(id, '(F14.10AF14.10AF14.10AF14.10)') x(i), ",", &
                        y(i), ",", xi(i), ",", yi(i)
                else
                    write(id, '(AF14.10AF14.10)') ",,", xi(i), ",", yi(i)
                end if
            end if
        end do
        close(id)
    end function

! ------------------------------------------------------------------------------

end module
