! curvefit_test_interp.f90

! Tests the interpolation routines.
module curvefit_test_interp
    use curvefit_core
    use curvefit_interp
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_linear_interp() result(rst)
        ! Parameters
        integer(i32), parameter :: n = 100
        integer(i32), parameter :: m = 1000

        ! Local Variables
        logical :: rst
        type(linear_interp) :: interp
        real(dp) :: dx, x(n), y(n), xi(m), yi(m)
        integer(i32) :: i, id

        ! Initialization
        rst = .true.
        dx = 0.1d0
        x(1) = 0.0d0
        y(1) = 0.0d0
        do i = 2, n
            x(i) = x(i-1) + dx
            y(i) = exp(-0.1d0 * x(i)) * sin(5.0d1 * x(i))
        end do

        xi(1) = -0.1d0
        dx = dx * (real(n, dp) / m)
        do i = 2, m
            xi(i) = xi(i-1) + dx
        end do

        ! Interpolate
        call interp%initialize(x, y, 1)
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
end module
