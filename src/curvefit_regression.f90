! curvefit_regression.f90

!> @brief \b curvefit_regression
!!
!! @par Purpose
!! To provide routines for perforing regression operations on sets of numerical
!! data.
module curvefit_regression
    use curvefit_core
    implicit none
    private

contains
! ******************************************************************************
! LOCAL REGRESSION - LOWESS
! ------------------------------------------------------------------------------
    !
    subroutine lowest(x, y, xs, ys, nleft, nright, w, userw, rw, ok)
        ! Arguments
        real(dp), intent(in), dimension(:) :: x, y, rw ! N ELEMENT
        real(dp), intent(in) :: xs
        real(dp), intent(out) :: ys
        integer(i32), intent(in) :: nleft, nright
        real(dp), intent(out), dimension(:) :: w ! N ELEMENT
        logical, intent(in) :: userw
        logical, intent(out) :: ok

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: p001 = 1.0d-3
        real(dp), parameter :: p999 = 0.999d0

        ! Local Variables
        integer(i32) :: j, n, nrt
        real(dp) :: range, h, h9, h1, a, b, c, r

        ! Initialization
        n = size(x)
        range = x(n) - x(1)
        h = max(xs - x(nleft), x(nright) - xs)
        h9 = p999 * h
        h1 = p001 * h
        a = zero

        ! Process
        do j = nleft, n
            w(j) = zero
            r = abs(x(j) - xs)
            if (r <= h9) then
                if (r > h1) then
                    w(j) = (one - (r / h)**3)**3
                else
                    w(j) = one
                end if
                if (userw) w(j) = rw(j) * w(j)
                a = a + w(j)
            else if (x(j) > xs) then
                exit
            end if
        end do

        nrt = j - 1
        if (a <= zero) then
            ok = .false.
        else
            ok = .true.
            w(nleft:nrt) = w(nleft:nrt) / a
            if (h > zero) then
                a = zero
                do j = nleft, nrt
                    a = a + w(j) * x(j)
                end do
                b = xs - a
                c = zero
                do j = nleft, nrt
                    c = c + w(j) * (x(j) - a)**2
                end do
                if (sqrt(c) > p001 * range) then
                    b = b / c
                    do j = nleft, nrt
                        w(j) = w(j) * (one + b * (x(j) - a))
                    end do
                end if
            end if
            ys = zero
            do j = nleft, nrt
                ys = ys + w(j) * y(j)
            end do
        end if
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
