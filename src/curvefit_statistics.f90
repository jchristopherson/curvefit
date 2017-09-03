! curvefit_statistics.f90

!> @brief \b curvefit_statistics
!!
!! @par Purpose
!! To provide a set of statistical routines for exploring curve fits of sets 
!! of numeric data.
module curvefit_statistics
    use curvefit_core
    implicit none
    private
    public :: mean
    public :: variance
    public :: standard_deviation

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Computes the mean of a data set.
    interface mean
        module procedure :: mean_dbl
    end interface

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    !> @brief Computes the sample variance of a data set.
    interface variance
        module procedure :: variance_dbl
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the corrected standard deviation of a data set.
    interface standard_deviation
        module procedure :: stdev_dbl
    end interface

contains
! ******************************************************************************
! MEAN, MEDIAN, ETC.
! ------------------------------------------------------------------------------
    !> @brief Computes the mean of a data set.
    !!
    !! @param[in] x The data set.
    !!
    !! @return The mean of x.
    !!
    !! @par Remarks
    !! To avoid overflow-type issues, Knuth's algorithm is employed.  A simple
    !! illustration of this algorithm can be found 
    !! [here](https://www.johndcook.com/blog/standard_deviation/).
    pure function mean_dbl(x) result(z)
        ! Arguments
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: z

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: i, n

        ! Process
        n = size(x)
        if (n == 0) then
            z = zero
        else
            z = x(1)
            do i = 2, n
                z = z + (x(i) - z) / i
            end do
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the median of a data set.
    !!
    !! @param[in] x The monotonic data set whose median is to be found.
    !!
    !! @return The median of @p x.
    ! function median_double(x) result(z)
    !     ! Arguments
    !     real(dp), intent(inout), dimension(:) :: x
    !     real(dp) :: z

    !     ! Parameters
    !     real(dp), parameter :: zero = 0.0d0
    !     real(dp), parameter :: half = 0.5d0

    !     ! Local Variables
    !     integer(i32) :: n, iflag, nmid, nmidp1

    !     ! Initialization
    !     n = size(x)
    !     nmid = n / 2
    !     nmidp1 = nmid + 1
    !     iflag = n - 2 * nmid

    !     ! Quick Return
    !     if (n == 0) then
    !         z = zero
    !         return
    !     end if

    !     ! Sort the array into ascending order
    !     call sort_array(.true., x)

    !     ! Compute the median
    !     if (iflag == 0) then
    !         z = half * (x(nmid) + x(nmidp1))
    !     else
    !         z = x(nmidp1)
    !     end if
    ! end function

! ------------------------------------------------------------------------------
    !> @brief Computes the sample variance of a data set.
    !!
    !! @param[in] x The data set.
    !!
    !! @return The variance of @p x.
    !!
    !! @par Remarks
    !! To avoid overflow-type issues, Knuth's algorithm is employed.  A simple
    !! illustration of this algorithm can be found 
    !! [here](https://www.johndcook.com/blog/standard_deviation/).
    pure function variance_dbl(x) result(v)
        ! Arguments
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: v

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, n
        real(dp) :: oldMean, newMean

        ! Process
        n = size(x)
        if (n <= 1) then
            v = zero
        else
            oldMean = x(1)
            v = zero
            do i = 2, n
                newMean = oldMean + (x(i) - oldMean) / i
                v = v + (x(i) - oldMean) * (x(i) - newMean)
                oldMean = newMean
            end do
            v = v / (n - 1)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the corrected standard deviation of a data set.
    !!
    !! @param[in] x The data set.
    !!
    !! @return The standard deviation of @p x.
    pure function stdev_dbl(x) result(s)
        ! Arguments
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: s

        ! Process
        s = sqrt(variance(x))
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
