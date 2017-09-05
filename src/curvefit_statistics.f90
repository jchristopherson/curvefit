! curvefit_statistics.f90

!> @brief \b curvefit_statistics
!!
!! @par Purpose
!! To provide a set of statistical routines for exploring curve fits of sets 
!! of numeric data.
module curvefit_statistics
    use curvefit_core
    use ferror, only : errors
    implicit none
    private
    public :: mean
    public :: variance
    public :: standard_deviation
    public :: confidence_interval
    public :: z_value
    public :: incomplete_gamma
    public :: incomplete_gamma_comp


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

! ------------------------------------------------------------------------------
    !> @brief Computes the confidence interval based upon a standard normal 
    !! distribution.
    interface confidence_interval
        module procedure :: conf_int
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the z-value (z-score) given a percentage of the area 
    !! under the standard normal distribution curve.
    interface z_value
        module procedure :: std_norm_dist_z_score
    end interface

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

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
    !> @brief Computes the confidence interval based upon a standard normal 
    !! distribution.
    !!
    !! @param[in] x The data set.
    !! @param[in] alpha The confidence level.  This value must lie between
    !! zero and one such that: 0 < alpha < 1.
    !!
    !! @return The confidence interval as the deviation from the mean.
    !!
    !! @par Remarks
    !! The confidence interval, assuming a standard normal distribution, is
    !! as follows: mu +/- z * s / sqrt(n), where mu = the mean, and s = the
    !! standard deviation.  This routine computes the z * s / sqrt(n) portion 
    !! leaving the computation of the mean to the user.
    function conf_int(x, alpha) result(ci)
        ! Arguments
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(in) :: alpha
        real(dp) :: ci

        ! Local Variables
        real(dp) :: n, sigma, z

        ! Ensure: 0 < alpha < 1

        ! Compute the standard deviation, and z-distribution value
        sigma = standard_deviation(x)
        n = real(size(x), dp)
        z = z_value(alpha)

        ! Compute the confidence interval - offset from the mean
        ci = z * sigma / sqrt(n)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the z-value (z-score) given a percentage of the area 
    !! under the standard normal distribution curve.
    !!
    !! @param[in] alpha The percentage of the area under the curve.  This value
    !!  must be between 0 and 1 such that: 0 < alpha < 1.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p alpha is does not satisfy:
    !!      0 < alpha < 1.
    !!
    !! @return The z-score or z-value by solving for z where: 
    !!  alpha = ERF(z / sqrt(2)), where ERF is the error function.
    function std_norm_dist_z_score(alpha, err) result(z)
        ! Supporting Modules
        use nonlin_types, only : fcn1var, fcn1var_helper, value_pair
        use nonlin_solve, only : brent_solver

        ! Arguments
        real(dp), intent(in) :: alpha
        real(dp) :: z
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: two = 2.0d0
        real(dp), parameter :: ten = 1.0d1

        ! Local Variables
        type(fcn1var_helper) :: obj
        procedure(fcn1var), pointer :: fcn
        type(brent_solver) :: solver
        type(value_pair) :: lim
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (alpha <= zero .or. alpha >= one) then
            call errmgr%report_error("std_norm_dist_z_score", &
                "The alpha parameter must lie between 0 and 1.", &
                CF_INVALID_INPUT_ERROR)
        end if

        ! Compute the solution
        fcn => zfun
        call obj%set_fcn(fcn)
        lim%x1 = zero
        lim%x2 = ten
        call solver%solve(obj, z, lim, err = errmgr)

    contains
        ! Compute the solution to: alpha = erf(z / sqrt(2)) for z
        function zfun(x) result(f)
            real(dp), intent(in) :: x
            real(dp) :: f
            f = alpha - erf(x / sqrt(two))
        end function
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ******************************************************************************
! GAMMA FUNCTION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the incomplete gamma function: 
    !! P(a,x) = 1 / gamma(a) * integrate(exp(-t) * t**(a - 1), t, 0, x).
    !!
    !! @param[in] a The coefficient.  This parameter must be positive-valued.
    !! @param[in] x The independent variable.  This parameter must be greater
    !!  than or equal to zero.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p x is negative, or if @p a is not
    !!      positive.
    !!
    !! @param return The value of the function at @p x.
    !!
    !! @par Remarks
    !! This implementation is based upon the Numerical Recipes implementation
    !! found in section 6.2 of the text (routine: gammp).
    function incomplete_gamma(a, x, err) result(g)
        ! Arguments
        real(dp), intent(in) :: a, x
        class(errors), intent(inout), optional, target :: err
        real(dp) :: g

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        g = zero

        ! Input Check
        if (x < zero) then
            call errmgr%report_error("incomplete_gamma", &
                "The independent variable (x) must be >= 0.", &
                CF_INVALID_INPUT_ERROR)
            return
        else if (a <= zero) then
            call errmgr%report_error("incomplete_gamma", &
                "The parameter a must be positive.", CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Process
        if (x < a + one) then
            g = inc_gamma_series(a, x)
        else
            g = one - inc_gamma_cf(a, x)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the complement of the incomplete gamma function: 
    !! Q(a,x) = 1 - P(a,x), where 
    !! P(a,x) = 1 / gamma(a) * integrate(exp(-t) * t**(a - 1), t, 0, x).
    !!
    !! @param[in] a The coefficient.  This parameter must be positive-valued.
    !! @param[in] x The independent variable.  This parameter must be greater
    !!  than or equal to zero.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p x is negative, or if @p a is not
    !!      positive.
    !!
    !! @param return The value of the function at @p x.
    !!
    !! @par Remarks
    !! This implementation is based upon the Numerical Recipes implementation
    !! found in section 6.2 of the text (routine: gammq).
    function incomplete_gamma_comp(a, x, err) result(g)
        ! Arguments
        real(dp), intent(in) :: a, x
        class(errors), intent(inout), optional, target :: err
        real(dp) :: g

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        g = zero

        ! Input Check
        if (x < zero) then
            call errmgr%report_error("incomplete_gamma_comp", &
                "The independent variable (x) must be >= 0.", &
                CF_INVALID_INPUT_ERROR)
            return
        else if (a <= zero) then
            call errmgr%report_error("incomplete_gamma_comp", &
                "The parameter a must be positive.", CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Process
        if (x < a + one) then
            g = one - inc_gamma_series(a, x)
        else
            g = inc_gamma_cf(a, x)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the incomplete gamma function: 
    !! P(a,x) = 1 / gamma(a) * integrate(exp(-t) * t**(a - 1), t, 0, x)
    !! is computed by its series representation.
    !!
    !! @param[in] a The parameter a.
    !! @param[in] x The parameter x.  This parameter must be greater than 0.
    !!
    !! @return The incomplete gamma function.
    !!
    !! @par Remarks
    !! This implementation is based upon the Numerical Recipes implementation
    !! found in section 6.2 of the text (routine: gser).
    pure function inc_gamma_series(a, x) result(g)
        ! Arguments
        real(dp), intent(in) :: a, x
        real(dp) :: g

        ! Parameters
        integer(i32), parameter :: itmax = 500
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: n
        real(dp) :: ap, del, sm, gln, eps

        ! Initialization
        eps = epsilon(eps)
        gln = log_gamma(a)
        ap = a
        sm = one / a
        del = sm
        g = zero
        do n = 1, itmax
            ap = ap + one
            del = del * x / ap
            sm = sm + del
            if (abs(del) < abs(sm) * eps) then
                g = sm * exp(-x + a * log(x) - gln)
                return
            end if
        end do

        ! If we're here, it didn't converge
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the incomplete gamma function: 
    !! Q(a,x) = 1 - P(a,x), where 
    !! P(a,x) = 1 / gamma(a) * integrate(exp(-t) * t**(a - 1), t, 0, x)
    !! is computed by Lentz's continued fraction approach.
    !!
    !! @param[in] a The parameter a.
    !! @param[in] x The parameter x.  This parameter must be greater than 0.
    !!
    !! @return The incomplete gamma function.
    !!
    !! @par Remarks
    !! This implementation is based upon the Numerical Recipes implementation
    !! found in section 6.2 of the text (routine: gcf).
    pure function inc_gamma_cf(a, x) result(g)
        ! Arguments
        real(dp), intent(in) :: a, x
        real(dp) :: g

        ! Parameters
        integer(i32), parameter :: itmax = 500
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: two = 2.0d0

        ! Local Variables
        integer(i32) :: i
        real(dp) :: an, b, c, d, del, h, gln, eps, fpmin

        ! Initialization
        eps = epsilon(eps)
        fpmin = tiny(fpmin) / eps
        gln = log_gamma(a)
        b = x + one - a
        c = one / fpmin
        d = one / b
        h = d
        g = zero
        do i = 1, itmax
            an = -i * (i - a)
            b = b + two
            d = an * d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an / c
            if (abs(c) < fpmin) c = fpmin
            d = one / d
            del = d * c
            h = h * del
            if (abs(del - one) < eps) then
                g = exp(-x + a * log(x) - gln) * h
                return
            end if
        end do

        ! If we're here, it didn't converge
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
