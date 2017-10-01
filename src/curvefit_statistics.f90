! curvefit_statistics.f90

!> @brief \b curvefit_statistics
!!
!! @par Purpose
!! To provide a set of statistical routines for exploring curve fits of sets 
!! of numeric data.
module curvefit_statistics
    use curvefit_core
    use ferror, only : errors
    use linalg_sorting, only : sort
    implicit none
    private
    public :: mean
    public :: median
    public :: variance
    public :: covariance
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
    !> @brief Computes the median of a data set.
    interface median
        module procedure :: median_dbl
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the sample variance of a data set.
    interface variance
        module procedure :: variance_dbl
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the covariance matrix of two data sets.
    interface covariance
        module procedure :: covariance_2sets
        module procedure :: covariance_mtx
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
    !> @brief Computes the incomplete gamma function: 
    !! P(a,x) = 1 / gamma(a) * integrate(exp(-t) * t**(a - 1), t, 0, x).
    interface incomplete_gamma
        module procedure :: incomplete_gamma_scalar
        module procedure :: incomplete_gamma_array
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the complement of the incomplete gamma function: 
    !! Q(a,x) = 1 - P(a,x), where 
    !! P(a,x) = 1 / gamma(a) * integrate(exp(-t) * t**(a - 1), t, 0, x).
    interface incomplete_gamma_comp
        module procedure :: incomplete_gamma_comp_scalar
        module procedure :: incomplete_gamma_comp_array
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the incomplete beta function:
    !! I(a,b) = 1 / B(a,b) * integrate(t**(a - 1) * (1 - t)**(b - 1), t, 0, x)
    interface incomplete_beta
        module procedure :: inc_beta_scalar
        module procedure :: inc_beta_array
    end interface


contains
! ******************************************************************************
! MEAN, MEDIAN, ETC.
! ------------------------------------------------------------------------------
    !> @brief Computes the mean of a data set.
    !!
    !! @param[in] x The data set.
    !!
    !! @return The mean of @p x.
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
    !! @param[in,out] x The data set whose median is to be found.  Ideally, the
    !!  data set should be monotonically increasing; however, if it is not, it
    !!  may be sorted by the routine, dependent upon the value of @p srt.  On
    !!  output, the array contents are unchanged; however, they may be sorted
    !!  into ascending order (dependent upon the value of @p srt).
    !! @param[in] srt An optional flag determining if @p x should be sorted. 
    !!  The default is to sort (true).
    !!
    !! @return The median of @p x.
    function median_dbl(x, srt) result(z)
        ! Arguments
        real(dp), intent(inout), dimension(:) :: x
        logical, intent(in), optional :: srt
        real(dp) :: z

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: half = 0.5d0

        ! Local Variables
        logical :: sortData
        integer(i32) :: n, iflag, nmid, nmidp1

        ! Initialization
        n = size(x)
        nmid = n / 2
        nmidp1 = nmid + 1
        iflag = n - 2 * nmid
        sortData = .true.
        if (present(srt)) sortData = srt

        ! Quick Return
        if (n == 0) then
            z = zero
            return
        end if

        ! Sort the array into ascending order
        if (sortData) call sort(x, .true.)

        ! Compute the median
        if (iflag == 0) then
            z = half * (x(nmid) + x(nmidp1))
        else
            z = x(nmidp1)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the sample variance of a data set.
    !!
    !! @param[in] x The data set.
    !!
    !! @return The variance of @p x.
    !!
    !! @par Remarks
    !! To avoid overflow-type issues, Welford's algorithm is employed.  A simple
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
            v = v / (n - 1.0d0)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the covariance matrix of two data sets.
    !!
    !! @param[in] x An N-element array containing the first data set.
    !! @param[in] y An N-element array containing the second data set.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
    !!
    !! @return The 2-by-2 covariance matrix.
    function covariance_2sets(x, y, err) result(c)
        ! Arguments
        real(dp), intent(in), dimension(:) :: x, y
        class(errors), intent(inout), optional, target :: err
        real(dp), dimension(2,2) :: c

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, n
        real(dp) :: oldMeanX, newMeanX, oldMeanY, newMeanY
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(x)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(y) /= n) then
            ! ERROR
            call errmgr%report_error("covariance_2sets", &
                "The input arrays must be the same size.", CF_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        c = zero
        if (n > 1) then
            oldMeanX = x(1)
            oldMeanY = y(1)
            do i = 2, n
                newMeanX = oldMeanX + (x(i) - oldMeanX) / i
                newMeanY = oldMeanY + (y(i) - oldMeanY) / i
                c(1,1) = c(1,1) + (x(i) - oldMeanX) * (x(i) - newMeanX)
                c(2,2) = c(2,2) + (y(i) - oldMeanY) * (y(i) - newMeanY)
                c(1,2) = c(1,2) + (x(i) - oldMeanX) * (y(i) - newMeanY)
                oldMeanX = newMeanX
                oldMeanY = newMeanY
            end do
            c(2,1) = c(1,2)
            c = c / (n - one)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the covariance matrix of N data sets of M observations.
    !!
    !! @param[in] x The M-by-N matrix.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!
    !! @return The N-by-N covariance matrix.
    function covariance_mtx(x, err) result(c)
        ! Arguments
        real(dp), intent(in), dimension(:,:) :: x
        class(errors), intent(inout), optional, target :: err
        real(dp), dimension(size(x, 2), size(x, 2)) :: c

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, k, m, n, flag
        real(dp), allocatable, dimension(:) :: oldMeans, newMeans
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(x, 1)
        n = size(x, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Local Memory Allocation
        allocate(oldMeans(n), stat = flag)
        if (flag == 0) allocate(newMeans(n), stat = flag)
        if (flag /= 0) then
            ! ERROR
            call errmgr%report_error("covariance_mtx", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Process
        c = zero
        oldMeans = x(1,:)
        do i = 2, m
            newMeans = oldMeans + (x(i,:) - oldMeans) / i
            do k = 1, n
                c(1:k,k) = c(1:k,k) + &
                    (x(i,1:k) - oldMeans(1:k)) * (x(i,k) - newMeans(k))
            end do
            oldMeans = newMeans
        end do
        do k = 2, n
            c(k,1:k-1) = c(1:k-1,k)
        end do
        c = c / (m - one)
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

! ******************************************************************************
! STUDENT'S T TEST ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes Student's T test decision for the null hypothesis that 
    !! the data vectors come from independent samples from normal distributions
    !! with equal means, and equal but unknown variances.
    !!
    !! @param[in] data1 The first data set.
    !! @param[in] data2 The second data set.
    !! @param[out] p
    !! 
    !! @return 
    function ttest_2sets(data1, data2, p) result(t)
        ! Arguments
        real(dp), intent(in), dimension(:) :: data1, data2
        real(dp), intent(out), optional :: p

        ! Parameter
        real(dp), parameter :: half = 0.5d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: two = 2.0d0

        ! Local Variables
        integer(i32) :: n1, n2
        real(dp) :: m1, m2, v1, v2, df, svar

        ! Initialization
        n1 = size(data1)
        n2 = size(data2)

        ! Compute the mean and variance of both data sets
        m1 = mean(data1)
        v1 = variance(data1)
        m2 = mean(data2)
        v2 = variance(data2)

        ! Process
        df = n1 + n2 - two
        svar = ((n - 1) * v1 + (n2 - 1) * v2) / df
        t = (m1 - m2) / sqrt(svar * (one / n1 + one / n2))
        if (present(p)) then
            p = incomplete_beta(half * df, half, df / (df + t**2))
        end if
    end function

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
    function incomplete_gamma_scalar(a, x, err) result(g)
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
            call errmgr%report_error("incomplete_gamma_scalar", &
                "The independent variable (x) must be >= 0.", &
                CF_INVALID_INPUT_ERROR)
            return
        else if (a <= zero) then
            call errmgr%report_error("incomplete_gamma_scalar", &
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
    !> @brief Computes the incomplete gamma function: 
    !! P(a,x) = 1 / gamma(a) * integrate(exp(-t) * t**(a - 1), t, 0, x).
    !!
    !! @param[in] a The coefficient.  This parameter must be positive-valued.
    !! @param[in] x An N-element array of independent variables.  All values
    !!  must be greater than or equal to zero.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p x is negative, or if @p a is not
    !!      positive.
    !!
    !! @param return The values of the function at @p x.
    !!
    !! @par Remarks
    !! This implementation is based upon the Numerical Recipes implementation
    !! found in section 6.2 of the text (routine: gammp).
    function incomplete_gamma_array(a, x, err) result(g)
        ! Arguments
        real(dp), intent(in) :: a
        real(dp), intent(in), dimension(:) :: x
        class(errors), intent(inout), optional, target :: err
        real(dp), dimension(size(x)) :: g

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        logical :: check
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(i32) :: i, n

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        g = zero
        n = size(x)

        ! Input Check
        check = .true.
        do i = 1, n
            if (x(i) < zero) then
                check = .false.
                exit
            end if
        end do
        if (.not.check) then
            call errmgr%report_error("incomplete_gamma_array", &
                "The independent variable (x) must be >= 0.", &
                CF_INVALID_INPUT_ERROR)
            return
        else if (a <= zero) then
            call errmgr%report_error("incomplete_gamma_array", &
                "The parameter a must be positive.", CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Process
        do concurrent (i = 1:n)
            if (x(i) < a + one) then
                g(i) = inc_gamma_series(a, x(i))
            else
                g(i) = one - inc_gamma_cf(a, x(i))
            end if
        end do
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
    function incomplete_gamma_comp_scalar(a, x, err) result(g)
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
            call errmgr%report_error("incomplete_gamma_comp_scalar", &
                "The independent variable (x) must be >= 0.", &
                CF_INVALID_INPUT_ERROR)
            return
        else if (a <= zero) then
            call errmgr%report_error("incomplete_gamma_comp_scalar", &
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
    !> @brief Computes the complement of the incomplete gamma function: 
    !! Q(a,x) = 1 - P(a,x), where 
    !! P(a,x) = 1 / gamma(a) * integrate(exp(-t) * t**(a - 1), t, 0, x).
    !!
    !! @param[in] a The coefficient.  This parameter must be positive-valued.
    !! @param[in] x An N-element array of independent variables.  All values
    !!  must be greater than or equal to zero.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p x is negative, or if @p a is not
    !!      positive.
    !!
    !! @param return The values of the function at @p x.
    !!
    !! @par Remarks
    !! This implementation is based upon the Numerical Recipes implementation
    !! found in section 6.2 of the text (routine: gammq).
    function incomplete_gamma_comp_array(a, x, err) result(g)
        ! Arguments
        real(dp), intent(in) :: a
        real(dp), intent(in), dimension(:) :: x
        class(errors), intent(inout), optional, target :: err
        real(dp), dimension(size(x)) :: g

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(i32) :: i, n
        logical :: check

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        g = zero
        n = size(x)

        ! Input Check
        check = .true.
        do i = 1, n
            if (x(i) < zero) then
                check = .false.
                exit
            end if
        end do
        if (.not.check) then
            call errmgr%report_error("incomplete_gamma_comp_array", &
                "The independent variable (x) must be >= 0.", &
                CF_INVALID_INPUT_ERROR)
            return
        else if (a <= zero) then
            call errmgr%report_error("incomplete_gamma_comp_array", &
                "The parameter a must be positive.", CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Process
        do concurrent (i = 1:n)
            if (x(i) < a + one) then
                g(i) = one - inc_gamma_series(a, x(i))
            else
                g(i) = inc_gamma_cf(a, x(i))
            end if
        end do
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

! ******************************************************************************
! BETA FUNCTIONS
! ------------------------------------------------------------------------------
    !> @brief Computes the incomplete beta function:
    !! I(a,b) = 1 / B(a,b) * integrate(t**(a - 1) * (1 - t)**(b - 1), t, 0, x)
    !!
    !! @param[in] a The parameter a.
    !! @param[in] b The parameter b.
    !! @param[in] x The parameter x.  This parameter must lie in the interval:
    !!  [0, 1].
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p x is not within its allowed 
    !!      range.
    function inc_beta_scalar(a, b, x, err) result(beta)
        ! Arguments
        real(dp), intent(in) :: a, b, x
        class(errors), intent(inout), optional, target :: err
        real(dp) :: beta

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: two = 2.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        if (x < zero .or. x > one) then
            call errmgr%report_error("inc_beta_scalar", &
                "The independent variable (x) must be >= 0, but not > 1.", &
                CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Process
        if (x == zero .or. x == one) then
            beta = zero
        else
            beta = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
                a * log(x) + b * log(one - x))
        end if
        if (x < (a + one) / (a + b + 2)) then
            beta = beta * inc_beta_cf(a, b, x) / a
        else
            beta = one - beta * inc_beta_cf(b, a, one - x) / b
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the incomplete beta function:
    !! I(a,b) = 1 / B(a,b) * integrate(t**(a - 1) * (1 - t)**(b - 1), t, 0, x)
    !!
    !! @param[in] a The parameter a.
    !! @param[in] b The parameter b.
    !! @param[in] x The parameter x.  This parameter must lie in the interval:
    !!  [0, 1].
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p x is not within its allowed 
    !!      range.
    function inc_beta_array(a, b, x, err) result(beta)
        ! Arguments
        real(dp), intent(in) :: a, b
        real(dp), intent(in), dimension(:) :: x
        class(errors), intent(inout), optional, target :: err
        real(dp), dimension(size(x)) :: beta

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: two = 2.0d0

        ! Local Variables
        logical :: check
        integer(i32) :: i, n
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(x)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        check = .true.
        do i = 1, n
            if (x(i) < zero .or. x(i) > one) then
                check = .false.
                exit
            end if
        end do
        if (.not.check) then
            call errmgr%report_error("inc_beta_array", &
                "The independent variable (x) must be >= 0, but not > 1.", &
                CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Process
        do i = 1, n
            if (x(i) == zero .or. x(i) == one) then
                beta(i) = zero
            else
                beta(i) = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
                    a * log(x(i)) + b * log(one - x(i)))
            end if
            if (x(i) < (a + one) / (a + b + 2)) then
                beta(i) = beta(i) * inc_beta_cf(a, b, x(i)) / a
            else
                beta(i) = one - beta(i) * inc_beta_cf(b, a, one - x(i)) / b
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    !> @brief Evaluates the incomplete beta function as a continued fraction.
    !!
    !! @param[in] a The parameter a.
    !! @param[in] b The parameter b.
    !! @param[in] x The independent variable.
    !!
    !! @return The result.
    function inc_beta_cf(a, b, x) result(beta)
        ! Arguments
        real(dp), intent(in) :: a, b, x
        real(dp) :: beta

        ! Parameters
        integer(i32), parameter :: itmax = 500
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, m, m2
        real(dp) :: eps, fpmin, aa, c, d, del, h, qab, qam, qap

        ! Initialization
        eps = epsilon(eps)
        fpmin = tiny(fpmin) / eps
        qab = a + b
        qap = a + one
        qam = a - one
        c = one
        d = one - qab * x / qap
        if (abs(d) < fpmin) d = fpmin
        d = one / d
        h = d
        do m = 1, itmax
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = one + aa * d
            if (abs(d) < fpmin) d = fpmin
            c = one + aa / c
            if (abs(c) < fpmin) c = fpmin
            d = one / d
            h = h * d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = one + aa * d
            if (abs(d) < fpmin) d = fpmin
            c = one + aa / c
            if (abs(c) < fpmin) c = fpmin
            d = one / d
            del = d * c
            h = h * del
            if (abs(del - one) < eps) exit
        end do
    end function

! ------------------------------------------------------------------------------
end module
