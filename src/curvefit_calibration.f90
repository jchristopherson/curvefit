! curvefit_calibration.f90

!> @brief \b curvefit_calibration
!!
!! @par Purpose
!! To provide routines for computing calibration performance metrics commonly 
!! used to assess the fitness of a calibration curve fit.
!!
!! @par References
!! - Wheeler, Anthony J., Ganji, Ahmad R., "Introduction to Engineering 
!!      Experimentation," Third Edition, Prentice Hall.
module curvefit_calibration
    use curvefit_core
    use ferror, only : errors
    use curvefit_interp, only : linear_interp
    implicit none
    private
    public :: seb_results
    public :: seb
    public :: nonlinearity
    public :: terminal_nonlinearity
    public :: hysteresis
    public :: split_ascend_descend

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Defines a container for static error band related information.
    type, bind(C) :: seb_results
        !> The static error band.
        real(dp) :: seb
        !> The static error band output, at full scale load.
        real(dp) :: output
        !> The slope of the static error band fit.
        real(dp) :: slope
    end type

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        function IDAMAX(n, dx, incx)
            use curvefit_core, only : dp, i32
            integer(i32), intent(in) :: n, incx
            real(dp), intent(in) :: dx(n)
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the static error band of a data set.
    interface seb
        module procedure :: seb_1
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the best-fit nonlinearity of a data set.
    interface nonlinearity
        module procedure :: bf_nonlin
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the terminal nonlinearity of a data set.
    interface terminal_nonlinearity
        module procedure :: term_nonlin
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the hysteresis in an ascending/descending data set.
    interface hysteresis
        module procedure :: hysteresis_1
        module procedure :: hysteresis_2
    end interface

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
    !> @brief Splits a data set into ascending and descending components.
    interface split_ascend_descend
        module procedure :: split_ascend_descend_1
    end interface

contains
! ------------------------------------------------------------------------------
    !> @brief Computes the static error band of a data set.
    !!
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] output An N-element array containing the values output by
    !!  the instrument as a result of the values given in @p applied.
    !! @param[in] fullscale The full scale measurement value for the instrument.
    !!  The units must be consistent with those of @p applied.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p applied and @p output are not the
    !!      same size.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p fullscale is sufficiently close
    !!      to zero to be considered zero.  Sufficiently close in this instance
    !!      is considered to be the square root of machine precision.
    !!
    !! @return The static error band information.
    function seb_1(applied, output, fullscale, err) result(rst)
        ! Arguments
        real(dp), intent(in), dimension(:) :: applied, output
        real(dp), intent(in) :: fullscale
        class(errors), intent(out), optional, target :: err
        type(seb_results) :: rst

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: i, j, npts
        real(dp), allocatable, dimension(:) :: ratio
        real(dp) :: arg, s, a, b, t, eps
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        eps = epsilon(eps)
        npts = size(applied)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(output) /= npts) then
            call errmgr%report_error("seb_1", "The output data array " // &
                "must be the same size as the applied data array.", &
                CF_ARRAY_SIZE_ERROR)
            return
        end if
        if (abs(fullscale) < sqrt(eps)) then
            call errmgr%report_error("seb_1", &
                "The full scale value is too close to zero to be valid.", &
                CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Process
        ratio = applied / fullscale
        t = zero
        b = zero
        do i = 1, npts - 1
            do j = 1, npts - 1
                if (i /= j) then
                    arg = ratio(j) + ratio(i)
                    if (arg /= zero) then
                        s = (output(j) + output(i)) / arg
                        a = abs((output(j) - s * ratio(j)) / s)
                    else
                        s = zero
                        a = zero
                    end if
                    if (a >= b) then
                        b = a
                        t = s
                    end if
                end if
            end do
        end do

        ! Collect the output
        rst%seb = b * fullscale
        rst%output = t
        rst%slope = fullscale / t
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the best-fit nonlinearity of a data set.
    !!
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] measured An N-element array containing the calibrated output
    !!  of the instrument as a result of the values given in @p applied.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p applied and @p measured are not the
    !!      same size.
    !!
    !! @return The nonlinearity error.
    function bf_nonlin(applied, measured, err) result(rst)
        ! Arguments
        real(dp), intent(in), dimension(:) :: applied, measured
        class(errors), intent(inout), optional, target :: err
        real(dp) :: rst

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: i, n
        real(dp) :: e
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        e = zero
        rst = zero
        n = size(applied)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(measured) /= n) then
            call errmgr%report_error("bf_nonlin", "The measured data " // &
                "array must be the same size as the applied data array.", &
                CF_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do i = 1, n
            e = measured(i) - applied(i)
            if (abs(e) > abs(rst)) rst = e
        end do
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the terminal nonlinearity of a data set.
    !!
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] measured An N-element array containing the calibrated output
    !!  of the instrument as a result of the values given in @p applied.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p applied and @p measured are not the
    !!      same size.
    !!
    !! @return The terminal nonlinearity error.
    function term_nonlin(applied, measured, err) result(rst)
        ! Arguments
        real(dp), intent(in), dimension(:) :: applied, measured
        class(errors), intent(inout), optional, target :: err
        real(dp) :: rst

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: factor = 1.0d-2

        ! Local Variables
        integer(i32) :: i, n, maxIndex, zeroIndex
        real(dp) :: zeroCheck, slope, e
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        e = zero
        rst = zero
        n = size(applied)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(measured) /= n) then
            call errmgr%report_error("term_nonlin", "The measured data " // &
                "array must be the same size as the applied data array.", &
                CF_ARRAY_SIZE_ERROR)
            return
        end if

        ! Locate the largest magnitude and initial zero values
        maxIndex = IDAMAX(n, applied, 1)
        zeroCheck = factor * abs(applied(maxIndex))
        zeroIndex = 1
        do i = 1, n
            if (abs(applied(i)) < zeroCheck) then
                zeroIndex = i
                exit
            end if
        end do

        ! Compute the slope of the line from zero to the terminal point
        slope = (applied(maxIndex) - applied(zeroIndex)) / &
            (measured(maxIndex) - measured(zeroIndex))

        ! Compute the nonlinearity error
        do i = 1, n
            e = slope * measured(i) - applied(i)
            if (abs(e) > abs(rst)) rst = e
        end do
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the hysteresis in an ascending/descending data set.
    !!
    !! @param[in] xascend An N-element array containing the ascending
    !!  calibration points.  This array must be monotonically increasing or
    !!  decreasing.
    !! @param[in] ascend An N-element array containing the sensor output to
    !!  the calibration points in @p xascend.
    !! @param[in] xdescend An M-element array containing the descending
    !!  calibration points.  This array must be monotonically increasing or
    !!  decreasing.
    !! @param[in] descend An M-element array containing the sensor output to
    !!  the calibration points in @p xdescend.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p xascend and @p ascend are not the
    !!      same size, or if @p xdescend and @p descend are not the same size.
    !!  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if the calibration data is not
    !!      monotonic in nature (either ascending or descending).
    !!
    !! @return The hysteresis error.
    !!
    !! @par Remarks
    !! In order to account for slight variations between similar ascending and
    !! descending points, the algorithm used performs a linear interpolation
    !! between data points.  The resulting interpolated value is then used to
    !! compute the reported hysteresis error.
    function hysteresis_1(xascend, ascend, xdescend, descend, err) result(rst)
        ! Arguments
        real(dp), intent(in), dimension(:) :: xascend, ascend, xdescend, descend
        class(errors), intent(inout), optional, target :: err
        real(dp) :: rst

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: i, na, nd
        real(dp) :: delta
        real(dp), allocatable, dimension(:) :: y
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        type(linear_interp) :: interp

        ! Initialization
        rst = zero
        na = size(xascend)
        nd = size(xdescend)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        if (size(ascend) /= na) then
            call errmgr%report_error("hysteresis_1", "The ascending data " // &
                "array must be the same size as the applied ascending array.", &
                CF_ARRAY_SIZE_ERROR)
            return
        end if
        if (size(descend) /= nd) then
            call errmgr%report_error("hysteresis_1", "The descending data " // &
                "array must be the same size as the applied descending array.",&
                CF_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        if (na >= nd) then
            ! Interpolate on the ascending data
            call interp%initialize(xascend, ascend, err = errmgr)
            if (errmgr%has_error_occurred()) return
            y = interp%interpolate(xdescend)

            ! Compute the delta between descend and y
            do i = 1, nd
                delta = abs(descend(i) - y(i))
                rst = max(delta, rst)
            end do
        else
            ! Interpolate on the descending data
            call interp%initialize(xdescend, descend, err = errmgr)
            if (errmgr%has_error_occurred()) return
            y = interp%interpolate(xascend)

            ! Compute the delta between ascend and y
            do i = 1, na
                delta = abs(ascend(i) - y(i))
                rst = max(delta, rst)
            end do
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the hysteresis in an ascending/descending data set.
    !!
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] measured An N-element array containing the calibrated output
    !!  of the instrument as a result of the values given in @p applied.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if the calibration data is not
    !!      monotonic in nature (either ascending or descending).
    !!
    !! @return The hysteresis error.
    !!
    !! @par Remarks
    !! In order to account for slight variations between similar ascending and
    !! descending points, the algorithm used performs a linear interpolation
    !! between data points.  The resulting interpolated value is then used to
    !! compute the reported hysteresis error.
    function hysteresis_2(applied, measured, err) result(rst)
        ! Arguments
        real(dp), intent(in), dimension(:) :: applied, measured
        class(errors), intent(inout), optional, target :: err
        real(dp) :: rst

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: n, na, nd, flag
        real(dp), allocatable, dimension(:) :: xa, ya, xd, yd
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        rst = zero
        n = size(applied)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(measured) /= n) then
            call errmgr%report_error("hysteresis_2", "The measured data " // &
                "array must be the same size as the applied data array.", &
                CF_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(xa(n), stat = flag)
        if (flag == 0) allocate(ya(n), stat = flag)
        if (flag == 0) allocate(xd(n), stat = flag)
        if (flag == 0) allocate(yd(n), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("hysteresis_2", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Split the input data set into ascending and descending parts
        call split_ascend_descend(applied, xa, xd, na, nd)

        ! Ensure there's descending data; else, there's no point in continuing
        if (nd == 0) return

        ! Copy over the sensor output data
        ya = measured(1:na)
        yd = measured(na:na+nd-1)

        ! Compute the hysteresis
        rst = hysteresis_1(xa(1:na), ya, xd(1:nd), yd, err)
    end function

! ------------------------------------------------------------------------------
    ! REPEATABILITY

! ------------------------------------------------------------------------------
    ! CROSSTALK
    
! ------------------------------------------------------------------------------

! ******************************************************************************
! SUPPORTING ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Splits a data set into ascending and descending components.
    !!
    !! @param[in] x An N-element array containing the data set to split.
    !! @param[out] ascend An array where the ascending points will be written.
    !!  Ensure this array is appropriately sized to accept all the ascending
    !!  points (it can be oversized).
    !! @param[out] descend An array where the descending points will be written.
    !!  Ensure this array is appropriately sized to accept all the descending
    !!  points (it can be oversized).
    !! @param[out] nascend The actual number of values written into @p ascend.
    !! @param[out] ndescend The actual number of values written into @p descend.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if either @p ascend or @p descend is
    !!      too small to actually accept all of the necessary data.
    !!
    !! @par Remarks
    !! The routine operates by finding the first occurrence where the data set
    !! is no longer monotonic, and then copies everything prior to that
    !! value, along with the the inflection value, into the output ascending
    !! data array.  The routine then searches for either a change in direction,
    !! or a value that matches the first value in the ascending data set within
    !! some tolerance to determine the bounds on the descending data set.  Once
    !! the bounds are determined, the descending data set is copied from the
    !! original array and placed in the output descending data array.  This
    !! then means that any remaining data in the original data set that lies
    !! after either of the aforementioned sets is ignored.
    !!
    !! @par Example
    !! @verbatim
    !! Given the following array X,
    !!  X:
    !!   0.0000000000000000
    !!   0.38905000686645508
    !!   0.77815997600555420
    !!   0.97268998622894287
    !!   1.1671400070190430
    !!   1.5559999942779541
    !!   1.9448399543762207
    !!   0.97259998321533203
    !!   -9.9999997473787516E-006
    !!
    !! This routine splits the array into the following ascending and
    !! descending arrays.
    !!
    !! ASCENDING:
    !!   0.0000000000000000
    !!   0.38905000686645508
    !!   0.77815997600555420
    !!   0.97268998622894287
    !!   1.1671400070190430
    !!   1.5559999942779541
    !!   1.9448399543762207
    !!
    !! DESCENDING:
    !!   1.9448399543762207
    !!   0.97259998321533203
    !!   -9.9999997473787516E-006
    !! @endverbatim
    subroutine split_ascend_descend_1(x, ascend, descend, nascend, ndescend, &
        err)
        ! Arguments
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: ascend, descend
        integer(i32), intent(out) :: nascend, ndescend
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: tol = 0.5d-2

        ! Local Variables
        logical :: is_ascend, check
        integer(i32) :: i, i1, i2, n
        real(dp) :: t
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        n = size(x)
        nascend = 0
        ndescend = 0
        if (n >= 2) then
            is_ascend = x(2) - x(1) > zero
        else
            is_ascend = .false.
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Determine where the data is no longer monotonic
        do i = 1, n
            if (i < n) then
                check = x(i+1) - x(i) < zero
            else
                check = is_ascend
            end if
            if (check .eqv. is_ascend) then
                nascend = i
                exit
            end if
        end do

        ! Ensure the receiving array is appropriately sized, and then copy over
        ! the correct items from X.
        if (size(ascend) < nascend) then
            write(errmsg, '(AI0AI0A)') "The receiving array for the " // &
                "ascending data set is too small.  The minimum required " // &
                "dimension is: ", nascend, ", but the current size is: ", &
                size(ascend), "."
            call errmgr%report_error("split_ascend_descend_1", trim(errmsg), &
                CF_ARRAY_SIZE_ERROR)
            return
        end if
        ascend(1:nascend) = x(1:nascend)

        ! If NASCEND == N, then the data set only included ascending data.
        if (nascend == n) then
            return
        end if

        ! As we now have the apex of the data set, scale the tolerance value
        ! based upon the largest magnitude value.  The tolerance value will be
        ! used to determine when we've returned back to the original starting
        ! point (useful for data sets that cross zero).
        t = tol * abs(x(nascend))

        ! Assuming we're here, it's now safe to start the descending search.
        ! Search until either we get to the same value, within T, of where we
        ! started, or the end of the data set.  Whichever comes first.  The
        ! other check is to determine if there is a switch in direction in the
        ! data.
        !i1 = nascend + 1
        i1 = nascend
        i2 = i1
        do i = i1, n
            if (i < n) then
                check = x(i+1) - x(i) > zero
            else
                check = .not.is_ascend
            end if
            if (abs(x(i) - x(1)) < t .or. check .eqv. is_ascend) then
                i2 = i
                exit
            end if
        end do
        if (i2 == i1) i2 = n
        ndescend = i2 - i1 + 1

        ! Ensure DESCEND is appropriately sized
        if (size(descend) < ndescend) then
            write(errmsg, '(AI0AI0A)') "The receiving array for the " // &
                "descending data set is too small.  The minimum required " // &
                "dimension is: ", nascend, ", but the current size is: ", &
                size(ascend), "."
            call errmgr%report_error("split_ascend_descend_1", trim(errmsg), &
                CF_ARRAY_SIZE_ERROR)
            return
        end if
        descend(1:ndescend) = x(i1:i2)
    end subroutine

! ------------------------------------------------------------------------------
end module
