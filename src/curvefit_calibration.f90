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
!!
!! @par Example
!! The following example illustrates the use of a selection of the routines in
!! this module.
!! @code{.f90}
!! program example
!!     use iso_fortran_env
!!     use curvefit_calibration
!!     use curvefit_regression, only : linear_least_squares
!!     implicit none
!!
!!     ! Local Variables
!!     real(real64), parameter :: fullscale = 5.0d2
!!     real(real64), dimension(11) :: applied, output, measured, applied_copy
!!     real(real64) :: hyst, gain, nlin
!!     type(seb_results) :: s
!!
!!     ! Initialization
!!     applied = [0.0d0, 1.0d2, 2.0d2, 3.0d2, 4.0d2, 5.0d2, 4.0d2, 3.0d2, &
!!       2.0d2, 1.0d2, 0.0d0]
!!     output = [0.0d0, 0.55983d0, 1.11975d0, 1.67982d0, 2.24005d0, &
!!       2.80039d0, 2.24023d0, 1.68021d0, 1.12026d0, 0.56021d0, 0.00006d0]
!!     applied_copy = applied
!!
!!     ! Determine a suitable calibration gain (the least squares routine modifies
!!     ! applied; hence, the need for the copy)
!!     gain = linear_least_squares(output, applied_copy)
!!
!!     ! Apply the calibration gain
!!     measured = gain * output
!!
!!     ! Compute the SEB
!!     s = seb(applied, output, fullscale)
!!
!!     ! Compute the best fit nonlinearity
!!     nlin = nonlinearity(applied, measured)
!!
!!     ! Compute the hysteresis
!!     hyst = hysteresis(applied, measured)
!!
!!     ! Display the results
!!     print '(AF9.5)', "Calibration Gain: ", gain
!!     print '(AF6.4)', "SEB: ", s%seb
!!     print '(AF7.5)', "SEB Output: ", s%output
!!     print '(AF7.4)', "Best Fit Nonlinearity: ", nlin
!!     print '(AF6.4)', "Hysteresis: ", hyst
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Calibration Gain: 178.55935
!! SEB: 0.0518
!! SEB Output: 2.80010
!! Best Fit Nonlinearity: -0.0582
!! Hysteresis: 0.0911
!! @endcode
module curvefit_calibration
    use, intrinsic :: iso_fortran_env, only : int32, real64
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
    public :: return_to_zero
    public :: repeatability
    public :: crosstalk
    public :: split_ascend_descend

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Defines a container for static error band related information.
    type, bind(C) :: seb_results
        !> The static error band.
        real(real64) :: seb
        !> The static error band output, at full scale load.
        real(real64) :: output
        !> The slope of the static error band fit.
        real(real64) :: slope
    end type

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        function IDAMAX(n, dx, incx)
            use iso_fortran_env
            integer(int32), intent(in) :: n, incx
            real(real64), intent(in) :: dx(n)
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
    !> @brief Computes the return to zero error in an ascending/descending data
    !! set.
    interface return_to_zero
        module procedure :: rtz_1
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the repeatability of a sequence of tests.
    interface repeatability
        module procedure :: repeat_1
    end interface

! ------------------------------------------------------------------------------
    !> @brief Computes the crosstalk errors for a multiple degree-of-freedom
    !! data set.
    interface crosstalk
        module procedure :: xtalk_1
    end interface

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
        real(real64), intent(in), dimension(:) :: applied, output
        real(real64), intent(in) :: fullscale
        class(errors), intent(out), optional, target :: err
        type(seb_results) :: rst

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, j, npts
        real(real64), allocatable, dimension(:) :: ratio
        real(real64) :: arg, s, a, b, t, eps
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
        real(real64), intent(in), dimension(:) :: applied, measured
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, n
        real(real64) :: e
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
        real(real64), intent(in), dimension(:) :: applied, measured
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: factor = 1.0d-2

        ! Local Variables
        integer(int32) :: i, n, maxIndex, zeroIndex
        real(real64) :: zeroCheck, slope, e
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
        real(real64), intent(in), dimension(:) :: xascend, ascend, xdescend, descend
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, na, nd
        real(real64) :: delta
        real(real64), allocatable, dimension(:) :: y
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
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p applied and @p measured are not the
    !!      same size.
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
        real(real64), intent(in), dimension(:) :: applied, measured
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: n, na, nd, flag
        real(real64), allocatable, dimension(:) :: xa, ya, xd, yd
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
    !> @brief Computes the return to zero error in an ascending/descending data
    !! set.
    !!
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] measured An N-element array containing the calibrated output
    !!  of the instrument as a result of the values given in @p applied.
    !! @param[in] tol An optional input that specifies the tolerance used in
    !!  finding the matching data points.  If no value is specified, the default
    !!  value of the square root of machine precision times the largest
    !!  magnitude value in @p xcal is used.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p applied and @p measured are not the
    !!      same size.
    !!
    !! @return The return to zero error.
    function rtz_1(applied, measured, tol, err) result(rst)
        ! Arguments
        real(real64), intent(in), dimension(:) :: applied, measured
        real(real64), intent(in), optional :: tol
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, i1, i2, n
        real(real64) :: t, eps
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(applied)
        rst = zero
        eps = epsilon(eps)
        t = zero
        if (present(tol)) t = tol
        if (t <= eps) then
            i1 = IDAMAX(n, applied, 1)
            t = sqrt(eps) * abs(applied(i1))
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        if (size(measured) /= n) then
            call errmgr%report_error("rtz_1", "The measured data " // &
                "array must be the same size as the applied data array.", &
                CF_ARRAY_SIZE_ERROR)
            return
        end if

        ! Determine the indices of the matching values - base the search on the
        ! applied values, not the measured values
        i1 = 0
        i2 = 0
        do i = 1, n
            if (abs(applied(i)) < t) then
                i1 = i
                exit
            end if
        end do
        if (i1 == 0) return
        do i = i1 + 1, n
            if (abs(applied(i)) < t) then
                i2 = i
                exit
            end if
        end do
        if (i2 == 0) return

        ! Compute the return to zero error
        rst = measured(i2) - measured(i1)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the repeatability of a sequence of tests.
    !!
    !! @param[in] applied An NPTS-by-NTEST matrix containing at least 2 columns
    !!  (tests) of NPTS values applied to the measurement instrument.
    !! @param[in] measured An NPTS-by-NTEST matrix containing the corresponding
    !!  calibrated output from the instrument.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p applied and @p measured are not the
    !!      same size.
    !!
    !! @return The largest magnitude deviation from the initial test.
    !!
    !! @par Remarks
    !! Repeatability is considered as the largest magnitude deviation of
    !! subsequent tests from the initial test.  Noting that it is very likely
    !! that consecutive test points will vary slightly, test 2 through test N
    !! are linearly interpolated such that their test points line up with those
    !! from test 1.
    function repeat_1(applied, measured, err) result(rst)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: applied, measured
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, npts, ntests, ip
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        type(linear_interp) :: interp
        logical :: ascending
        real(real64), allocatable, dimension(:) :: y

        ! Initialization
        npts = size(applied, 1)
        ntests = size(applied, 2)
        rst = zero
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(measured, 1) /= npts .or. size(measured, 2) /= ntests) then
            call errmgr%report_error("repeat_1", "The measured data array " // &
                "must be the same size as the applied data array.", &
                CF_ARRAY_SIZE_ERROR)
            return
        end if

        ! Quick Return
        if (npts < 2 .or. ntests < 2) return

        ! Determine where the inflection point is within the data such that
        ! ascending and descending segments may be considered seperately
        ascending = applied(2,1) - applied(1,1) > zero
        ip = npts
        do i = 2, npts
            if (ascending) then
                if (applied(i,1) < applied(i-1,1)) then
                    ip = i - 1
                    exit
                end if
            else
                if (applied(i,1) > applied(i-1,1)) then
                    ip = i - 1
                    exit
                end if
            end if
        end do

        ! Cycle over each test, and use the initial test (column 1) as the
        ! reference
        do i = 2, ntests
            ! Determine how far subsequent tests move from the initial test
            ! for the ascending data
            call interp%initialize(applied(1:ip,i), measured(1:ip,i), &
                err = errmgr)
            if (errmgr%has_error_occurred()) return
            y = abs(interp%interpolate(applied(1:ip,1)) - measured(1:ip,1))
            rst = max(rst, maxval(y))

            ! Consider the descending data
            if (ip < npts) then
                call interp%initialize(applied(ip:npts,i), &
                    measured(ip:npts,i), err = errmgr)
                if (errmgr%has_error_occurred()) return
                y = abs(interp%interpolate(applied(ip:npts,1)) - &
                    measured(ip:npts,1))
                rst = max(rst, maxval(y))
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the crosstalk errors for a multiple degree-of-freedom
    !! data set.
    !!
    !! @param[in] xerr An NPTS-by-NDOF matrix containing the measurement error
    !!  values (computed such that XERR = X MEASURED - X APPLIED).
    !! @param[in] indices A 2*NDOF element array containing row indices defining
    !!  the rows where each degree-of-freedom was applied in the data set
    !!  @p xerr.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p indices is not 2*NDOF in size.
    !!  - CF_ARRAY_INDEX_ERROR: Occurs if any of the entries in @p indices are
    !!      outside the row bounds of @p xerr.
    !!
    !! @return A NDOF-by-NDOF matrix containing the crosstalk errors such that
    !!  each loaded degree of freedom is represented by its own row, and each
    !!  responding degree of freedom is represented by its own column.
    !!
    !! @par Usage
    !! The following program computes the crosstalk errors for a 2 DOF system.
    !! The applied data is as follows (there are 34 data points for each DOF).
    !! @verbatim
    !! 0	            0
    !! 3000	        0
    !! 6000	        0
    !! 7500	        0
    !! 9000	        0
    !! 12000	        0
    !! 15000	        0
    !! 7500	        0
    !! 0	            0
    !! 0	            0
    !! -3000	        0
    !! -6000	        0
    !! -7500	        0
    !! -9000	        0
    !! -12000	        0
    !! -15000	        0
    !! -7500	        0
    !! 0	            0
    !! --------------------
    !! 0	            0
    !! 0	            67.79087067
    !! 0	            135.5817413
    !! 0	            203.3726196
    !! 0	            271.1634827
    !! 0	            338.9543762
    !! 0	            203.3726196
    !! 0	            0
    !! 0	            0
    !! 0	            -67.79087067
    !! 0	            -135.5817413
    !! 0	            -203.3726196
    !! 0	            -271.1634827
    !! 0	            -338.9543762
    !! 0	            -203.3726196
    !! 0	            0
    !! @endverbatim
    !! The data output from the instrument under test is as follows.
    !! @verbatim
    !! 0	                0
    !! 0.389050007	    1.22E-03
    !! 0.778159976	    2.59E-03
    !! 0.972689986	    2.90E-03
    !! 1.167140007	    3.14E-03
    !! 1.555999994	    3.38E-03
    !! 1.944839954	    3.56E-03
    !! 0.972599983	    4.77E-03
    !! -1E-05	            -1.00E-05
    !! 0	                0
    !! -0.388886005	    2.10E-04
    !! -0.777750015	    5.10E-04
    !! -0.972150028	    6.90E-04
    !! -1.166540027	    8.80E-04
    !! -1.555330038	    1.30E-03
    !! -1.944100022	    1.78E-03
    !! -0.971710026	    5.80E-04
    !! 4.00E-05	        3.00E-05
    !! ----------------------------
    !! 0	                0
    !! -4.40E-04	        0.271560013
    !! -1.30E-03	        0.543290019
    !! -2.40E-03	        0.815069973
    !! -3.82E-03	        1.086820006
    !! -5.28E-03	        1.358809948
    !! -2.57E-03	        0.815530002
    !! 1.50E-04	        1.00E-05
    !! 0	                0
    !! 1.44E-03	        -0.271450013
    !! 3.06E-03	        -0.543120027
    !! 4.46E-03	        -0.814930022
    !! 5.67E-03	        -1.086799979
    !! 6.88E-03	        -1.35879004
    !! 4.51E-03	        -0.815479994
    !! -2.00E-05	        0
    !! @endverbatim
    !! The code to compute the crosstalk errors given the above raw data is then
    !! as follows.
    !! @code{.f90}
    !! program main
    !!   ! Parameters
    !!   integer(int32), parameter :: npts = 34
    !!   integer(int32), parameter :: ndof = 2
    !!
    !!   ! Local Variables
    !!   integer(int32) :: indices(2*ndof)
    !!   real(real64), dimension(npts, ndof) :: xin, xout, xerr, xmeas
    !!   real(real64), dimension(ndof, npts) :: xint, xmeast
    !!   real(real64), dimension(ndof, ndof) :: c, ans, xt
    !!
    !!   ! Initialization
    !!   xin = reshape([0.0, 3000.0, 6000.0, 7500.0, 9000.0, 12000.0, &
    !!       15000.0, 7500.0, 0.0, 0.0, -3000.0, -6000.0, -7500.0, -9000.0, &
    !!       -12000.0, -15000.0, -7500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
    !!       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
    !!       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
    !!       0.0, 0.0, 0.0, 67.7908728, 135.5817456, 203.3726184, 271.1634912, &
    !!       338.954364, 203.3726184, 0.0, 0.0, -67.7908728, -135.5817456, &
    !!       -203.3726184, -271.1634912, -338.954364, -203.3726184, 0.0], &
    !!       [npts, ndof])
    !!   xout = reshape([0.0, 0.38905, 0.77816, 0.97269, 1.16714, 1.556, &
    !!       1.94484, 0.9726, -1.0e-5, 0.0, -0.388886, -0.77775, -0.97215, &
    !!       -1.16654, -1.55533, -1.9441, -0.97171, 4.0e-5, 0.0, -0.00044, &
    !!       -0.0013, -0.0024, -0.00382, -0.00528, -0.00257, 0.00015, 0.0, &
    !!       0.00144, 0.00306, 0.00446, 0.00567, 0.00688, 0.00451, -2.0e-5, &
    !!       0.0, 0.00122, 0.00259, 0.0029, 0.00314, 0.00338, 0.00356, 0.00477,&
    !!       -1.0e-5, 0.0, 0.00021, 0.00051, 0.00069, 0.00088, 0.0013, 0.00178,&
    !!       0.00058, 3.0e-5, 0.0, 0.27156, 0.54329, 0.81507, 1.08682, 1.35881,&
    !!       0.81553, 1.0e-5, 0.0, -0.27145, -0.54312, -0.81493, -1.0868, &
    !!       -1.35879, -0.81548, 0.0], [npts, ndof])
    !!
    !!   ! Compute the calibration gains
    !!   xint = transpose(xin)
    !!   xmeast = transpose(xout)
    !!   c = linear_least_squares(xmeast, xint)
    !!   xmeas = matmul(xout, transpose(c))
    !!   xerr = xmeas - xin
    !!
    !!   ! The indices are
    !!   indices = [1, 17, 18, 34]
    !!
    !!   ! Compute the crosstalk matrix
    !!   xt = crosstalk(xerr, indices)
    !! end program
    !! @endcode
    !! The least squares fit generates the following matrix of calibration
    !! gains.
    !! @verbatim
    !! 7713.710427	    33.5206917
    !! -0.214743728	    249.4768498
    !! @endverbatim
    !! The resulting measured values are then as follows.
    !! @verbatim
    !! 0	                0
    !! 3001.05999	        0.220815702
    !! 6002.58754	        0.479040049
    !! 7503.146099	    0.514603781
    !! 9003.085297	    0.532721294
    !! 12002.64668	    0.509090486
    !! 15002.05157	    0.470495427
    !! 7502.514526	    0.981144827
    !! -0.077472309	    -0.002492621
    !! 0	                0
    !! -2999.74699	    0.135900968
    !! -5999.321307	    0.294250127
    !! -7498.860677	    0.380902151
    !! -8998.322469	    0.470046784
    !! -11997.32196	    0.658317276
    !! -14996.16495	    0.861552092
    !! -7495.47032	    0.353365205
    !! 0.30955403	        7.48E-03
    !! ----------------------------
    !! 0	                0
    !! 5.708846869	    67.74803112
    !! 8.183633648	    135.5385616
    !! 8.808803389	    203.3416047
    !! 6.96458466	        271.1372518
    !! 4.81985707	        338.9927591
    !! 7.512893885	    203.4564077
    !! 1.157391826	    2.46E-03
    !! 0	                0
    !! 2.008550989	    -67.72080332
    !! 5.398195074	    -135.4965304
    !! 7.086130239	    -203.3071324
    !! 7.306450061	    -271.1326527
    !! 7.522743936	    -338.9881361
    !! 7.453379661	    -203.4443484
    !! -0.154274205	    4.29E-06
    !! @endverbatim
    !! The crosstalk error matrix is then as follows.
    !! @verbatim
    !! 0	                0.981144827
    !! 8.808803389	    0
    !! @endverbatim
    function xtalk_1(xerr, indices, err) result(xt)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: xerr
        integer(int32), intent(in), dimension(:) :: indices
        real(real64), dimension(size(xerr, 2), size(xerr, 2)) :: xt
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, j, n, npts, ndof, first, last, ind
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        npts = size(xerr, 1)
        ndof = size(xerr, 2)
        n = size(indices) ! must be even
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        xt = zero

        ! Input Check
        if (size(indices) /= 2 * ndof) then
            call errmgr%report_error("xtalk_1", "The indices array is " // &
                "not compatible in size with the supplied data array.", &
                CF_ARRAY_SIZE_ERROR)
            return
        end if
        if (maxval(indices) > npts .or. minval(indices) < 1) then
            call errmgr%report_error("xtalk_1", "There are entries in " // &
                "indices that exceed the bounds of the input data array.", &
                CF_ARRAY_INDEX_ERROR)
            return
        end if


        ! Process
        do i = 1, ndof ! Cycle over each DOF
            ! Locate the starting and stopping row indices for the DOF
            first = indices(2 * i - 1)
            last = indices(2 * i)

            ! Compute the error for each DOF over the specified range
            do j = 1, ndof
                if (i /= j) then
                    npts = last - first + 1
                    ind = IDAMAX(npts, xerr(first:last,j), 1) + first - 1
                    xt(i,j) = xerr(ind,j)
                else
                    xt(i,j) = zero
                end if
            end do
        end do
    end function

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
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: ascend, descend
        integer(int32), intent(out) :: nascend, ndescend
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: tol = 0.5d-2

        ! Local Variables
        logical :: is_ascend, check
        integer(int32) :: i, i1, i2, n
        real(real64) :: t
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
