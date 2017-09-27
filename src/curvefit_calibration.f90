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
    implicit none
    private
    public :: seb_results
    public :: seb
    public :: nonlinearity
    public :: terminal_nonlinearity

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
    ! HYSTERESIS

! ------------------------------------------------------------------------------
    ! REPEATABILITY

! ------------------------------------------------------------------------------
    ! CROSSTALK
    
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
