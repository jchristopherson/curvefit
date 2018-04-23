! curvefit_c_binding.f90

!> @brief \b curvefit_c_binding
!!
!! @par Purpose
!! Provides C bindings to the curvefit library.
module curvefit_c_binding
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use curvefit_core
    use curvefit_interp
    use curvefit_statistics
    use curvefit_regression
    use curvefit_calibration
    use ferror, only : errors
    use ferror_c_binding, only : errorhandler, get_errorhandler
    use nonlin_types, only : iteration_behavior
    use nonlin_c_binding, only : solver_control
    implicit none

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        !> @brief Describes a routine for finding the coefficients of a function
        !! of one variable.
        !!
        !! @param[in] x The independent variable.
        !! @param[in] n The number of coefficients in @p c.
        !! @param[in] c An array of function coefficients.
        !!
        !! @result The value of the function at @p x.
        function creg_fcn(x, n, c) result(f)
            use curvefit_core, only : dp, i32
            real(real64), intent(in), value :: x
            integer(int32), intent(in), value :: n
            real(real64), intent(in) :: c(n)
            real(real64) :: f
        end function
    end interface

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief A C compatible type encapsulating a linear_interp object.
    type, bind(C) :: c_linear_interp
        !> @brief A pointer to the linear_interp object.
        type(c_ptr) :: ptr
        !> @brief The size of the linear_interp object, in bytes.
        integer(int32) :: n
    end type

! ------------------------------------------------------------------------------
    !> @brief A C compatible type encapsulating a polynomial_interp object.
    type, bind(C) :: c_polynomial_interp
        !> @brief A pointer to the polynomial_interp object.
        type(c_ptr) :: ptr
        !> @brief The size of the polynomial_interp object, in bytes.
        integer(int32) :: n
    end type

! ------------------------------------------------------------------------------
    !> @brief A C compatible type encapsulating a spline_interp object.
    type, bind(C) :: c_spline_interp
        !> @brief A pointer to the spline_interp object.
        type(c_ptr) :: ptr
        !> @brief The size of the spline_interp object, in bytes.
        integer(int32) :: n
    end type

! ------------------------------------------------------------------------------
    !> @brief A C compatible type encapsulating a lowess_smoothing object.
    type, bind(C) :: c_lowess_smoothing
        !> @brief A pointer to the lowess_smoothing object.
        type(c_ptr) :: ptr
        !> @brief The size of the lowess_smoothing object, in bytes.
        integer(int32) :: n
    end type

! ------------------------------------------------------------------------------
    !> @brief A C compatible type encapsulating a nonlinear_regression object.
    type, bind(C) :: c_nonlinear_regression
        !> @brief A pointer to the nonlinear_regression object.
        type(c_ptr) :: ptr
        !> @brief The size of the nonlinear_regression object, in bytes.
        integer(int32) :: n
    end type

! ------------------------------------------------------------------------------
    !> @brief A type for helping to interface between a C function pointer, and
    !! the nonlinear_regression type.
    type, extends(nonlinear_regression) :: cnonlin_reg_helper
        private
        !> A pointer to the target creg_fcn routine.
        procedure(creg_fcn), pointer, nopass :: m_cfcn => null()
    contains
        !> @brief Executes the routine containing the function to evaluate.
        procedure, public :: fcn => crh_fcn
        !> @brief Tests if the pointer to the function containing the equation
        !! to solve has been assigned.
        procedure, public :: is_fcn_defined => crh_is_fcn_defined
        !> @brief Establishes a pointer to the routine containing the equations
        !! to solve.
        procedure, public :: set_cfcn => crh_set_fcn
    end type

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

contains
! ******************************************************************************
! CURVEFIT_CORE MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Tests to see if an array is montonically increasing or decreasing.
    !!
    !! @param[in] n The number of elements in the array.
    !! @param[in] x The array to test.
    !!
    !! @return Returns true if @p x is monotonic; else, false.
    pure function is_monotonic_c(n, x) result(rst) &
            bind(C, name = "is_monotonic")
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        logical(c_bool) :: rst
        rst = is_monotonic(x)
    end function

! ******************************************************************************
! LINEAR INTERPOLATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Retrieves the linear_interp object from the C compatible
    !! c_linear_interp data structure.
    !!
    !! @param[in] obj The C compatible c_linear_interp object.
    !! @param[out] li A pointer to the resulting linear_interp object.  This
    !!  pointer can be NULL dependent upon the state of @p obj.
    subroutine get_linear_interp(obj, li)
        ! Arguments
        type(c_linear_interp), intent(in), target :: obj
        type(linear_interp), intent(out), pointer :: li

        ! Local Variables
        type(c_ptr) :: testptr

        ! Process
        testptr = c_loc(obj)
        nullify(li)
        if (.not.c_associated(testptr)) return
        if (.not.c_associated(obj%ptr)) return
        if (obj%n == 0) return
        call c_f_pointer(obj%ptr, li)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Initializes a new c_linear_interp object.
    !!
    !! @param[out] obj The c_linear_interp object to initialize.
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the x-components of each data
    !!  point.  This array must be monotonic (ascending or descending only).
    !! @param[in] y An N-element array containing the y-components of each data
    !!  point.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
    !!      increasing or decreasing.
    subroutine lininterp_init_c(obj, n, x, y, err) &
            bind(C, name = "alloc_linear_interp")
        ! Arguments
        type(c_linear_interp), intent(out) :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n), y(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr
        type(linear_interp), pointer :: li

        ! Process
        allocate(li)
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            call li%initialize(x, y, 1, eptr)
        else
            call li%initialize(x, y, 1)
        end if
        obj%ptr = c_loc(li)
        obj%n = sizeof(li)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Frees resources held by a c_linear_interp object.
    !!
    !! @param[in,out] obj The c_linear_interp object.
    subroutine lininterp_free_c(obj) bind(C, name = "free_linear_interp")
        ! Arguments
        type(c_linear_interp), intent(inout), target :: obj

        ! Local Variables
        type(linear_interp), pointer :: li

        ! Process
        call get_linear_interp(obj, li)
        if (associated(li)) deallocate(li)
        obj%n = 0
        obj%ptr = c_null_ptr
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Performs a linear interpolation to determine the points @p y that
    !! for the requested indendent variable values in @p x.
    !!
    !! @param[in] obj The c_linear_interp object.
    !! @param[in] n The number of points to interpolate.
    !! @param[in] x An N-element array containing the values of the independent
    !!  variable at which to interpolate.
    !! @param[out] y An N-element array where the interpolated values can be
    !!  written.
    subroutine lininterp_interp_c(obj, n, x, y) &
            bind(C, name = "linear_interpolate")
        ! Arguments
        type(c_linear_interp), intent(in), target :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64), intent(out) :: y(n)

        ! Local Variables
        type(linear_interp), pointer :: li

        ! Process
        call get_linear_interp(obj, li)
        if (.not.associated(li)) return
        y = li%interpolate(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the number of points used by the interpolation object.
    !!
    !! @param[in] obj The c_linear_interp object.
    !!
    !! @return The number of points.
    function lininterp_get_pt_count_c(obj) result(n) &
            bind(C, name = "linear_interp_get_point_count")
        ! Arguments
        type(c_linear_interp), intent(in), target :: obj
        integer(int32) :: n

        ! Local Variables
        type(linear_interp), pointer :: li

        ! Process
        n = 0
        call get_linear_interp(obj, li)
        if (.not.associated(li)) return
        n = li%get_count()
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets a copy of the data points stored by the interpolation object.
    !!
    !! @param[in] obj The c_linear_interp object.
    !! @param[in] n The size of the buffer arrays.
    !! @param[out] x An N-element array where the x-coordinate data will be
    !!  written.
    !! @param[out] y An N-element array where the y-coordinate data will be
    !!  written.
    !!
    !! @par Remarks
    !! If @p n is different than the actual number of points that exist, the
    !! lesser of the two values will be utilized.  The interpolation object
    !! can be queried to determine the quantity of stored points.
    subroutine lininterp_get_pts_c(obj, n, x, y) &
            bind(C, name = "linear_interp_get_points")
        ! Arguments
        type(c_linear_interp), intent(in), target :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(out) :: x(n), y(n)

        ! Local Variables
        integer(int32) :: i, npts
        type(linear_interp), pointer :: li

        ! Process
        call get_linear_interp(obj, li)
        if (.not.associated(li)) return
        npts = min(n, li%get_count())
        do concurrent (i = 1:npts)
            x(i) = li%get_x(i)
            y(i) = li%get_y(i)
        end do
    end subroutine

! ******************************************************************************
! POLYNOMIAL INTERPOLATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Retrieves the polynomial_interp object from the C compatible
    !! c_polynomial_interp data structure.
    !!
    !! @param[in] obj The C compatible c_polynomial_interp object.
    !! @param[out] interp A pointer to the resulting polynomial_interp object.
    !!  This pointer can be NULL dependent upon the state of @p obj.
    subroutine get_polynomial_interp(obj, interp)
        ! Arguments
        type(c_polynomial_interp), intent(in), target :: obj
        type(polynomial_interp), intent(out), pointer :: interp

        ! Local Variables
        type(c_ptr) :: testptr

        ! Process
        testptr = c_loc(obj)
        nullify(interp)
        if (.not.c_associated(testptr)) return
        if (.not.c_associated(obj%ptr)) return
        if (obj%n == 0) return
        call c_f_pointer(obj%ptr, interp)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Initializes a new c_polynomial_interp object.
    !!
    !! @param[out] obj The c_polynomial_interp object to initialize.
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the x-components of each data
    !!  point.  This array must be monotonic (ascending or descending only).
    !! @param[in] y An N-element array containing the y-components of each data
    !!  point.
    !! @param[in] order The order of the interpolating polynomial.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
    !!      increasing or decreasing.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
    subroutine polyinterp_init_c(obj, n, x, y, order, err) &
            bind(C, name = "alloc_polynomial_interp")
        ! Arguments
        type(c_polynomial_interp), intent(out) :: obj
        integer(int32), intent(in), value :: n, order
        real(real64), intent(in) :: x(n), y(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr
        type(polynomial_interp), pointer :: interp

        ! Process
        allocate(interp)
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            call interp%initialize(x, y, order, eptr)
        else
            call interp%initialize(x, y, order)
        end if
        obj%ptr = c_loc(interp)
        obj%n = sizeof(interp)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Frees resources held by a c_polynomial_interp object.
    !!
    !! @param[in,out] obj The c_polynomial_interp object.
    subroutine polyinterp_free_c(obj) bind(C, name = "free_polynomial_interp")
        ! Arguments
        type(c_polynomial_interp), intent(inout), target :: obj

        ! Local Variables
        type(polynomial_interp), pointer :: interp

        ! Process
        call get_polynomial_interp(obj, interp)
        if (associated(interp)) deallocate(interp)
        obj%n = 0
        obj%ptr = c_null_ptr
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Performs a polynomial interpolation to determine the points @p y
    !! that for the requested indendent variable values in @p x.
    !!
    !! @param[in] obj The c_polynomial_interp object.
    !! @param[in] n The number of points to interpolate.
    !! @param[in] x An N-element array containing the values of the independent
    !!  variable at which to interpolate.
    !! @param[out] y An N-element array where the interpolated values can be
    !!  written.
    subroutine polyinterp_interp_c(obj, n, x, y) &
            bind(C, name = "polynomial_interpolate")
        ! Arguments
        type(c_polynomial_interp), intent(in), target :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64), intent(out) :: y(n)

        ! Local Variables
        type(polynomial_interp), pointer :: interp

        ! Process
        call get_polynomial_interp(obj, interp)
        if (.not.associated(interp)) return
        y = interp%interpolate(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the number of points used by the interpolation object.
    !!
    !! @param[in] obj The c_polynomial_interp object.
    !!
    !! @return The number of points.
    function polyinterp_get_pt_count_c(obj) result(n) &
            bind(C, name = "polynomial_interp_get_point_count")
        ! Arguments
        type(c_polynomial_interp), intent(in), target :: obj
        integer(int32) :: n

        ! Local Variables
        type(polynomial_interp), pointer :: interp

        ! Process
        n = 0
        call get_polynomial_interp(obj, interp)
        if (.not.associated(interp)) return
        n = interp%get_count()
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets a copy of the data points stored by the interpolation object.
    !!
    !! @param[in] obj The c_polynomial_interp object.
    !! @param[in] n The size of the buffer arrays.
    !! @param[out] x An N-element array where the x-coordinate data will be
    !!  written.
    !! @param[out] y An N-element array where the y-coordinate data will be
    !!  written.
    !!
    !! @par Remarks
    !! If @p n is different than the actual number of points that exist, the
    !! lesser of the two values will be utilized.  The interpolation object
    !! can be queried to determine the quantity of stored points.
    subroutine polyinterp_get_pts_c(obj, n, x, y) &
            bind(C, name = "polynomial_interp_get_points")
        ! Arguments
        type(c_polynomial_interp), intent(in), target :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(out) :: x(n), y(n)

        ! Local Variables
        integer(int32) :: i, npts
        type(polynomial_interp), pointer :: interp

        ! Process
        call get_polynomial_interp(obj, interp)
        if (.not.associated(interp)) return
        npts = min(n, interp%get_count())
        do concurrent (i = 1:npts)
            x(i) = interp%get_x(i)
            y(i) = interp%get_y(i)
        end do
    end subroutine

! ******************************************************************************
! SPLINE INTERPOLATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Retrieves the spline_interp object from the C compatible
    !! c_spline_interp data structure.
    !!
    !! @param[in] obj The C compatible c_spline_interp object.
    !! @param[out] interp A pointer to the resulting spline_interp object.
    !!  This pointer can be NULL dependent upon the state of @p obj.
    subroutine get_spline_interp(obj, interp)
        ! Arguments
        type(c_spline_interp), intent(in), target :: obj
        type(spline_interp), intent(out), pointer :: interp

        ! Local Variables
        type(c_ptr) :: testptr

        ! Process
        testptr = c_loc(obj)
        nullify(interp)
        if (.not.c_associated(testptr)) return
        if (.not.c_associated(obj%ptr)) return
        if (obj%n == 0) return
        call c_f_pointer(obj%ptr, interp)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Initializes a new c_spline_interp object.
    !!
    !! @param[out] obj The c_spline_interp object to initialize.
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the x-components of each data
    !!  point.  This array must be monotonic (ascending or descending only).
    !! @param[in] y An N-element array containing the y-components of each data
    !!  point.
    !! @param[in] ibcbeg An input that defines the nature of the
    !!  boundary condition at the beginning of the spline.  If an invalid
    !!  parameter is used, the code defaults to SPLINE_QUADRATIC_OVER_INTERVAL.
    !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
    !!      initial interval.  No value is required for @p ybcbeg.
    !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
    !!      initial point is provided in @p ybcbeg.
    !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
    !!      initial point is provided in @p ybcbeg.
    !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
    !!      at x(2).  No value is required for @p ybcbeg.
    !! @param[in] ybcbeg If needed, the value of the initial point boundary
    !!  condition.  If not needed, this parameter is ignored.
    !! @param[in] ibcend An input that defines the nature of the
    !!  boundary condition at the end of the spline.  If an invalid
    !!  parameter is used, the code defaults to SPLINE_QUADRATIC_OVER_INTERVAL.
    !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
    !!      final interval.  No value is required for @p ybcend.
    !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
    !!      initial point is provided in @p ybcend.
    !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
    !!      initial point is provided in @p ybcend.
    !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
    !!      at x(n-1).  No value is required for @p ybcend.
    !! @param[in] ybcend If needed, the value of the final point boundary
    !!  condition.  If not needed, this parameter is ignored.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically
    !!      increasing or decreasing.
    subroutine splineinterp_init_c(obj, n, x, y, ibcbeg, ybcbeg, ibcend, &
            ybcend, err) bind(C, name = "alloc_spline_interp")
        ! Arguments
        type(c_spline_interp), intent(out) :: obj
        integer(int32), intent(in), value :: n, ibcbeg, ibcend
        real(real64), intent(in), value :: ybcbeg, ybcend
        real(real64), intent(in) :: x(n), y(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr
        type(spline_interp), pointer :: interp

        ! Process
        allocate(interp)
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            call interp%initialize_spline(x, y, ibcbeg, ybcbeg, ibcend, &
                ybcend, eptr)
        else
            call interp%initialize_spline(x, y, ibcbeg, ybcbeg, ibcend, &
                ybcend)
        end if
        obj%ptr = c_loc(interp)
        obj%n = sizeof(interp)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Frees resources held by a c_spline_interp object.
    !!
    !! @param[in,out] obj The c_spline_interp object.
    subroutine splineinterp_free_c(obj) bind(C, name = "free_spline_interp")
        ! Arguments
        type(c_spline_interp), intent(inout), target :: obj

        ! Local Variables
        type(spline_interp), pointer :: interp

        ! Process
        call get_spline_interp(obj, interp)
        if (associated(interp)) deallocate(interp)
        obj%n = 0
        obj%ptr = c_null_ptr
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Performs a spline interpolation to determine the points @p y
    !! that for the requested indendent variable values in @p x.
    !!
    !! @param[in] obj The c_spline_interp object.
    !! @param[in] n The number of points to interpolate.
    !! @param[in] x An N-element array containing the values of the independent
    !!  variable at which to interpolate.
    !! @param[out] y An N-element array where the interpolated values can be
    !!  written.
    subroutine splineinterp_interp_c(obj, n, x, y) &
            bind(C, name = "spline_interpolate")
        ! Arguments
        type(c_spline_interp), intent(in), target :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64), intent(out) :: y(n)

        ! Local Variables
        type(spline_interp), pointer :: interp

        ! Process
        call get_spline_interp(obj, interp)
        if (.not.associated(interp)) return
        y = interp%interpolate(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the number of points used by the interpolation object.
    !!
    !! @param[in] obj The c_spline_interp object.
    !!
    !! @return The number of points.
    function splineinterp_get_pt_count_c(obj) result(n) &
            bind(C, name = "spline_interp_get_point_count")
        ! Arguments
        type(c_spline_interp), intent(in), target :: obj
        integer(int32) :: n

        ! Local Variables
        type(spline_interp), pointer :: interp

        ! Process
        n = 0
        call get_spline_interp(obj, interp)
        if (.not.associated(interp)) return
        n = interp%get_count()
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets a copy of the data points stored by the interpolation object.
    !!
    !! @param[in] obj The c_spline_interp object.
    !! @param[in] n The size of the buffer arrays.
    !! @param[out] x An N-element array where the x-coordinate data will be
    !!  written.
    !! @param[out] y An N-element array where the y-coordinate data will be
    !!  written.
    !!
    !! @par Remarks
    !! If @p n is different than the actual number of points that exist, the
    !! lesser of the two values will be utilized.  The interpolation object
    !! can be queried to determine the quantity of stored points.
    subroutine splineinterp_get_pts_c(obj, n, x, y) &
            bind(C, name = "spline_interp_get_points")
        ! Arguments
        type(c_spline_interp), intent(in), target :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(out) :: x(n), y(n)

        ! Local Variables
        integer(int32) :: i, npts
        type(spline_interp), pointer :: interp

        ! Process
        call get_spline_interp(obj, interp)
        if (.not.associated(interp)) return
        npts = min(n, interp%get_count())
        do concurrent (i = 1:npts)
            x(i) = interp%get_x(i)
            y(i) = interp%get_y(i)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the interpolated first derivative.
    !!
    !! @param[in] obj The c_spline_interp object.
    !! @param[in] n The number of points to interpolate.
    !! @param[in] x An N-element array containing the values of the independent
    !!  variable at which to interpolate.
    !! @param[out] y An N-element array where the interpolated values can be
    !!  written.
    subroutine splineinterp_diff1_c(obj, n, x, y) &
            bind(C, name = "spline_interp_diff1")
        ! Arguments
        type(c_spline_interp), intent(in), target :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64), intent(out) :: y(n)

        ! Local Variables
        type(spline_interp), pointer :: interp

        ! Process
        call get_spline_interp(obj, interp)
        if (.not.associated(interp)) return
        y = interp%first_derivative(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the interpolated second derivative.
    !!
    !! @param[in] obj The c_spline_interp object.
    !! @param[in] n The number of points to interpolate.
    !! @param[in] x An N-element array containing the values of the independent
    !!  variable at which to interpolate.
    !! @param[out] y An N-element array where the interpolated values can be
    !!  written.
    subroutine splineinterp_diff2_c(obj, n, x, y) &
            bind(C, name = "spline_interp_diff2")
        ! Arguments
        type(c_spline_interp), intent(in), target :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64), intent(out) :: y(n)

        ! Local Variables
        type(spline_interp), pointer :: interp

        ! Process
        call get_spline_interp(obj, interp)
        if (.not.associated(interp)) return
        y = interp%second_derivative(x)
    end subroutine

! ******************************************************************************
! STATISTICS ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the mean of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !!
    !! @return The mean of @p x.
    pure function mean_c(n, x) result(z) bind(C, name = "mean")
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64) :: z
        z = mean(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the median of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in,out] x The data set whose median is to be found.  Ideally, the
    !!  data set should be monotonically increasing; however, if it is not, it
    !!  may be sorted by the routine, dependent upon the value of @p srt.  On
    !!  output, the array contents are unchanged; however, they may be sorted
    !!  into ascending order (dependent upon the value of @p srt).
    !! @param[in] srt A logical flag determining if @p x should be sorted.
    !!
    !! @return The median of @p x.
    function median_c(n, x, srt) result(z) bind(C, name = "median")
        integer(int32), intent(in), value :: n
        real(real64), intent(inout) :: x(n)
        logical(c_bool), intent(in), value :: srt
        real(real64) :: z
        z = median(x, logical(srt))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the sample variance of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !!
    !! @return The variance of @p x.
    !!
    !! @par Remarks
    !! To avoid overflow-type issues, Welford's algorithm is employed.  A simple
    !! illustration of this algorithm can be found
    !! [here](https://www.johndcook.com/blog/standard_deviation/).
    pure function variance_c(n, x) result(v) bind(C, name = "variance")
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64) :: v
        v = variance(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the covariance matrix of N data sets of M observations.
    !!
    !! @param[in] m The number of observations.
    !! @param[in] n The number of data sets.
    !! @param[in] x The M-by-N matrix of data.
    !! @param[out] c The N-by-N matrix where the resulting covariance matrix
    !!  will be written.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    subroutine covariance_c(m, n, x, c, err) bind(C, name = "covariance")
        ! Arguments
        integer(int32), intent(in), value :: m, n
        real(real64), intent(in) :: x(m,n)
        real(real64), intent(out) :: c(n,n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            c = covariance(x, eptr)
        else
            c = covariance(x)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the corrected standard deviation of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !!
    !! @return The standard deviation of @p x.
    pure function stdev_c(n, x) result(s) bind(C, name = "standard_deviation")
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64) :: s
        s = standard_deviation(x)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the confidence interval based upon a standard normal
    !! distribution.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the data set.
    !! @param[in] alpha The confidence level.  This value must lie between
    !! zero and one such that: 0 < alpha < 1.
    !! @param[in] use_t Set to true to use the t-distribution in the event of
    !!  an unknown true standard deviation; else, set to true to use a normal
    !!  distribution.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p alpha is does not satisfy:
    !!      0 < alpha < 1.
    !!
    !! @return The confidence interval as the deviation from the mean.
    !!
    !! @par Remarks
    !! The confidence interval, assuming a standard normal distribution, is
    !! as follows: mu +/- z * s / sqrt(n), where mu = the mean, and s = the
    !! standard deviation.  This routine computes the z * s / sqrt(n) portion
    !! leaving the computation of the mean to the user.
    function conf_int_c(n, x, alpha, use_t, err) result(c) &
            bind(C, name = "confidence_interval")
        ! Arguments
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64), intent(in), value :: alpha
        logical(c_bool), intent(in), value :: use_t
        type(errorhandler), intent(inout) :: err
        real(real64) :: c

        ! Local Variables
        type(errors), pointer :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            c = confidence_interval(x, alpha, logical(use_t), eptr)
        else
            c = confidence_interval(x, alpha, logical(use_t))
        end if
    end function

! ******************************************************************************
! REGRESSION & SMOOTHING ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Applies a moving average to smooth a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in,out] x On input, the signal to smooth.  On output, the smoothed
    !!  signal.
    !! @param[in] npts The size of the averaging window.  This value must be
    !!  at least 2, but no more than the number of elements in @p x.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p npts is less than 2, or greater
    !!      than the length of @p x.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    subroutine moving_average_c(n, x, npts, err) &
            bind(C, name = "moving_average")
        ! Arguments
        integer(int32), intent(in), value :: n, npts
        real(real64), intent(inout) :: x(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            call moving_average(x, npts, eptr)
        else
            call moving_average(x, npts)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Employs a least squares fit to determine the coefficient A in the
    !! linear system: Y = A * X.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the independent variable data.
    !! @param[in,out] y An N-element array containing the dependent variable
    !!  data corresponding to @p x.  On output, the contents of this array are
    !!  overwritten as it is used for storage purposes by the algorithm.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are different sizes.
    !!
    !! @return The scalar coefficient A.
    function linlsq_1var_c(n, x, y, err) result(a) &
            bind(C, name = "least_squares_fit_1var")
        ! Arguments
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n)
        real(real64), intent(inout) :: y(n)
        type(errorhandler), intent(inout) :: err
        real(real64) :: a

        ! Local Variables
        type(errors), pointer :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            a = linear_least_squares(x, y, eptr)
        else
            a = linear_least_squares(x, y)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Employs a least squares fit to determine the coefficient A in the
    !! linear system: Y = A * X.
    !!
    !! @param[in] m The number of dependent variables.
    !! @param[in] n The number of independent variables.
    !! @param[in,out] x An N-by-NPTS matrix containing the P data points of the
    !!  N independent variables.
    !! @param[in] y An M-by-NPTS matrix containing the P data points of the M
    !!  dependent variables.
    !! @param[out] a The M-by-N matrix where the resulting coefficient matrix A
    !!  will be written.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if any of the matrix dimensions are not
    !!      compatiable.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!
    !! @par Remarks
    !! The algorithm attempts to compute the coefficient matrix A as follows.
    !! Y * X**T = A * X * X**T
    !! Y * X**T * INV(X * X**T) = A
    !! This does require that X * X**T does not result in a singular matrix.  To
    !! handle the situation where X * X**T is singular, the Moore-Penrose
    !! pseudo-inverse, computed by means of singular value decomposition, is
    !! utilized to still arrive at a solution that, at minimum, has a minimum
    !! Euclidean norm of its residual.
    !! Let: PINV(X) = X**T * INV(X * X**T),
    !! Then: A = Y * PINV(X)
    subroutine linlsq_nvar_c(m, n, npts, x, y, a, err) &
            bind(C, name = "least_squares_fit_nvar")
        ! Arguments
        integer(int32), intent(in), value :: m, n, npts
        real(real64), intent(inout) :: x(n,npts)
        real(real64), intent(in) :: y(m,npts)
        real(real64), intent(out) :: a(m,n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            a = linear_least_squares(x, y, err = eptr)
        else
            a = linear_least_squares(x, y)
        end if
    end subroutine

! ******************************************************************************
! LOWESS SMOOTHING ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Retrieves the lowess_smoothing object from the C compatible
    !! c_lowess_smoothing data structure.
    !!
    !! @param[in] obj The C compatible c_lowess_smoothing object.
    !! @param[out] ptr A pointer to the resulting lowess_smoothing object.  This
    !!  pointer can be NULL dependent upon the state of @p obj.
    subroutine get_lowess_smoothing(obj, ptr)
        ! Arguments
        type(c_lowess_smoothing), intent(in), target :: obj
        type(lowess_smoothing), intent(out), pointer :: ptr

        ! Local Variables
        type(c_ptr) :: testptr

        ! Process
        testptr = c_loc(obj)
        nullify(ptr)
        if (.not.c_associated(testptr)) return
        if (.not.c_associated(obj%ptr)) return
        if (obj%n == 0) return
        call c_f_pointer(obj%ptr, ptr)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Initializes a new c_lowess_smoothing object.
    !!
    !! @param[out] obj The c_lowess_smoothing object.
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element array containing the x-coordinate data.
    !!  Ideally, the data set should be monotonically increasing; however, if
    !!  it is not, it may be sorted by the routine, dependent upon the value
    !!  of @p srt.
    !! @param[in] y An N-element array containing the y-coordinate data.
    !! @param[in] srt A logical flag determining if @p x should be sorted.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    subroutine lowess_init_c(obj, n, x, y, srt, err) &
            bind(C, name = "alloc_lowess")
        ! Arguments
        type(c_lowess_smoothing), intent(out) :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: x(n), y(n)
        logical(c_bool), intent(in), value :: srt
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr
        type(lowess_smoothing), pointer :: ptr

        ! Process
        allocate(ptr)
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            call ptr%initialize(x, y, logical(srt), eptr)
        else
            call ptr%initialize(x, y, logical(srt))
        end if
        obj%ptr = c_loc(ptr)
        obj%n = sizeof(ptr)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Frees resources held by a c_lowess_smoothing object.
    !!
    !! @param[in,out] obj The c_lowess_smoothing object.
    subroutine lowess_free_c(obj) bind(C, name = "free_lowess")
        ! Arguments
        type(c_lowess_smoothing), intent(inout), target :: obj

        ! Local Variables
        type(lowess_smoothing), pointer :: ptr

        ! Process
        call get_lowess_smoothing(obj, ptr)
        if (associated(ptr)) deallocate(ptr)
        obj%n = 0
        obj%ptr = c_null_ptr
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Performs the actual smoothing operation.
    !!
    !! @param[in,out] obj The c_lowess_smoothing object.
    !! @param[in] f Specifies the amount of smoothing.  More specifically, this
    !! value is the fraction of points used to compute each value.  As this
    !! value increases, the output becomes smoother.  Choosing a value in the
    !! range of 0.2 to 0.8 usually results in a good fit.  As such, a reasonable
    !! starting point, in the absence of better information, is a value of 0.5.
    !! @param[in] n The size of the buffer @p y.  Ideally, this parameter is
    !!  equal to the number of points stored in @p obj; however, the routine
    !!  will only traverse the minimum of the this parameter or the number of
    !!  points stored in @p obj.
    !! @param[out] y An N-element array to which the smoothed data will be
    !!  written.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_NO_DATA_DEFINED_ERROR: Occurs if no data has been defined.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    subroutine lowess_smooth_c(obj, f, n, y, err) &
            bind(C, name = "lowess_smooth")
        ! Arguments
        type(c_lowess_smoothing), intent(inout) :: obj
        real(real64), intent(in), value :: f
        integer(int32), intent(in), value :: n
        real(real64), intent(out) :: y(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr
        type(lowess_smoothing), pointer :: ptr
        integer(int32) :: npts

        ! Process
        call get_errorhandler(err, eptr)
        call get_lowess_smoothing(obj, ptr)
        if (.not.associated(ptr)) return
        npts = min(n, ptr%get_count())
        if (associated(eptr)) then
            y(1:npts) = ptr%smooth(f, eptr)
        else
            y(1:npts) = ptr%smooth(f)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the number of points used by the lowess_smoothing object.
    !!
    !! @param[in] obj The c_lowess_smoothing object.
    !!
    !! @return The number of points.
    function lowess_get_pt_count_c(obj) result(n) &
            bind(C, name = "lowess_get_point_count")
        ! Arguments
        type(c_lowess_smoothing), intent(in) :: obj
        integer(int32) :: n

        ! Local Variables
        type(lowess_smoothing), pointer :: ptr

        ! Process
        n = 0
        call get_lowess_smoothing(obj, ptr)
        if (.not.associated(ptr)) return
        n = ptr%get_count()
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets a copy of the data points stored by the lowess_smoothing
    !! object.
    !!
    !! @param[in] obj The c_lowess_smoothing object.
    !! @param[in] n The size of the buffer arrays.
    !! @param[out] x An N-element array where the x-coordinate data will be
    !!  written.
    !! @param[out] y An N-element array where the y-coordinate data will be
    !!  written.
    !!
    !! @par Remarks
    !! If @p n is different than the actual number of points that exist, the
    !! lesser of the two values will be utilized.  The lowess_smoothing object
    !! can be queried to determine the quantity of stored points.
    subroutine lowess_get_pts_c(obj, n, x, y) &
            bind(C, name = "lowess_get_points")
        ! Arguments
        type(c_lowess_smoothing), intent(in) :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(out) :: x(n), y(n)

        ! Local Variables
        integer(int32) :: i, npts
        type(lowess_smoothing), pointer :: ptr

        ! Process
        call get_lowess_smoothing(obj, ptr)
        if (.not.associated(ptr)) return
        npts = min(n, ptr%get_count())
        do concurrent (i = 1:npts)
            x(i) = ptr%get_x(i)
            y(i) = ptr%get_y(i)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the residuals from each data point.
    !!
    !! @param[in] this The c_lowess_smoothing object.
    !! @param[in] n The number of elements available in the buffer array @p x.
    !! @param[out] x An N-element array where the residual data should be
    !!  written.
    subroutine lowess_get_residual_c(obj, n, x) &
            bind(C, name = "lowess_get_residuals")
        ! Arguments
        type(c_lowess_smoothing), intent(in) :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(out) :: x(n)

        ! Local Variables
        type(lowess_smoothing), pointer :: ptr

        ! Process
        call get_lowess_smoothing(obj, ptr)
        if (.not.associated(ptr)) return
        call ptr%get_residuals(x)
    end subroutine

! ******************************************************************************
! NONLINEAR REGRESSION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Retrieves the cnonlin_reg_helper object from the C compatible
    !! c_nonlinear_regression data structure.
    !!
    !! @param[in] obj The C compatible c_nonlinear_regression object.
    !! @param[out] ptr A pointer to the resulting cnonlin_reg_helper object.
    !!  This pointer can be NULL dependent upon the state of @p obj.
    subroutine get_nonlinear_regression(obj, ptr)
        ! Arguments
        type(c_nonlinear_regression), intent(in), target :: obj
        type(cnonlin_reg_helper), intent(out), pointer :: ptr

        ! Local Variables
        type(c_ptr) :: testptr

        ! Process
        testptr = c_loc(obj)
        nullify(ptr)
        if (.not.c_associated(testptr)) return
        if (.not.c_associated(obj%ptr)) return
        if (obj%n == 0) return
        call c_f_pointer(obj%ptr, ptr)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Initializes a new c_nonlinear_regression object.
    !!
    !! @param[out] obj The c_nonlinear_regression object.
    !! @param[in] n The number of data points.
    !! @param[in] x An N-element containing the independent variable values of
    !!  the data set.
    !! @param[in] y  An N-element array of the dependent variables corresponding
    !!  to @p x.
    !! @param[in] fcn A pointer to the function whose coefficients are to be
    !!  determined.
    !! @param[in] ncoeff The number of coefficients in the function defined in
    !!  @p fcn.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p ncoeff is less than or equal to
    !!      zero.
    subroutine nlr_init_c(obj, n, x, y, fcn, ncoeff, err) &
            bind(C, name = "alloc_nonlinear_regression")
        ! Arguments
        type(c_nonlinear_regression), intent(out) :: obj
        integer(int32), intent(in), value :: n, ncoeff
        real(real64), intent(in) :: x(n), y(n)
        type(c_funptr), intent(in), value :: fcn
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr
        type(cnonlin_reg_helper), pointer :: ptr
        procedure(creg_fcn), pointer :: fptr
        procedure(reg_fcn), pointer :: nullfptr

        ! Process
        nullfptr => null()
        call c_f_procpointer(fcn, fptr)
        allocate(ptr)
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            call ptr%initialize(x, y, nullfptr, ncoeff, eptr)
        else
            call ptr%initialize(x, y, nullfptr, ncoeff, eptr)
        end if
        call ptr%set_cfcn(fptr)
        obj%ptr = c_loc(ptr)
        obj%n = sizeof(ptr)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Frees resources held by a c_nonlinear_regression object.
    !!
    !! @param[in,out] obj The c_nonlinear_regression object.
    subroutine nlr_free_c(obj) bind(C, name = "free_nonlinear_regression")
        ! Arguments
        type(c_nonlinear_regression), intent(inout), target :: obj

        ! Local Variables
        type(cnonlin_reg_helper), pointer :: ptr

        ! Process
        call get_nonlinear_regression(obj, ptr)
        if (associated(ptr)) deallocate(ptr)
        obj%n = 0
        obj%ptr = c_null_ptr
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the solution to the nonlinear regression problem using
    !! the Levenberg-Marquardt method.
    !!
    !! @param[in,out] obj The c_nonlinear_regression object.
    !! @param[in] n The number of coefficients to determine.
    !! @param[in,out] c On input, an array containing initial estimates of the
    !!  coefficients.  On output, the comptued coefficient values.
    !! @param[out] ib An output parameter that allows the caller to obtain
    !!  iteration performance statistics.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if the number of equations is less than
    !!      than the number of variables.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      correctly.
    !!  - CF_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
    !!      the allowed number of iterations.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - CF_TOLERANCE_TOO_SMALL_ERROR: Occurs if the requested tolerance is
    !!      to small to be practical for the problem at hand.
    subroutine nlr_solve_c(obj, n, c, ib, err) &
            bind(C, name = "nonlinreg_solve")
        ! Arguments
        type(c_nonlinear_regression), intent(inout) :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(inout) :: c(n)
        type(iteration_behavior), intent(out) :: ib
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr
        type(cnonlin_reg_helper), pointer :: cptr

        ! Process
        call get_nonlinear_regression(obj, cptr)
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            call cptr%solve(c, ib = ib, err = eptr)
        else
            call cptr%solve(c, ib = ib)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the number of points used by the c_nonlinear_regression
    !! object.
    !!
    !! @param[in] obj The c_nonlinear_regression object.
    !!
    !! @return The number of points.
    function nlr_get_pt_count_c(obj) result(n) &
            bind(C, name = "nonlinreg_get_point_count")
        ! Arguments
        type(c_nonlinear_regression), intent(in) :: obj
        integer(int32) :: n

        ! Local Variables
        type(cnonlin_reg_helper), pointer :: cptr

        ! Process
        n = 0
        call get_nonlinear_regression(obj, cptr)
        if (.not.associated(cptr)) return
        n = cptr%get_count()
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets a copy of the data points stored by the
    !! c_nonlinear_regression object.
    !!
    !! @param[in] obj The c_nonlinear_regression object.
    !! @param[in] n The size of the buffer arrays.
    !! @param[out] x An N-element array where the x-coordinate data will be
    !!  written.
    !! @param[out] y An N-element array where the y-coordinate data will be
    !!  written.
    !!
    !! @par Remarks
    !! If @p n is different than the actual number of points that exist, the
    !! lesser of the two values will be utilized.  The c_nonlinear_regression
    !! object can be queried to determine the quantity of stored points.
    subroutine nlr_get_pts_c(obj, n, x, y) &
            bind(C, name = "nonlinreg_get_points")
        ! Arguments
        type(c_nonlinear_regression), intent(in) :: obj
        integer(int32), intent(in), value :: n
        real(real64), intent(out) :: x(n), y(n)

        ! Local Variables
        integer(int32) :: i, npts
        type(cnonlin_reg_helper), pointer :: ptr

        ! Process
        call get_nonlinear_regression(obj, ptr)
        if (.not.associated(ptr)) return
        npts = min(n, ptr%get_count())
        do concurrent (i = 1:npts)
            x(i) = ptr%get_x(i)
            y(i) = ptr%get_y(i)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the nonlinear regression solver solution control parameters.
    !!
    !! @param[in] obj The c_nonlinear_regression object.
    !! @param[out] cntrl The solver_control object that, on output, will contain
    !!  the current solver control parameters.
    subroutine nlr_get_solver_params_c(obj, cntrl) &
            bind(C, name = "nonlinreg_get_solver_params")
        ! Arguments
        type(c_nonlinear_regression), intent(in) :: obj
        type(solver_control), intent(out) :: cntrl

        ! Local Variables
        type(cnonlin_reg_helper), pointer :: ptr

        ! Process
        call get_nonlinear_regression(obj, ptr)
        if (.not.associated(ptr)) return
        cntrl%max_evals = ptr%get_max_fcn_evals()
        cntrl%fcn_tolerance = ptr%get_fcn_tolerance()
        cntrl%var_tolerance = ptr%get_var_tolerance()
        cntrl%grad_tolerance = ptr%get_gradient_tolerance()
        cntrl%print_status = ptr%get_print_status()
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Sets  the nonlinear regression solver solution control parameters.
    !!
    !! @param[in,out] obj The c_nonlinear_regression object.
    !! @param[in] cntrl The solver_control object that contains the current
    !!  solver control parameters.
    subroutine nlr_set_solver_params_c(obj, cntrl) &
            bind(C, name = "nonlinreg_set_solver_params")
        ! Arguments
        type(c_nonlinear_regression), intent(inout) :: obj
        type(solver_control), intent(in) :: cntrl

        ! Local Variables
        type(cnonlin_reg_helper), pointer :: ptr

        ! Process
        call get_nonlinear_regression(obj, ptr)
        if (.not.associated(ptr)) return
        call ptr%set_max_fcn_evals(cntrl%max_evals)
        call ptr%set_fcn_tolerance(cntrl%fcn_tolerance)
        call ptr%set_var_tolerance(cntrl%var_tolerance)
        call ptr%set_gradient_tolerance(cntrl%grad_tolerance)
        call ptr%set_print_status(logical(cntrl%print_status))
    end subroutine

! ******************************************************************************
! CNONLIN_REG_HELPER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Computes the residual between the supplied data set, and the
    !! function value given a set of coefficients.
    !!
    !! @param[in] this The cnonlin_reg_helper object.
    !! @param[in] x An N-element array containing the N coefficients.
    !! @param[out] f An M-element array that, on output, contains the residual
    !!  at each of the M data points.
    subroutine crh_fcn(this, x, f)
        ! Arguments
        class(cnonlin_reg_helper), intent(in) :: this
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f

        ! Local Variables
        integer(int32) :: i, n, ncoeff

        ! Compute the value of the function at each value of x (get_x)
        n = this%get_equation_count()
        ncoeff = this%get_variable_count()
        do i = 1, n
            f(i) = this%get_y(i) - this%m_cfcn(this%get_x(i), ncoeff, x)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the function containing the equation to
    !! solve has been assigned.
    !!
    !! @param[in] this The cnonlin_reg_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function crh_is_fcn_defined(this) result(x)
        class(cnonlin_reg_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_cfcn)
    end function

! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine containing the equations to
    !! solve.
    !!
    !! @param[in,out] this The cnonlin_reg_helper object.
    !! @param[in] fcn The function pointer.
    subroutine crh_set_fcn(this, fcn)
        class(cnonlin_reg_helper), intent(inout) :: this
        procedure(creg_fcn), intent(in), pointer :: fcn
        this%m_cfcn => fcn
    end subroutine

! ******************************************************************************
! CALIBRATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Computes the static error band of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] output An N-element array containing the values output by
    !!  the instrument as a result of the values given in @p applied.
    !! @param[in] fullscale The full scale measurement value for the instrument.
    !!  The units must be consistent with those of @p applied.
    !! @param[out] rst An seb_results object where the calculation results will
    !!  be written.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p fullscale is sufficiently close
    !!      to zero to be considered zero.  Sufficiently close in this instance
    !!      is considered to be the square root of machine precision.
    subroutine seb_c(n, applied, output, fullscale, rst, err) &
            bind(C, name = "seb")
        ! Arguments
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: applied(n), output(n)
        real(real64), intent(in), value :: fullscale
        type(seb_results), intent(out) :: rst
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            rst = seb(applied, output, fullscale, eptr)
        else
            rst = seb(applied, output, fullscale)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the best-fit nonlinearity of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] measured An N-element array containing the calibrated output
    !!  of the instrument as a result of the values given in @p applied.
    !!
    !! @return The nonlinearity error.
    function nonlin_c(n, applied, measured) result(rst) &
            bind(C, name = "nonlinearity")
        ! Arguments
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: applied(n), measured(n)
        real(real64) :: rst

        ! Process
        rst = nonlinearity(applied, measured)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the terminal nonlinearity of a data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] measured An N-element array containing the calibrated output
    !!  of the instrument as a result of the values given in @p applied.
    !!
    !! @return The nonlinearity error.
    function term_nonlin_c(n, applied, measured) result(rst) &
            bind(C, name = "terminal_nonlinearity")
        ! Argument
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: applied(n), measured(n)
        real(real64) :: rst

        ! Process
        rst = terminal_nonlinearity(applied, measured)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the hysteresis in an ascending/descending data set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] measured An N-element array containing the calibrated output
    !!  of the instrument as a result of the values given in @p applied.
    !!
    !! @return The hysteresis error.
    function hysteresis_c(n, applied, measured) result(rst) &
            bind(C, name = "hysteresis")
        ! Arguments
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: applied(n), measured(n)
        real(real64) :: rst

        ! Process
        rst = hysteresis(applied, measured)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the return to zero error in an ascending/descending data
    !! set.
    !!
    !! @param[in] n The number of data points.
    !! @param[in] applied An N-element array containing the values applied to
    !!  the measurement instrument.
    !! @param[in] measured An N-element array containing the calibrated output
    !!  of the instrument as a result of the values given in @p applied.
    !! @param[in] tol An input argument that specifies the tolerance used in
    !!  finding the matching zero data point.
    !!
    !! @return The return to zero error.
    function rtz_c(n, applied, measured, tol) result(rst) &
            bind(C, name = "return_to_zero")
        ! Arguments
        integer(int32), intent(in), value :: n
        real(real64), intent(in) :: applied(n), measured(n)
        real(real64), intent(in), value :: tol
        real(real64) :: rst

        ! Process
        rst = return_to_zero(applied, measured, tol)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the repeatability of a sequence of tests.
    !!
    !! @param[in] npts The number of data points per test.
    !! @param[in] ntests The number of tests.
    !! @param[in] applied An NPTS-by-NTEST matrix containing at least 2 columns
    !!  (tests) of NPTS values applied to the measurement instrument.
    !! @param[in] measured An NPTS-by-NTEST matrix containing the corresponding
    !!  calibrated output from the instrument.
    !!
    !! @return The largest magnitude deviation from the initial test.
    !!
    !! @par Remarks
    !! Repeatability is considered as the largest magnitude deviation of
    !! subsequent tests from the initial test.  Noting that it is very likely
    !! that consecutive test points will vary slightly, test 2 through test N
    !! are linearly interpolated such that their test points line up with those
    !! from test 1.
    function repeat_c(npts, ntests, applied, measured) result(rst) &
            bind(C, name = "repeatability")
        ! Arguments
        integer(int32), intent(in), value :: npts, ntests
        real(real64), intent(in) :: applied(npts, ntests), measured(npts, ntests)
        real(real64) :: rst

        ! Process
        rst = repeatability(applied, measured)
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the crosstalk errors for a multiple degree-of-freedom
    !! data set.
    !!
    !! @param[in] npts The number of data points in each degree of freedom.
    !! @param[in] ndof The number of degrees of freedom.
    !! @param[in] xerr An NPTS-by-NDOF matrix containing the measurement error
    !!  values (computed such that XERR = X MEASURED - X APPLIED).
    !! @param[in] indices A 2*NDOF element array containing row indices defining
    !!  the rows where each degree-of-freedom was applied in the data set
    !!  @p xerr.
    !! @param[out] xt An NDOF-by-NDOF matrix that, on output, will contain the
    !!  crosstalk errors such that each loaded degree of freedom is represented
    !!  by its own row, and each responding degree of freedom is represented by
    !!  its own column.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - CF_ARRAY_INDEX_ERROR: Occurs if any of the entries in @p indices are
    !!      outside the row bounds of @p xerr.
    subroutine xtalk_c(npts, ndof, xerr, indices, xt, err) &
            bind(C, name = "crosstalk")
        ! Arguments
        integer(int32), intent(in), value :: npts, ndof
        real(real64), intent(in) :: xerr(npts, ndof)
        integer(int32), intent(in) :: indices(2*ndof)
        real(real64), intent(out) :: xt(ndof, ndof)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            xt = crosstalk(xerr, indices, eptr)
        else
            xt = crosstalk(xerr, indices)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Splits a data set into ascending and descending components.
    !!
    !! @param[in] n The number of data points in @p x.
    !! @param[in] x An N-element array containing the data set to split.
    !! @param[in] na The capacity of @p ascend.
    !! @param[out] ascend An array where the ascending points will be written.
    !!  Ensure this array is appropriately sized to accept all the ascending
    !!  points (it can be oversized).
    !! @param[in] nd The capacity of @p descend.
    !! @param[out] descend An array where the descending points will be written.
    !!  Ensure this array is appropriately sized to accept all the descending
    !!  points (it can be oversized).
    !! @param[out] nascend The actual number of values written into @p ascend.
    !! @param[out] ndescend The actual number of values written into @p descend.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
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
    subroutine split_c(n, x, na, ascend, nd, descend, nascend, ndescend, err) &
            bind(C, name = "split_ascend_descend")
        ! Arguments
        integer(int32), intent(in), value :: n, na, nd
        real(real64), intent(in) :: x(n)
        real(real64), intent(out) :: ascend(na), descend(nd)
        integer(int32), intent(out) :: nascend, ndescend
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr

        ! Process
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            call split_ascend_descend(x, ascend, descend, nascend, ndescend, &
                eptr)
        else
            call split_ascend_descend(x, ascend, descend, nascend, ndescend)
        end if
    end subroutine

! ------------------------------------------------------------------------------
end module
