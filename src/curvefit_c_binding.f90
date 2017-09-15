! curvefit_c_binding.f90

!> @breif \b curvefit_c_binding
!!
!! @par Purpose
!! Provides C bindings to the curvefit library.
module curvefit_c_binding
    use, intrinsic :: iso_c_binding
    use curvefit_core
    use curvefit_interp
    use curvefit_statistics
    use curvefit_regression
    use ferror, only : errors
    use ferror_c_binding, only : errorhandler, get_errorhandler
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
            real(dp), intent(in), value :: x
            integer(i32), intent(in), value :: n
            real(dp), intent(in) :: c(n)
            real(dp) :: f
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
        integer(i32) :: n
    end type

! ------------------------------------------------------------------------------
    !> @brief A C compatible type encapsulating a polynomial_interp object.
    type, bind(C) :: c_polynomial_interp
        !> @brief A pointer to the polynomial_interp object.
        type(c_ptr) :: ptr
        !> @brief The size of the polynomial_interp object, in bytes.
        integer(i32) :: n
    end type

! ------------------------------------------------------------------------------
    !> @brief A C compatible type encapsulating a spline_interp object.
    type, bind(C) :: c_spline_interp
        !> @brief A pointer to the spline_interp object.
        type(c_ptr) :: ptr
        !> @brief The size of the spline_interp object, in bytes.
        integer(i32) :: n
    end type

! ------------------------------------------------------------------------------
    !> @brief A C compatible type encapsulating a lowess_smoothing object.
    type, bind(C) :: c_lowess_smoothing
        !> @brief A pointer to the lowess_smoothing object.
        type(c_ptr) :: ptr
        !> @brief The size of the lowess_smoothing object, in bytes.
        integer(i32) :: n
    end type

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n), y(n)
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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: y(n)

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
        integer(i32) :: n

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
        integer(i32), intent(in), value :: n
        real(dp), intent(out) :: x(n), y(n)

        ! Local Variables
        integer(i32) :: i, npts
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
        integer(i32), intent(in), value :: n, order
        real(dp), intent(in) :: x(n), y(n)
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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: y(n)

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
        integer(i32) :: n

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
        integer(i32), intent(in), value :: n
        real(dp), intent(out) :: x(n), y(n)

        ! Local Variables
        integer(i32) :: i, npts
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
        integer(i32), intent(in), value :: n, ibcbeg, ibcend
        real(dp), intent(in), value :: ybcbeg, ybcend
        real(dp), intent(in) :: x(n), y(n)
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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: y(n)

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
        integer(i32) :: n

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
        integer(i32), intent(in), value :: n
        real(dp), intent(out) :: x(n), y(n)

        ! Local Variables
        integer(i32) :: i, npts
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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: y(n)

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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: y(n)

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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp) :: z
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
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: x(n)
        logical(c_bool), intent(in), value :: srt
        real(dp) :: z
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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp) :: v
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
        integer(i32), intent(in), value :: m, n
        real(dp), intent(in) :: x(m,n)
        real(dp), intent(out) :: c(n,n)
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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp) :: s
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
    !!
    !! @return The confidence interval as the deviation from the mean.
    !!
    !! @par Remarks
    !! The confidence interval, assuming a standard normal distribution, is
    !! as follows: mu +/- z * s / sqrt(n), where mu = the mean, and s = the
    !! standard deviation.  This routine computes the z * s / sqrt(n) portion 
    !! leaving the computation of the mean to the user.
    function conf_int_c(n, x, alpha) result(c) &
            bind(C, name = "confidence_interval")
        integer(i32), intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(in), value :: alpha
        real(dp) :: c
        c = confidence_interval(x, alpha)
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
        integer(i32), intent(in), value :: n, npts
        real(dp), intent(inout) :: x(n)
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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(inout) :: y(n)
        type(errorhandler), intent(inout) :: err
        real(dp) :: a

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
        integer(i32), intent(in), value :: m, n, npts
        real(dp), intent(inout) :: x(n,npts)
        real(dp), intent(in) :: y(m,npts)
        real(dp), intent(out) :: a(m,n)
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
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n), y(n)
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
        real(dp), intent(in), value :: f
        integer(i32), intent(in), value :: n
        real(dp), intent(out) :: y(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr
        type(lowess_smoothing), pointer :: ptr
        integer(i32) :: npts

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
        integer(i32) :: n

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
        integer(i32), intent(in), value :: n
        real(dp), intent(out) :: x(n), y(n)

        ! Local Variables
        integer(i32) :: i, npts
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
        integer(i32), intent(in), value :: n
        real(dp), intent(out) :: x(n)

        ! Local Variables
        type(lowess_smoothing), pointer :: ptr 

        ! Process
        call get_lowess_smoothing(obj, ptr)
        if (.not.associated(ptr)) return
        call ptr%get_residuals(x)
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
