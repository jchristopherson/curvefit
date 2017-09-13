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
    !> @brief A C compatiable type encapsulating a linear_interp object.
    type, bind(C) :: c_linear_interp
        !> @brief A pointer to the linear_interp object.
        type(c_ptr) :: ptr
        !> @brief The size of the linear_interp object, in bytes.
        integer(i32) :: n
    end type


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
            x(i) = this%get_x(i)
            y(i) = this%get_y(i)
        end do
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
