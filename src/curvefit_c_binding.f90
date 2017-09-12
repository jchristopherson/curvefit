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
            use curvefit_core, only : dp
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
    !

! ------------------------------------------------------------------------------
    !> @brief Initializes a new c_linear_interp object.
    !!
    !! @param[out] obj The c_linear_interp object to initialize.
    !! @param[in] n
    !! @param[in] x
    !! @param[in] y
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    subroutine lininterp_init_c(obj, n, x, y, err) bind(C, name = "init_linear_interp")
        ! Arguments
        type(c_linear_interp), intent(out) :: obj
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n), y(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        type(errors), pointer :: eptr
        type(linear_interp), pointer :: optr

        ! Process
        call get_errorhandler(err, eptr)
        if (associated(eptr)) then
            call optr%initialize(x, y, 1, eptr)
        else
            call optr%initialize(x, y, 1)
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
