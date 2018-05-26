! curvefit_core.f90

!> @mainpage
!!
!! @section intro_sec Introduction
!! CURVEFIT is a library for fitting functions to sets of data.
!! @image html lowess_example_1.png

!> @brief \b curvefit_core
!!
!! @par Purpose
!! To provide core types and routines for the CURVEFIT library.
module curvefit_core
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use nonlin_constants
    implicit none
    private
    public :: CF_ARRAY_SIZE_ERROR
    public :: CF_OUT_OF_MEMORY_ERROR
    public :: CF_NO_DATA_DEFINED_ERROR
    public :: CF_INVALID_INPUT_ERROR
    public :: CF_NONMONOTONIC_ARRAY_ERROR
    public :: CF_INVALID_OPERATION_ERROR
    public :: CF_CONVERGENCE_ERROR
    public :: CF_TOLERANCE_TOO_SMALL_ERROR
    public :: CF_ARRAY_INDEX_ERROR
    public :: is_monotonic
    public :: reg_fcn

! ******************************************************************************
! ERROR FLAGS
! ------------------------------------------------------------------------------
    !> An error flag denoting an improperly sized array.
    integer(int32), parameter :: CF_ARRAY_SIZE_ERROR = NL_ARRAY_SIZE_ERROR
    !> An error denoting that there is insufficient memory available.
    integer(int32), parameter :: CF_OUT_OF_MEMORY_ERROR = NL_OUT_OF_MEMORY_ERROR
    !> An error denoting that no data has been defined.
    integer(int32), parameter :: CF_NO_DATA_DEFINED_ERROR = 300
    !> An error flag denoting an invalid input.
    integer(int32), parameter :: CF_INVALID_INPUT_ERROR = NL_INVALID_INPUT_ERROR
    !> An error flag denoting a non-monotonic array was given when a monotonic
    !! array was expected.
    integer(int32), parameter :: CF_NONMONOTONIC_ARRAY_ERROR = 301
    !> An error resulting from an invalid operation.
    integer(int32), parameter :: CF_INVALID_OPERATION_ERROR = &
        NL_INVALID_OPERATION_ERROR
    !> An error resulting from a lack of convergence.
    integer(int32), parameter :: CF_CONVERGENCE_ERROR = NL_CONVERGENCE_ERROR
    !> An error indicating the user-requested tolerance is too small to be
    !! practical for the problem at hand.
    integer(int32), parameter :: CF_TOLERANCE_TOO_SMALL_ERROR = &
        NL_TOLERANCE_TOO_SMALL_ERROR
    !> An error indicating an array index was out of bounds.
    integer(int32), parameter :: CF_ARRAY_INDEX_ERROR = 302

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Tests to see if an array is montonically increasing or decreasing.
    interface is_monotonic
        module procedure :: is_monotonic_dbl
        module procedure :: is_monotonic_i32
    end interface

! ------------------------------------------------------------------------------
    interface
        !> @brief Describes a routine for finding the coefficients of a function
        !! of one variable.
        !!
        !! @param[in] x The independent variable.
        !! @param[in] c An array of function coefficients.
        !!
        !! @result The value of the function at @p x.
        function reg_fcn(x, c) result(f)
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: c
            real(real64) :: f
        end function
    end interface


contains
! ------------------------------------------------------------------------------
    !> @brief Tests to see if an array is montonically increasing or decreasing.
    !!
    !! @param[in] x The array to test.
    !!
    !! @return Returns true if @p x is monotonic; else, false.
    pure function is_monotonic_dbl(x) result(rst)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x
        logical :: rst

        ! Process
        integer(int32) :: i, n
        logical :: ascend
        rst = .true.
        n = size(x)
        ascend = x(n) > x(1)
        if (ascend) then
            do i = 2, n
                if (x(i) <= x(i-1)) then
                    rst = .false.
                    exit
                end if
            end do
        else
            do i = 2, n
                if (x(i) >= x(i-1)) then
                    rst = .false.
                    exit
                end if
            end do
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Tests to see if an array is montonically increasing or decreasing.
    !!
    !! @param[in] x The array to test.
    !!
    !! @return Returns true if @p x is monotonic; else, false.
    pure function is_monotonic_i32(x) result(rst)
        ! Arguments
        integer(int32), intent(in), dimension(:) :: x
        logical :: rst

        ! Process
        integer(int32) :: i, n
        logical :: ascend
        rst = .true.
        n = size(x)
        ascend = x(n) > x(1)
        if (ascend) then
            do i = 2, n
                if (x(i) <= x(i-1)) then
                    rst = .false.
                    exit
                end if
            end do
        else
            do i = 2, n
                if (x(i) >= x(i-1)) then
                    rst = .false.
                    exit
                end if
            end do
        end if
    end function

! ------------------------------------------------------------------------------


end module
