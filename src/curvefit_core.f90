! curvefit_core.f90

!> @mainpage
!!
!! @section intro_sec Introduction
!! CURVEFIT is a library for fitting functions to sets of data.
!!
!! @author Jason Christopherson
!! @version 1.0.0

!> @brief \b curvefit_core
!!
!! @par Purpose
!! To provide core types and routines for the CURVEFIT library.
module curvefit_core
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use nonlin_types, only : NL_ARRAY_SIZE_ERROR, NL_OUT_OF_MEMORY_ERROR, &
        NL_INVALID_INPUT_ERROR
    implicit none
    private
    public :: dp
    public :: i32
    public :: CF_ARRAY_SIZE_ERROR
    public :: CF_OUT_OF_MEMORY_ERROR
    public :: CF_NO_DATA_DEFINED_ERROR
    public :: CF_INVALID_INPUT_ERROR

! ******************************************************************************
! NUMERIC TYPE CONSTANTS
! ------------------------------------------------------------------------------
    !> @brief Defines a double-precision (64-bit) floating-point type.
    integer, parameter :: dp = real64
    !> @brief Defines a 32-bit signed integer type.
    integer, parameter :: i32 = int32

! ******************************************************************************
! ERROR FLAGS
! ------------------------------------------------------------------------------
    !> An error flag denoting an improperly sized array.
    integer, parameter :: CF_ARRAY_SIZE_ERROR = NL_ARRAY_SIZE_ERROR
    !> An error denoting that there is insufficient memory available.
    integer, parameter :: CF_OUT_OF_MEMORY_ERROR = NL_OUT_OF_MEMORY_ERROR
    !> An error denoting that no data has been defined.
    integer, parameter :: CF_NO_DATA_DEFINED_ERROR = 300
    !> An error flag denoting an invalid input.
    integer, parameter :: CF_INVALID_INPUT_ERROR = NL_INVALID_INPUT_ERROR

end module
