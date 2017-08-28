! curvefit_core.f90

!> @mainpage
!!
!! @section intro_sec Introduction
!! CURVEFIT is a ...
!!
!! @author Jason Christopherson
!! @version 1.0.0

!> @brief \b curvefit_core
!!
!! @par Purpose
!! To provide core types and routines for the CURVEFIT library.
module curvefit_core
    use, intrinsic :: iso_fortran_env, only : int32, real64
    implicit none
    private
    public :: dp
    public :: i32

! ******************************************************************************
! NUMERIC TYPE CONSTANTS
! ------------------------------------------------------------------------------
    !> @brief Defines a double-precision (64-bit) floating-point type.
    integer, parameter :: dp = real64
    !> @brief Defines a 32-bit signed integer type.
    integer, parameter :: i32 = int32

end module
