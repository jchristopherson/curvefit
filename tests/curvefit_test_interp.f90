! curvefit_test_interp.f90

! Tests the interpolation routines.
module curvefit_test_interp
    use curvefit_core
    use curvefit_interp
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_linear_interp() result(rst)
        ! Local Variables
        logical :: rst

        ! Initialization
        rst = .true.
    end function

! ------------------------------------------------------------------------------
end module
