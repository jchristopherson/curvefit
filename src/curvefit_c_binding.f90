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

contains
end module
