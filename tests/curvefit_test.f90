! curvefit_test.f90

! The testing application for the CURVEFIT library.
program main
    use curvefit_test_interp
    implicit none

    ! Local Variables
    logical :: rst, overall

    ! Initialization
    overall = .true.

    ! Interpolation Tests
    rst = test_linear_interp()
    if (.not.rst) overall = .false.

    ! End
    if (overall) then
        print '(A)', "CURVEFIT TEST STATUS: PASS"
    else
        print '(A)', "CURVEFIT TEST STATUS: FAILED"
    end if
end program
