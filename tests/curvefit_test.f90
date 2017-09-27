! curvefit_test.f90

! The testing application for the CURVEFIT library.
program main
    use curvefit_test_interp
    use curvefit_test_statistics
    use curvefit_test_regression
    use curvefit_test_calibration
    implicit none

    ! Local Variables
    logical :: rst, overall

    ! Initialization
    overall = .true.

    ! Interpolation Tests
    rst = test_linear_interp()
    if (.not.rst) overall = .false.

    rst = test_poly_interp()
    if (.not.rst) overall = .false.

    rst = test_spline_interp()
    if (.not.rst) overall = .false.

    ! Statistics Tests
    rst = test_z_value()
    if (.not.rst) overall = .false.

    rst = test_mean()
    if (.not.rst) overall = .false.

    rst = test_var()
    if (.not.rst) overall = .false.

    rst = test_confidence_interval()
    if (.not.rst) overall = .false.

    rst = test_inc_gamma()
    if (.not.rst) overall = .false.

    rst = test_inc_gamma_comp()
    if (.not.rst) overall = .false.

    rst = test_covariance()
    if (.not.rst) overall = .false.

    rst = test_covariance_2()
    if (.not.rst) overall = .false.

    rst = test_median()
    if (.not.rst) overall = .false.

    ! Regression Tests
    rst = test_lowess()
    if (.not.rst) overall = .false.

    ! Calibration Tests
    rst = test_seb()
    if (.not.rst) overall = .false.

    ! End
    if (overall) then
        print '(A)', "CURVEFIT TEST STATUS: PASS"
    else
        print '(A)', "CURVEFIT TEST STATUS: FAILED"
    end if
end program
