! curvefit_cal_example.f90

program example
    use curvefit_calibration
    use curvefit_core, only : dp, i32
    use curvefit_regression, only : linear_least_squares
    implicit none

    ! Local Variables
    real(dp), parameter :: fullscale = 5.0d2
    real(dp), dimension(11) :: applied, output, measured, applied_copy
    real(dp) :: hyst, gain, nlin
    type(seb_results) :: s

    ! Initialization
    applied = [0.0d0, 1.0d2, 2.0d2, 3.0d2, 4.0d2, 5.0d2, 4.0d2, 3.0d2, &
        2.0d2, 1.0d2, 0.0d0]
    output = [0.0d0, 0.55983d0, 1.11975d0, 1.67982d0, 2.24005d0, &
        2.80039d0, 2.24023d0, 1.68021d0, 1.12026d0, 0.56021d0, 0.00006d0]
    applied_copy = applied
    
    ! Determine a suitable calibration gain (the least squares routine modifies 
    ! applied; hence, the need for the copy)
    gain = linear_least_squares(output, applied_copy)

    ! Apply the calibration gain
    measured = gain * output

    ! Compute the SEB
    s = seb(applied, output, fullscale)

    ! Compute the best fit nonlinearity
    nlin = nonlinearity(applied, measured)

    ! Compute the hysteresis
    hyst = hysteresis(applied, measured)

    ! Display the results
    print '(AF9.5)', "Calibration Gain: ", gain
    print '(AF6.4)', "SEB: ", s%seb
    print '(AF7.5)', "SEB Output: ", s%output
    print '(AF6.4)', "Best Fit Nonlinearity: ", nlin
    print '(AF6.4)', "Hysteresis: ", hyst
end program