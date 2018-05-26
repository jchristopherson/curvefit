! curvefit_stats_example_1.f90

program example
    use iso_fortran_env
    use curvefit_statistics

    ! Compute the T Value for a 95% confidence interval and 20 samples
    tval20 = t_value(0.95d0, 20)
    print '(AF8.4)', "95% Confidence Interval T Value (20 samples): ", tval20

    ! Compute the T value for a 95% confidence interval and 30 samples
    tval30 = t_value(0.95d0, 30)
    print '(AF8.4)', "95% Confidence Interval T Value (30 samples): ", tval30

    ! Compute the T value for a 95% confidence interval and 40 samples
    tval40 = t_value(0.95d0, 40)
    print '(AF8.4)', "95% Confidence Interval T Value (40 samples): ", tval40

    ! Compute the Z Value for a 95% confidence interval
    zval = z_value(0.95d0)
    print '(AF8.4)', "95% Confidence Interval Z Value: ", zval
end program
