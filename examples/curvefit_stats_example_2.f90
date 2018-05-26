! curvefit_stats_example_2.f90

program example
    use iso_fortran_env
    use curvefit_statistics
    implicit none

    ! Parameters
    integer(int32), parameter :: n = 100

    ! Variables
    real(real64) :: x(n), avg, md, std, ci

    ! Generate a random population
    call random_number(x)

    ! Compute the mean, median, standard deviation, and 95% confidence interval
    ! of the data set
    avg = mean(x)
    md = median(x)
    std = standard_deviation(x)
    ci = confidence_interval(x, 0.95d0)

    ! Print the output
    print '(AF8.4)', "Mean: ", avg
    print '(AF8.4)', "Median: ", md
    print '(AF8.4)', "Standard Deviation: ", std
    print '(AF8.4)', "95% Confidence Interval: ", ci
end program
