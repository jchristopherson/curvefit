! curvefit_test_statistics.f90

module curvefit_test_statistics
    use curvefit_core
    use curvefit_statistics
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_z_value() result(rst)
        ! Local Variables
        logical :: rst
        real(dp) :: c(4), z_ans(4), z
        integer(i32) :: i

        ! Initialization
        rst = .true.
        c = [0.99d0, 0.98d0, 0.95d0, 0.90d0]
        z_ans = [2.576d0, 2.326d0, 1.96d0, 1.645d0]

        ! Compute each z test value
        do i = 1, size(c)
            z = z_value(c(i))
            if (abs(z - z_ans(i)) > 1.0d-3) then
                rst = .false.
                print '(AF5.3AF5.3A)', "Test Failed: Expected a z-value of ", &
                    z_ans(i), ", but computed ", z, "."
                exit
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    function test_mean() result(rst)
        ! Local Variables
        logical :: rst
        real(dp) :: x(5), mu

        ! Parameters
        real(dp), parameter :: ans = 10.0d0
        real(dp), parameter :: tol = 1.0d-8

        ! Initialization
        rst = .true.
        x = [9.0d0, 10.0d0, 11.0d0, 7.0d0, 13.0d0]

        ! Compute the mean
        mu = mean(x)
        if (abs(mu - ans) > tol) then
            rst = .false. 
            print '(AF5.3AF5.3A)', "Test Failed: Expected a mean value of ", &
                ans, ", but computed a value of ", mu, "."
        end if
    end function

! ------------------------------------------------------------------------------
    function test_var() result(rst)
        ! Local Variables
        logical :: rst
        real(dp) :: x(5), v

        ! Parameters
        real(dp), parameter :: ans = 5.0d0
        real(dp), parameter :: tol = 1.0d-8

        ! Initialization
        rst = .true.
        x = [9.0d0, 10.0d0, 11.0d0, 7.0d0, 13.0d0]

        ! Compute the variance
        v = variance(x)
        if (abs(v - ans) > tol) then
            rst = .false.
            print '(AF5.3AF5.3A)', "Test Failed: Expected a variance of ", &
                ans, ", but computed a value of ", v, "."
        end if
    end function

! ------------------------------------------------------------------------------
    function test_confidence_interval() result(rst)
        ! Parameters
        integer(i32), parameter :: n = 1000
        real(dp), parameter :: alpha = 0.95d0
        real(dp), parameter :: z = 1.96d0   ! 95% confidence
        real(dp), parameter :: tol = 1.0d-3 ! Do to precision of Z

        ! Local Variables
        logical :: rst
        real(dp) :: x(n), c, ans

        ! Initialization
        rst = .true.
        call random_number(x)
        ans = z * standard_deviation(x) / sqrt(real(size(x), dp))

        ! Compute the confidence interval
        c = confidence_interval(x, alpha)
        if (abs(c - ans) > tol) then
            rst = .false.
            print '(AF5.3AF5.3A)', &
                "Test Failed: Expected a confidence inverval of ", ans, &
                ", but computed a value of ", c, "."
        end if
    end function

! ------------------------------------------------------------------------------
end module
