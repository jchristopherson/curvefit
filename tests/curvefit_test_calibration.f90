! curvefit_test_calibration.f90

module curvefit_test_calibration
    use curvefit_core
    use curvefit_calibration
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_seb() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(dp), parameter :: fullscale = 5.0d2
        real(dp), parameter :: sebout_ans = 2.80010d0
        real(dp), parameter :: sebpfs_ans = 0.0104d0
        real(dp), parameter :: tol = 1.0d-4

        ! Local Variables
        real(dp), dimension(11) :: applied, measured
        type(seb_results) :: s

        ! Initialization
        rst = .true.
        applied = [0.0d0, 1.0d2, 2.0d2, 3.0d2, 4.0d2, 5.0d2, 4.0d2, 3.0d2, &
            2.0d2, 1.0d2, 0.0d0]
        measured = [0.0d0, 0.55983d0, 1.11975d0, 1.67982d0, 2.24005d0, &
            2.80039d0, 2.24023d0, 1.68021d0, 1.12026d0, 0.56021d0, 0.00006d0]
        
        ! Test
        s = seb(applied, measured, fullscale)
        if (abs(1.0d2 * s%seb / fullscale - sebpfs_ans) > tol) then
            rst = .false.
            print '(AF7.5AF7.5A)', "Test Failed: Expected an SEB of: ", &
                sebpfs_ans, "%, but computed an SEB of: ", &
                1.0d2 * s%seb / fullscale, "%."
        end if
        if (abs(s%output - sebout_ans) > tol) then
            rst = .false.
            print '(AF7.5AF7.5A)', "Test Failed: Expected an SEB Output of: ", &
                sebout_ans, " mV/V, but computed an SEB Output of: ", &
                s%output, "mV/V"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_nonlin() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(dp), parameter :: ans = -0.062122611d0
        real(dp), parameter :: slope = 178.5558182d0
        real(dp), parameter :: tol = 1.0d-4

        ! Local Variables
        real(dp), dimension(11) :: applied, output, measured
        real(dp) :: x

        ! Initialization
        rst = .true.
        applied = [0.0d0, 1.0d2, 2.0d2, 3.0d2, 4.0d2, 5.0d2, 4.0d2, 3.0d2, &
            2.0d2, 1.0d2, 0.0d0]
        output = [0.0d0, 0.55983d0, 1.11975d0, 1.67982d0, 2.24005d0, &
            2.80039d0, 2.24023d0, 1.68021d0, 1.12026d0, 0.56021d0, 0.00006d0]
        
        ! Apply the calibration
        measured = slope * output

        ! Test
        x = nonlinearity(applied, measured)
        if (abs(x - ans) > tol) then
            rst = .false.
            print '(AF7.5AF7.5)', &
                "Test Failed: Expected a nonlinearity error of: ", ans, &
                ", but computed: ", x
            return
        end if
    end function

! ------------------------------------------------------------------------------
    function test_term_nonlin() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(dp), parameter :: ans = -0.073918276d0
        real(dp), parameter :: slope = 178.5465596d0
        real(dp), parameter :: tol = 1.0d-4

        ! Local Variables
        real(dp), dimension(11) :: applied, output, measured
        real(dp) :: x

        ! Initialization
        rst = .true.
        applied = [0.0d0, 1.0d2, 2.0d2, 3.0d2, 4.0d2, 5.0d2, 4.0d2, 3.0d2, &
            2.0d2, 1.0d2, 0.0d0]
        output = [0.0d0, 0.55983d0, 1.11975d0, 1.67982d0, 2.24005d0, &
            2.80039d0, 2.24023d0, 1.68021d0, 1.12026d0, 0.56021d0, 0.00006d0]
        
        ! Apply the calibration
        measured = slope * output

        ! Test
        x = terminal_nonlinearity(applied, measured)
        if (abs(x - ans) > tol) then
            rst = .false.
            print '(AF7.5AF7.5)', &
                "Test Failed: Expected a terminal nonlinearity error of: ", &
                ans, ", but computed: ", x
            return
        end if
    end function

! ------------------------------------------------------------------------------
    function test_hysteresis() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(dp), parameter :: ans = 0.091063467
        real(dp), parameter :: slope = 178.5465596d0
        real(dp), parameter :: tol = 1.0d-4

        ! Local Variables
        real(dp), dimension(11) :: applied, output, measured
        real(dp) :: x

        ! Initialization
        rst = .true.
        applied = [0.0d0, 1.0d2, 2.0d2, 3.0d2, 4.0d2, 5.0d2, 4.0d2, 3.0d2, &
            2.0d2, 1.0d2, 0.0d0]
        output = [0.0d0, 0.55983d0, 1.11975d0, 1.67982d0, 2.24005d0, &
            2.80039d0, 2.24023d0, 1.68021d0, 1.12026d0, 0.56021d0, 0.00006d0]
        
        ! Apply the calibration
        measured = slope * output

        ! Test
        x = hysteresis(applied, measured)
        if (abs(x - ans) > tol) then
            rst = .false.
            print '(AF7.5AF7.5)', &
                "Test Failed: Expected a hysteresis error of: ", &
                ans, ", but computed: ", x
            return
        end if
    end function

! ------------------------------------------------------------------------------
end module
