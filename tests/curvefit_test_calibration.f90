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
        print '(AF8.5A)', "SEB: ", 1.0d2 * s%seb / fullscale, "%"
        print '(AF7.5A)', "Output: ", s%output, " mV/V"
        print '(AF10.5A)', "Slope: ", s%slope, " N/mV/V"

    end function

! ------------------------------------------------------------------------------
end module
