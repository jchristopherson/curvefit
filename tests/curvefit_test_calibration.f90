! curvefit_test_calibration.f90

module curvefit_test_calibration
    use curvefit_core
    use curvefit_calibration
    use curvefit_regression, only : linear_least_squares
    use test_core
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
    function test_seb_2() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(dp), parameter :: fullscale = 6.0e3
        real(dp), parameter :: sebpfs_ans = 0.0055d0
        real(dp), parameter :: tol = 1.0d-4

        ! Local Variables
        real(dp), dimension(8) :: applied, measured
        type(seb_results) :: s

        ! Initialization
        rst = .true.
        applied = [0.0d0, 1157.972206d0, 2366.64471d0, 3574.718125d0, &
            4783.30376d0, 5991.492905d0, 3574.610382d0, -0.001208501d0]
        measured = [0.0d0, 1157.719872d0, 2366.12397d0, 3574.031281d0, &
            4782.673423d0, 5990.380077d0, 3574.473523d0, 0.273213503d0]

        ! Test
        s = seb(applied, measured, fullscale)
        ! print '(AF7.5)', "SEB: ", 1.0d2 * s%seb / fullscale
        if (abs(1.0d2 * s%seb / fullscale - sebpfs_ans) > tol) then
            rst = .false.
            print '(AF7.5AF7.5A)', "Test Failed: Expected an SEB of: ", &
                sebpfs_ans, "%, but computed an SEB of: ", &
                1.0d2 * s%seb / fullscale, "%."
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
        real(dp), parameter :: ans = 0.091063467d0
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
    function test_return_to_zero() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(dp), parameter :: ans = 0.010713349d0
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
        x = return_to_zero(applied, measured)
        if (abs(x - ans) > tol) then
            rst = .false.
            print '(AF7.5AF7.5)', &
                "Test Failed: Expected a return to zero error of: ", &
                ans, ", but computed: ", x
            return
        end if
    end function

! ------------------------------------------------------------------------------
    function test_repeatability() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(i32), parameter :: npts = 15
        integer(i32), parameter :: ntests = 3
        real(dp), parameter :: slope = 2736.031907d0
        real(dp), parameter :: ans = 0.30096351d0
        real(dp), parameter :: tol = 1.0d-4

        ! Local Variables
        real(dp), dimension(npts, ntests) :: applied, output, measured
        real(dp) :: x

        ! Initialization
        rst = .true.
        applied = reshape([0.0d0, 100.0d0, 300.0d0, 600.0d0, 900.0d0, &
            1200.0d0, 1500.0d0, 1800.0d0, 2100.0d0, 2400.0d0, 2700.0d0, &
            3000.0d0, 3400.0d0, 1500.0d0, 0.0d0, 0.0d0, 100.0d0, 300.0d0, &
            600.0d0, 900.0d0, 1200.0d0, 1500.0d0, 1800.0d0, 2100.0d0, &
            2400.0d0, 2700.0d0, 3000.0d0, 3400.0d0, 1500.0d0, 0.0d0, 0.0d0, &
            100.0d0, 300.0d0, 600.0d0, 900.0d0, 1200.0d0, 1500.0d0, 1800.0d0, &
            2100.0d0, 2400.0d0, 2700.0d0, 3000.0d0, 3400.0d0, 1500.0d0, &
            0.0d0], [npts, ntests])
        output = reshape([0.0d0, 0.03654d0, 0.10961d0, 0.2192d0, 0.32883d0, &
            0.43845d0, 0.5481d0, 0.65775d0, 0.76743d0, 0.87711d0, 0.9868d0, &
            1.09646d0, 1.2427d0, 0.5481d0, 0.0d0, 0.0d0, 0.03655d0, 0.10965d0, &
            0.21925d0, 0.32888d0, 0.43852d0, 0.54818d0, 0.65784d0, 0.76752d0, &
            0.8772d0, 0.98689d0, 1.09656d0, 1.24281d0, 0.54814d0, 0.00001d0, &
            0.0d0, 0.03655d0, 0.10961d0, 0.21921d0, 0.32884d0, 0.43847d0, &
            0.54813d0, 0.65779d0, 0.76746d0, 0.87715d0, 0.98682d0, 1.09648d0, &
            1.24274d0, 0.54813d0, 0.00001d0], [npts, ntests])
        measured = slope * output

        ! Test
        x = repeatability(applied, measured)
        if (abs(x - ans) > tol) then
            rst = .false.
            print '(AF7.5AF7.5)', &
                "Test Failed: Expected a repeatability error of: ", &
                ans, ", but computed: ", x
            return
        end if
    end function

! ------------------------------------------------------------------------------
    function test_crosstalk() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(i32), parameter :: npts = 34
        integer(i32), parameter :: ndof = 2
        real(dp), parameter :: tol = 1.0d-4

        ! Local Variables
        integer(i32) :: indices(2*ndof)
        real(dp), dimension(npts, ndof) :: xin, xout, xerr, xmeas
        real(dp), dimension(ndof, npts) :: xint, xmeast
        real(dp), dimension(ndof, ndof) :: c, ans, xt

        ! Initialization
        rst = .true.
        xin = reshape([0.0, 3000.0, 6000.0, 7500.0, 9000.0, 12000.0, &
            15000.0, 7500.0, 0.0, 0.0, -3000.0, -6000.0, -7500.0, -9000.0, &
            -12000.0, -15000.0, -7500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 67.7908728, 135.5817456, 203.3726184, 271.1634912, &
            338.954364, 203.3726184, 0.0, 0.0, -67.7908728, -135.5817456, &
            -203.3726184, -271.1634912, -338.954364, -203.3726184, 0.0], &
            [npts, ndof])
        xout = reshape([0.0, 0.38905, 0.77816, 0.97269, 1.16714, 1.556, &
            1.94484, 0.9726, -1.0e-5, 0.0, -0.388886, -0.77775, -0.97215, &
            -1.16654, -1.55533, -1.9441, -0.97171, 4.0e-5, 0.0, -0.00044, &
            -0.0013, -0.0024, -0.00382, -0.00528, -0.00257, 0.00015, 0.0, &
            0.00144, 0.00306, 0.00446, 0.00567, 0.00688, 0.00451, -2.0e-5, &
            0.0, 0.00122, 0.00259, 0.0029, 0.00314, 0.00338, 0.00356, 0.00477, &
            -1.0e-5, 0.0, 0.00021, 0.00051, 0.00069, 0.00088, 0.0013, 0.00178, &
            0.00058, 3.0e-5, 0.0, 0.27156, 0.54329, 0.81507, 1.08682, 1.35881, &
            0.81553, 1.0e-5, 0.0, -0.27145, -0.54312, -0.81493, -1.0868, &
            -1.35879, -0.81548, 0.0], [npts, ndof])
        
        ! Compute the calibration gains
        xint = transpose(xin)
        xmeast = transpose(xout)
        c = linear_least_squares(xmeast, xint)
        xmeas = matmul(xout, transpose(c))
        xerr = xmeas - xin

        ! The solution is:
        ans(1,1) = 0.0d0
        ans(2,1) = nonlinearity(xin(18:34,1), xmeas(18:34,1))

        ans(1,2) = nonlinearity(xin(1:17,2), xmeas(1:17,2))
        ans(2,2) = 0.0d0

        ! The indices are
        indices = [1, 17, 18, 34]

        ! Compute the crosstalk matrix
        xt = crosstalk(xerr, indices)

        ! Test
        if (.not.is_mtx_equal(xt, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: The crosstalk error test failed."
        end if
    end function

! ------------------------------------------------------------------------------
end module
