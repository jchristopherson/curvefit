! curvefit_regression.f90

!> @brief \b curvefit_regression
!!
!! @par Purpose
!! To provide routines for perforing regression operations, and other data
!! smoothing operations on sets of numerical data.
module curvefit_regression
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use curvefit_core
    use ferror, only : errors
    use nonlin_core, only : vecfcn_helper, iteration_behavior
    use nonlin_least_squares, only : least_squares_solver
    use curvefit_statistics, only : mean
    use linalg_core
    use nonlin_polynomials
    implicit none
    private
    public :: moving_average
    public :: lowess_smoothing
    public :: nonlinear_regression
    public :: linear_least_squares

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Applies a moving average to smooth a data set.
    !!
    !! @par Example
    !! The following example illustrates the use of a moving average to smooth
    !! a noisy data set.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use curvefit_regression
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: n = 1000
    !!     integer(int32), parameter :: nAvg = 20
    !!     real(real64), parameter :: xMin = 0.0d0
    !!     real(real64), parameter :: xMax = 1.0d0
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: d1, d2
    !!
    !!     ! Local Variables
    !!     real(real64) :: x(n), y(n), yAvg(n)
    !!
    !!     ! Initialization
    !!     x = linspace(xMin, xMax, n)
    !!     call random_number(y)
    !!     y = (y - 0.5d0) + &
    !!         sin(15.0d0 * x) + &
    !!         0.5d0 * sin(25.0d0 * x) + &
    !!         0.1d0 * sin(75.0d0 * x)
    !!     yAvg = y
    !!
    !!     ! Apply the moving average to the data set
    !!     call moving_average(yAvg, nAvg)
    !!
    !!     ! Plot the results
    !!     call plt%initialize()
    !!     call plt%set_font_size(14)
    !!
    !!     call d1%set_name("Noisy Data")
    !!     call d1%set_line_color(CLR_BLUE)
    !!     call d1%define_data(x, y)
    !!
    !!     call d2%set_name("Averaged")
    !!     call d2%set_line_color(CLR_RED)
    !!     call d2%set_line_width(3.0)
    !!     call d2%define_data(x, yAvg)
    !!
    !!     call plt%push(d1)
    !!     call plt%push(d2)
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @image html moving_average_example.png
    interface moving_average
        module procedure :: moving_average_1
    end interface

! ------------------------------------------------------------------------------
    !> @brief Employs a least squares fit to determine the coefficient A in the
    !! linear system: Y = A * X, where A can either be a scalar, or a matrix.
    interface linear_least_squares
        module procedure :: linear_least_squares_1var
        module procedure :: linear_least_squares_nvar
    end interface

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Defines a type for computing a smoothing of an X-Y data set using
    !! a robust locally weighted scatterplot smoothing (LOWESS) algorithm.
    !!
    !! @par Example
    !! The following example utilizes LOWESS smoothing on a data set.  The
    !! example compares various levels of smoothing, and compares with
    !! nonlinear regression.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use curvefit_regression
    !!     use curvefit_core
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: n = 100
    !!     real(real64), parameter :: maxX = 1.0d0
    !!     real(real64), parameter :: minX = 0.0d0
    !!     type(color), parameter :: orange = color(255, 102, 0)
    !!
    !!     ! Local Variables
    !!     integer(int32) :: i
    !!     real(real64) :: x(n), y(n), yr(n), ys(n), ys2(n), cnl(5), ynl(n)
    !!     type(lowess_smoothing) :: fit
    !!     type(nonlinear_regression) :: solver
    !!     procedure(reg_fcn), pointer :: fcn
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: d1, d2, d3, d4, d5
    !!
    !!     ! Initialization
    !!     x = linspace(minX, maxX, n)
    !!     y = 0.5d0 * sin(2.0d1 * x) + cos(5.0d0 * x) * exp(-0.1d0 * x)
    !!     call random_number(yr)
    !!     yr = y + (yr - 0.5d0)
    !!
    !!     ! Generate the fit
    !!     call fit%initialize(x, yr)
    !!     ys = fit%smooth(0.2d0)
    !!     ys2 = fit%smooth(0.8d0)
    !!
    !!     ! For comparison purposes, consider a nonlinear regression fit.  As we know
    !!     ! the coefficients, they provide a very good starting guess.
    !!     cnl = [0.5d0, 2.0d0, 20.0d0, 5.0d0, -0.1d0]
    !!     fcn => nrfun
    !!     call solver%initialize(x, yr, fcn, size(cnl))
    !!     call solver%solve(cnl)
    !!     do i = 1, n
    !!         ynl(i) = fcn(x(i), cnl)
    !!     end do
    !!
    !!     ! Display the computed coefficients
    !!     print '(A)', "f(x) = c0 * sin(c1 * x) + c2 * cos(c3 * x) * exp(c4 * x):"
    !!     print '(AF12.10)', "c0: ", cnl(1)
    !!     print '(AF13.10)', "c1: ", cnl(2)
    !!     print '(AF12.10)', "c2: ", cnl(3)
    !!     print '(AF12.10)', "c3: ", cnl(4)
    !!     print '(AF13.10)', "c4: ", cnl(5)
    !!
    !!     ! Plot the data
    !!     call plt%initialize()
    !!     call plt%set_font_size(14)
    !!
    !!     call d1%set_name("Original Signal")
    !!     call d1%set_line_color(CLR_BLUE)
    !!     call d1%set_line_width(2.0)
    !!     call d1%define_data(x, y)
    !!
    !!     call d2%set_name("Raw Data")
    !!     call d2%set_line_color(CLR_BLACK)
    !!     call d2%set_draw_line(.false.)
    !!     call d2%set_draw_markers(.true.)
    !!     call d2%set_marker_style(MARKER_EMPTY_CIRCLE)
    !!     call d2%set_marker_scaling(2.0)
    !!     call d2%set_line_width(2.0)
    !!     call d2%define_data(x, yr)
    !!
    !!     call d3%set_name("Smoothed (f = 0.2)")
    !!     call d3%set_line_color(CLR_GREEN)
    !!     call d3%set_line_width(2.0)
    !!     call d3%set_line_style(LINE_DASHED)
    !!     call d3%define_data(x, ys)
    !!
    !!     call d4%set_name("Smoothed (f = 0.8)")
    !!     call d4%set_line_color(CLR_MAGENTA)
    !!     call d4%set_line_width(2.0)
    !!     call d4%set_line_style(LINE_DASH_DOTTED)
    !!     call d4%define_data(x, ys2)
    !!
    !!     call d5%set_name("Nonlinear Regression")
    !!     call d5%set_line_color(orange)
    !!     call d5%set_line_width(2.0)
    !!     call d5%set_line_style(LINE_DASH_DOT_DOT)
    !!     call d5%define_data(x, ynl)
    !!
    !!     call plt%push(d1)
    !!     call plt%push(d2)
    !!     call plt%push(d3)
    !!     call plt%push(d4)
    !!     call plt%push(d5)
    !!
    !!     call plt%draw()
    !!
    !! contains
    !!     function nrfun(xp, c) result(fn)
    !!         real(real64), intent(in) :: xp
    !!         real(real64), intent(in), dimension(:) :: c
    !!         real(real64) :: fn
    !!         fn = c(1) * sin(c(2) * xp) + c(3) * cos(c(4) * xp) * exp(c(5) * xp)
    !!     end function
    !!
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @image html lowess_example_1.png
    !!
    !! @par
    !! Notice, as the data points are generated randomly, it is possible to
    !! obtain a different output than illustrated above.  The following output
    !! was generated by simply running the same program again.
    !! @image html lowess_example_2.png
    type lowess_smoothing
        private
        !> N-element array of x data points - sorted into ascending order.
        real(real64), allocatable, dimension(:) :: m_x
        !> N-element array of y data points.
        real(real64), allocatable, dimension(:) :: m_y
        !> N-element array containing the robustness weights for each data
        !! point.
        real(real64), allocatable, dimension(:) :: m_weights
        !> N-element array containing the residuals (Y - YS)
        real(real64), allocatable, dimension(:) :: m_residuals
        !> Scaling parameter used to define the nature of the linear
        !! interpolations used by the algorithm.
        real(real64) :: m_delta
        !> Tracks whether or not ls_init has been called
        logical :: m_init = .false.
    contains
        !> @brief Initializes the lowess_smoothing object.
        procedure, public :: initialize => ls_init
        !> @brief Performs the actual smoothing operation.
        procedure, public :: smooth => ls_smooth
        !> @brief Gets the number of stored data points.
        procedure, public :: get_count => ls_get_num_pts
        !> @brief Gets the x component of the requested data point.
        procedure, public :: get_x => ls_get_x
        !> @brief Gets the y component of the requested data point.
        procedure, public :: get_y => ls_get_y
        !> @brief Gets the residuals from each data point.
        procedure, public :: get_residuals => ls_get_residual
    end type

! ------------------------------------------------------------------------------
    !> @brief A type for supporting nonlinear regression calculations.
    !!
    !! @par Example
    !! The following example illustrates the use of the nonlinear_regression
    !! type to fit an equation to a noisy set of data.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use curvefit_regression
    !!     use curvefit_core
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: npts = 100
    !!     real(real64), parameter :: xMin = 0.0d0
    !!     real(real64), parameter :: xMax = 1.0d0
    !!     real(real64), parameter :: a = 2.0d0
    !!     real(real64), parameter :: w = 30.0d0
    !!     real(real64), parameter :: d = 0.2d0
    !!     real(real64), parameter :: p = 0.0d0
    !!
    !!     ! Local variables
    !!     procedure(reg_fcn), pointer :: fun
    !!     real(real64) :: x(npts), y(npts), c(4), yfit(npts)
    !!     type(nonlinear_regression) :: solver
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: d1, d2
    !!     integer(int32) :: i
    !!
    !!     ! Define the data to fit
    !!     x = linspace(xMin, xMax, npts)
    !!     call random_number(y)
    !!     y = a * exp(-w * d * x) * sin(w * x - p) + (y - 0.5d0)
    !!
    !!     ! Set up the solver.  Assume an equation of the form:
    !!     !   A exp(-w d x) sin(w * x - p)
    !!     fun => fcn
    !!     c = [1.0d0, 20.0d0, 0.0d0, 0.0d0]
    !!     call solver%initialize(x, y, fun, size(c))
    !!
    !!     ! Solve
    !!     call solver%solve(c)
    !!
    !!     ! Compute the fitted form of the routine
    !!     do i = 1, size(x)
    !!         yfit(i) = fun(x(i), c)
    !!     end do
    !!
    !!     ! Print out the coefficients
    !!     print '(AF8.4AF8.4)', "A = ", c(1), ", Actual = ", a
    !!     print '(AF8.4AF8.4)', "w = ", c(2), ", Actual = ", w
    !!     print '(AF8.4AF8.4)', "d = ", c(3), ", Actual = ", d
    !!     print '(AF8.4AF8.4)', "p = ", c(4), ", Actual = ", p
    !!
    !!     ! Plot the results
    !!     call plt%initialize()
    !!     call plt%set_font_size(14)
    !!
    !!     call d1%set_name("Raw Data")
    !!     call d1%set_line_color(CLR_BLACK)
    !!     call d1%set_draw_line(.false.)
    !!     call d1%set_draw_markers(.true.)
    !!     call d1%set_marker_style(MARKER_X)
    !!     call d1%define_data(x, y)
    !!
    !!     call d2%set_name("Fit")
    !!     call d2%set_line_color(CLR_BLUE)
    !!     call d2%set_line_width(2.0)
    !!     call d2%define_data(x, yfit)
    !!
    !!     call plt%push(d1)
    !!     call plt%push(d2)
    !!     call plt%draw()
    !!
    !! contains
    !!     function fcn(x, c) result(f)
    !!         real(real64), intent(in) :: x
    !!         real(real64), intent(in), dimension(:) :: c
    !!         real(real64) :: f
    !!
    !!         f = c(1) * exp(-c(2) * c(3) * x) * sin(c(2) * x - c(4))
    !!     end function
    !! end program
    !! @endcode
    !! The above program produces the following outputs.
    !! @code{.txt}
    !! A =   1.9637, Actual =   2.0000
    !! w =  30.6889, Actual =  30.0000
    !! d =   0.1947, Actual =   0.2000
    !! p =   0.1400, Actual =   0.0000
    !! @endcode
    !! @image html nlreg_example_2.png
    type, extends(vecfcn_helper) :: nonlinear_regression
        private
        !> A pointer to the routine containing the function of interest.
        procedure(reg_fcn), pointer, nopass :: m_rfcn => null()
        !> The x data points.
        real(real64), allocatable, dimension(:) :: m_x
        !> The y data points.
        real(real64), allocatable, dimension(:) :: m_y
        !> The number of coefficients in the function of interest
        integer(int32) :: m_ncoeff = 0
        !> Tracks whether or not nr_init has been called
        logical :: m_init = .false.
        !> The Levenberg-Marquardt solver
        type(least_squares_solver) :: m_solver
    contains
        !> @brief Initializes the nonlinear_regression object.
        procedure, public :: initialize => nr_init
        !> @brief Computes the residual between the supplied data set, and the
        !! function value given a set of coefficients.
        procedure, public :: fcn => nr_fcn
        !> @brief Determines if the function has been defined.
        procedure, public :: is_fcn_defined => nr_is_fcn_defined
        !> @brief Gets the number of equations required to solve the regression
        !! problem.
        procedure, public :: get_equation_count => nr_get_eqn_count
        !> @brief Gets the number of variables (coefficients).
        procedure, public :: get_variable_count => nr_get_var_count
        !> @brief Computes the solution to the nonlinear regression problem
        !! using the Levenberg-Marquardt method.
        procedure, public :: solve => nr_solve
        !> @brief Gets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: get_max_fcn_evals => nr_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: set_max_fcn_evals => nr_set_max_eval
        !> @brief Gets the convergence on function value tolerance.
        procedure, public :: get_fcn_tolerance => nr_get_fcn_tol
        !> @brief Sets the convergence on function value tolerance.
        procedure, public :: set_fcn_tolerance => nr_set_fcn_tol
        !> @brief Gets the convergence on change in variable tolerance.
        procedure, public :: get_var_tolerance => nr_get_var_tol
        !> @brief Sets the convergence on change in variable tolerance.
        procedure, public :: set_var_tolerance => nr_set_var_tol
        !> @brief Gets the convergence on slope of the gradient vector
        !! tolerance.
        procedure, public :: get_gradient_tolerance => nr_get_grad_tol
        !> @brief Sets the convergence on slope of the gradient vector
        !! tolerance.
        procedure, public :: set_gradient_tolerance => nr_set_grad_tol
        !> @brief Gets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: get_print_status => nr_get_print_status
        !> @brief Sets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: set_print_status => nr_set_print_status
        !> @brief Gets the number of stored data points.
        procedure, public :: get_count => nr_get_num_pts
        !> @brief Gets the x component of the requested data point.
        procedure, public :: get_x => nr_get_x
        !> @brief Gets the y component of the requested data point.
        procedure, public :: get_y => nr_get_y
    end type


contains
! ******************************************************************************
! MISC. PUBLIC ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Applies a moving average to smooth a data set.
    !!
    !! @param[in,out] x On input, the signal to smooth.  On output, the smoothed
    !!  signal.
    !! @param[in] npts The size of the averaging window.  This value must be
    !!  at least 2, but no more than the number of elements in @p x.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p npts is less than 2, or greater
    !!      than the length of @p x.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    subroutine moving_average_1(x, npts, err)
        ! Arguments
        real(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(in) :: npts
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: i, n, flag
        real(real64), allocatable, dimension(:) :: buffer

        ! Initialization
        n = size(x)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (npts < 2 .or. npts > n) then
            call errmgr%report_error("moving_average_1", &
                "The averaging size window must be at least 2, and less " // &
                "than the size of the data set.", CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(buffer(npts), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("moving_average_1", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Process
        buffer = zero
        do i = 1, n
            ! Index the buffer
            buffer(2:npts) = buffer(1:npts-1)

            ! Add a new sample value to the buffer
            buffer(1) = x(i)

            ! Compute the mean
            x(i) = mean(buffer)
        end do
    end subroutine

! ******************************************************************************
! LOCAL REGRESSION - LOWESS
! ------------------------------------------------------------------------------
    !> @brief A support routine for the LOWESS library used to compute the
    !! smoothing of a desired value from a data set.
    !!
    !! @param[in] x An N-element containing the independent variable values of
    !!  the data set.  This array must be in a monotonically increasing order.
    !! @param[in] y  An N-element array of the dependent variables corresponding
    !!  to @p x.
    !! @param[in] xs The value of the independent variable at which the
    !!  smoothing is computed.
    !! @param[out] ys The fitted value.
    !! @param[in] nleft The index of the first point which should be considered
    !!  in computing the fit.
    !! @param[in] nright The index of the last point which should be considered
    !!  in computing the fit.
    !! @param[out] w An N-element array that, on output, contains the weights
    !!  for @p y in the expression for @p ys.
    !! @param[in] userw  If true, a robust fit is carried out using the weights
    !!  in @p rw.  If false, the values in @p rw are not used.
    !! @param[in] rw An N-element array containing the robustness weights.
    !! @param[out] ok Returns true if the calculations were performed; however,
    !!  returns false if the weights are all zero-valued.
    !!
    !! @par Remarks
    !! This routines is an implementation of the LOWEST routine from the LOWESS
    !! library.  A link to this library, along with a basic description of the
    !! algorithm is available
    !! [here](https://en.wikipedia.org/wiki/Local_regression).  For a detailed
    !! understanding of the algorithm, see the [paper]
    !! (http://www.aliquote.org/cours/2012_biomed/biblio/Cleveland1979.pdf) by
    !! William Cleveland.
    subroutine lowest(x, y, xs, ys, nleft, nright, w, userw, rw, ok)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x, y, rw ! N ELEMENT
        real(real64), intent(in) :: xs
        real(real64), intent(out) :: ys
        integer(int32), intent(in) :: nleft, nright
        real(real64), intent(out), dimension(:) :: w ! N ELEMENT
        logical, intent(in) :: userw
        logical, intent(out) :: ok

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: p001 = 1.0d-3
        real(real64), parameter :: p999 = 0.999d0

        ! Local Variables
        integer(int32) :: j, n, nrt
        real(real64) :: range, h, h9, h1, a, b, c, r

        ! Initialization
        n = size(x)
        range = x(n) - x(1)
        h = max(xs - x(nleft), x(nright) - xs)
        h9 = p999 * h
        h1 = p001 * h
        a = zero

        ! Process
        do j = nleft, n
            w(j) = zero
            r = abs(x(j) - xs)
            if (r <= h9) then
                if (r > h1) then
                    w(j) = (one - (r / h)**3)**3
                else
                    w(j) = one
                end if
                if (userw) w(j) = rw(j) * w(j)
                a = a + w(j)
            else if (x(j) > xs) then
                exit
            end if
        end do

        nrt = j - 1
        if (a <= zero) then
            ok = .false.
        else
            ok = .true.
            w(nleft:nrt) = w(nleft:nrt) / a
            if (h > zero) then
                a = zero
                do j = nleft, nrt
                    a = a + w(j) * x(j)
                end do
                b = xs - a
                c = zero
                do j = nleft, nrt
                    c = c + w(j) * (x(j) - a)**2
                end do
                if (sqrt(c) > p001 * range) then
                    b = b / c
                    do j = nleft, nrt
                        w(j) = w(j) * (one + b * (x(j) - a))
                    end do
                end if
            end if
            ys = zero
            do j = nleft, nrt
                ys = ys + w(j) * y(j)
            end do
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes a smoothing of an X-Y data set using a robust locally
    !! weighted scatterplot smoothing (LOWESS) algorithm.  Fitted values are
    !! computed at each of the supplied x values.
    !!
    !! @param[in] x An N-element containing the independent variable values of
    !!  the data set.  This array must be in a monotonically increasing order.
    !! @param[in] y  An N-element array of the dependent variables corresponding
    !!  to @p x.
    !! @param[in] f Specifies the amount of smoothing.  More specifically, this
    !! value is the fraction of points used to compute each value.  As this
    !! value increases, the output becomes smoother.  Choosing a value in the
    !! range of 0.2 to 0.8 usually results in a good fit.  As such, a reasonable
    !! starting point, in the absence of better information, is a value of 0.5.
    !! @param[in] nsteps The number of iterations in the robust fit.  If set to
    !!  zero, a nonrobust fit is returned.  Seeting this parameter equal to 2
    !!  should serve most purposes.
    !! @param[in] delta A nonnegative parameter which may be used to save
    !!  computations.  If N is less than 100, set delta equal to 0.0.  If N is
    !!  larger than 100, set delta = range(x) / k, where k determines the
    !!  interpolation window used by the linear weighted regression
    !!  computations.
    !! @param[out] ys An N-element array that, on output, contains the fitted
    !!  values.
    !! @param[out] rw  An N-element array that, on output, contains the
    !!  robustness weights given to each data point.
    !! @param[out] rs An N-element array that, on output, contains the residual
    !!  @p y - @p ys.
    !!
    !! @par Remarks
    !! This routines is an implementation of the LOWESS routine from the LOWESS
    !! library.  A link to this library, along with a basic description of the
    !! algorithm is available
    !! [here](https://en.wikipedia.org/wiki/Local_regression).  For a detailed
    !! understanding of the algorithm, see the [paper]
    !! (http://www.aliquote.org/cours/2012_biomed/biblio/Cleveland1979.pdf) by
    !! William Cleveland.
    subroutine lowess(x, y, f, nsteps, delta, ys, rw, res)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x, y
        real(real64), intent(in) :: f
        integer(int32), intent(in) :: nsteps
        real(real64), intent(in) :: delta
        real(real64), intent(out), dimension(:) :: ys, rw, res

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: three = 3.0d0
        real(real64), parameter :: p001 = 1.0d-3
        real(real64), parameter :: p999 = 0.999d0

        ! Local Variables
        logical :: ok
        integer(int32) :: iter, i, j, n, nleft, nright, ns, last, m1, m2
        real(real64) :: d1, d2, denom, alpha, cut, eps, cmad, c1, c9, r

        ! Initialization
        n = size(x)
        ns = max(min(int(f * real(n, real64), int32), n), 2)
        eps = epsilon(eps)

        ! Quick Return
        if (n < 2) then
            ys = y
            return
        end if

        ! Process
        do iter = 1, nsteps + 1
            nleft = 1
            nright = ns
            last = 0
            i = 1
            do
                do while (nright < n)
                    d1 = x(i) - x(nleft)
                    d2 = x(nright+1) - x(i)
                    if (d1 <= d2) exit
                    nleft = nleft + 1
                    nright = nright + 1
                end do

                call lowest(x, y, x(i), ys(i), nleft, nright, res, iter > 1, &
                    rw, ok)
                if (.not.ok) ys(i) = y(i)
                if (last < i - 1) then
                    denom = x(i) - x(last)
                    do j = last + 1, i - 1
                        alpha = (x(j) - x(last)) / denom
                        ys(j) = alpha * ys(i) + (one - alpha) * ys(last)
                    end do
                end if
                last = i
                cut = x(last) + delta
                do i = last + 1, n
                    if (x(i) > cut) exit
                    if (abs(x(i) - x(last)) < eps) then
                        ys(i) = ys(last)
                        last = i
                    end if
                end do
                i = max(last + 1, i - 1)

                if (last >= n) exit
            end do

            res = y - ys
            if (iter > nsteps) exit
            rw = abs(res)
            call sort(rw, .true.)
            m1 = 1 + n / 2
            m2 = n - m1 + 1
            cmad = three * (rw(m1) + rw(m2))
            c9 = p999 * cmad
            c1 = p001 * cmad
            do i = 1, n
                r = abs(res(i))
                if (r <= c1) then
                    rw(i) = one
                else if (r > c9) then
                    rw(i) = zero
                else
                    rw(i) = (one - (r / cmad)**2)**2
                end if
            end do
        end do
    end subroutine

! ******************************************************************************
! LOWESS_SMOOTHING MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Initializes the lowess_smoothing object.
    !!
    !! @param[in,out] this The lowess_smoothing object.
    !! @param[in] x An N-element containing the independent variable values of
    !!  the data set.  This array must be in a monotonically increasing order.
    !!  The routine is capable of sorting the array into ascending order,
    !!  dependent upon the value of @p srt.  If sorting is performed, this
    !!  routine will also shuffle @p y to match.
    !! @param[in] y  An N-element array of the dependent variables corresponding
    !!  to @p x.
    !! @param[in] srt An optional flag determining if @p x should be sorted.
    !!  The default is to sort (true).
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    subroutine ls_init(this, x, y, srt, err)
        ! Arguments
        class(lowess_smoothing), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: x, y
        logical, intent(in), optional :: srt
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, n, flag
        integer(int32), allocatable, dimension(:) :: indices
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        logical :: sortData

        ! Initialization
        this%m_init = .false.
        n = size(x)
        sortData = .true.
        if (present(srt)) sortData = srt
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(y) /= n) then
            call errmgr%report_error("ls_init", &
                "Input array sizes must match.", CF_ARRAY_SIZE_ERROR)
            return
        end if

        ! Memory Allocations
        if (allocated(this%m_x)) deallocate(this%m_x)
        if (allocated(this%m_y)) deallocate(this%m_y)
        if (allocated(this%m_weights)) deallocate(this%m_weights)
        if (allocated(this%m_residuals)) deallocate(this%m_residuals)
        allocate(this%m_x(n), stat = flag)
        if (flag == 0) allocate(this%m_y(n), stat = flag)
        if (flag == 0) allocate(this%m_weights(n), stat = flag)
        if (flag == 0) allocate(this%m_residuals(n), stat = flag)
        if (flag == 0 .and. sortData) allocate(indices(n), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("ls_init", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Copy over the data
        if (sortData) then
            do concurrent (i = 1:n)
                this%m_x(i) = x(i)
                indices(i) = i
            end do
            call sort(this%m_x, indices, .true.)
            this%m_y = y(indices)
        else
            do concurrent (i = 1:n)
                this%m_x(i) = x(i)
                this%m_y(i) = y(i)
            end do
        end if

        ! Additional Initialization
        this%m_delta = zero
        if (n > 100) then
        end if
        this%m_init = .true.
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Performs the actual smoothing operation.
    !!
    !! @param[in,out] this The lowess_smoothing object.
    !! @param[in] f Specifies the amount of smoothing.  More specifically, this
    !! value is the fraction of points used to compute each value.  As this
    !! value increases, the output becomes smoother.  Choosing a value in the
    !! range of 0.2 to 0.8 usually results in a good fit.  As such, a reasonable
    !! starting point, in the absence of better information, is a value of 0.5.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_NO_DATA_DEFINED_ERROR: Occurs if no data has been defined.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!
    !! @return The smoothed data points.
    function ls_smooth(this, f, err) result(ys)
        ! Arguments
        class(lowess_smoothing), intent(inout) :: this
        real(real64), intent(in) :: f
        class(errors), intent(inout), optional, target :: err
        real(real64), allocatable, dimension(:) :: ys

        ! Local Variables
        integer(int32) :: n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Input Check
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        if (.not.this%m_init) then
            ! ERROR
            call errmgr%report_error("ls_smooth", &
                "No data has been defined.", CF_NO_DATA_DEFINED_ERROR)
            return
        end if
        n = size(this%m_x)

        ! Process
        allocate(ys(n), stat = flag)
        if (flag /= 0) then
            ! ERROR
            call errmgr%report_error("ls_smooth", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if
        call lowess(this%m_x, this%m_y, f, 2, this%m_delta, ys, &
            this%m_weights, this%m_residuals)
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the number of stored data points.
    !!
    !! @param[in] this The lowess_smoothing object.
    !!
    !! @return The number of data points.
    pure function ls_get_num_pts(this) result(n)
        class(lowess_smoothing), intent(in) :: this
        integer(int32) :: n
        if (this%m_init) then
            n = size(this%m_x)
        else
            n = 0
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the x component of the requested data point.
    !!
    !! @param[in] this The lowess_smoothing object.
    !! @param[in] ind The one-based index of the data point to retrieve.
    !!
    !! @return The x component of the requested data point.
    pure function ls_get_x(this, ind) result(x)
        class(lowess_smoothing), intent(in) :: this
        integer(int32), intent(in) :: ind
        real(real64) :: x
        if (this%m_init) then
            x = this%m_x(ind)
        else
            x = 0.0d0
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the y component of the requested data point.
    !!
    !! @param[in] this The lowess_smoothing object.
    !! @param[in] ind The one-based index of the data point to retrieve.
    !!
    !! @return The y component of the requested data point.
    pure function ls_get_y(this, ind) result(y)
        class(lowess_smoothing), intent(in) :: this
        integer(int32), intent(in) :: ind
        real(real64) :: y
        if (this%m_init) then
            y = this%m_y(ind)
        else
            y = 0.0d0
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the residuals from each data point.
    !!
    !! @param[in] this The lowess_smoothing object.
    !! @param[out] x An N-element array where the residual data should be
    !!  written.
    subroutine ls_get_residual(this, x)
        ! Arguments
        class(lowess_smoothing), intent(in) :: this
        real(real64), intent(out), dimension(:) :: x

        ! Local Variables
        integer(int32) :: n

        ! Process
        n = min(size(x), this%get_count())
        x(1:n) = this%m_residuals(1:n)
    end subroutine

! ******************************************************************************
! nonlinear_regression MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Initializes the nonlinear_regression object.
    !!
    !! @param[in,out] this The nonlinear_regression object.
    !! @param[in] x An N-element containing the independent variable values of
    !!  the data set.
    !! @param[in] y  An N-element array of the dependent variables corresponding
    !!  to @p x.
    !! @param[in] fcn A pointer to the function whose coefficients are to be
    !!  determined.
    !! @param[in] ncoeff The number of coefficients in the function defined in
    !!  @p fcn.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p ncoeff is less than or equal to
    !!      zero.
    subroutine nr_init(this, x, y, fcn, ncoeff, err)
        ! Arguments
        class(nonlinear_regression), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: x, y
        procedure(reg_fcn), pointer, intent(in) :: fcn
        integer(int32), intent(in) :: ncoeff
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: i, n, flag

        ! Initialization
        this%m_init = .false.
        n = size(x)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(y) /= n) then
            call errmgr%report_error("nr_init", &
                "The input arrays must be the same size.", CF_ARRAY_SIZE_ERROR)
            return
        end if
        if (ncoeff <= 0) then
            call errmgr%report_error("nr_init", &
                "The number of equation coefficients must be positive.", &
                CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Allocate memory
        if (allocated(this%m_x)) deallocate(this%m_x)
        if (allocated(this%m_y)) deallocate(this%m_y)
        allocate(this%m_x(n), stat = flag)
        if (flag == 0) allocate(this%m_y(n), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("nr_init", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Process
        do concurrent (i = 1:n)
            this%m_x(i) = x(i)
            this%m_y(i) = y(i)
        end do
        this%m_rfcn => fcn
        this%m_ncoeff = ncoeff
        this%m_init = .true.
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the residual between the supplied data set, and the
    !! function value given a set of coefficients.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @param[in] x An N-element array containing the N coefficients.
    !! @param[out] f An M-element array that, on output, contains the residual
    !!  at each of the M data points.
    subroutine nr_fcn(this, x, f)
        ! Arguments
        class(nonlinear_regression), intent(in) :: this
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f

        ! Local Variables
        integer(int32) :: i, n

        ! Compute the value of the function at each value of m_x
        n = size(this%m_x)
        do i = 1, n
            f(i) = this%m_y(i) - this%m_rfcn(this%m_x(i), x)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Determines if the function has been defined.
    !!
    !! @param[in] this The nonlinear_regression object.
    !!
    !! @return Returns true if the function has been defined; else, false.
    pure function nr_is_fcn_defined(this) result(x)
        class(nonlinear_regression), intent(in) :: this
        logical :: x
        x = associated(this%m_rfcn)
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the number of equations required to solve the regression
    !! problem.
    !!
    !! @param[in] this The nonlinear_regression object.
    !!
    !! @return The number of equations.
    pure function nr_get_eqn_count(this) result(n)
        class(nonlinear_regression), intent(in) :: this
        integer(int32) :: n
        if (allocated(this%m_x)) then
            n = size(this%m_x)
        else
            n = 0
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the number of variables (coefficients).
    !!
    !! @param[in] this The nonlinear_regression object.
    !!
    !! @return The number of variables.
    pure function nr_get_var_count(this) result(n)
        class(nonlinear_regression), intent(in) :: this
        integer(int32) :: n
        n = this%m_ncoeff
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the solution to the nonlinear regression problem using
    !! the Levenberg-Marquardt method.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @param[in,out] c On input, an array containing initial estimates of the
    !!  coefficients.  On output, the comptued coefficient values.
    !! @param[out] res An optional output array, whose size corresponds to the
    !!  number of data points, that can be used to retrieve the residual error
    !!  at each data point.
    !! @param[out] ib An optional output, that if provided, allows the
    !!  caller to obtain iteration performance statistics.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if the number of equations is less than
    !!      than the number of variables.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      correctly.
    !!  - CF_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
    !!      the allowed number of iterations.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - CF_TOLERANCE_TOO_SMALL_ERROR: Occurs if the requested tolerance is
    !!      to small to be practical for the problem at hand.
    subroutine nr_solve(this, c, res, ib, err)
        ! Arguments
        class(nonlinear_regression), intent(inout) :: this
        real(real64), intent(inout), dimension(:) :: c
        real(real64), intent(out), dimension(:), target, optional :: res
        type(iteration_behavior), intent(out), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: n, flag
        real(real64), allocatable, target, dimension(:) :: f
        real(real64), pointer, dimension(:) :: fptr

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (.not.this%m_init) then
        end if
        if (size(c) /= this%get_variable_count()) then
        end if

        ! Local Memory Allocation
        n = this%get_equation_count()
        if (present(res)) then
            if (size(res) /= n) then
            end if
            fptr => res
        else
            allocate(f(n), stat = flag)
            if (flag /= 0) then
            end if
            fptr => f
        end if

        ! Compute the solution
        call this%m_solver%solve(this, c, fptr, ib = ib, err = errmgr)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return The maximum number of function evaluations.
    pure function nr_get_max_eval(this) result(n)
        class(nonlinear_regression), intent(in) :: this
        integer(int32) :: n
        n = this%m_solver%get_max_fcn_evals()
    end function

! --------------------
    !> @brief Sets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in,out] this The nonlinear_regression object.
    !! @param[in] n The maximum number of function evaluations.
    subroutine nr_set_max_eval(this, n)
        class(nonlinear_regression), intent(inout) :: this
        integer(int32), intent(in) :: n
        call this%m_solver%set_max_fcn_evals(n)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on function value tolerance.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return The tolerance value.
    pure function nr_get_fcn_tol(this) result(x)
        class(nonlinear_regression), intent(in) :: this
        real(real64) :: x
        x = this%m_solver%get_fcn_tolerance()
    end function

! --------------------
    !> @brief Sets the convergence on function value tolerance.
    !!
    !! @param[in,out] this The nonlinear_regression object.
    !! @param[in] x The tolerance value.
    subroutine nr_set_fcn_tol(this, x)
        class(nonlinear_regression), intent(inout) :: this
        real(real64), intent(in) :: x
        call this%m_solver%set_fcn_tolerance(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on change in variable tolerance.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return The tolerance value.
    pure function nr_get_var_tol(this) result(x)
        class(nonlinear_regression), intent(in) :: this
        real(real64) :: x
        x = this%m_solver%get_var_tolerance()
    end function

! --------------------
    !> @brief Sets the convergence on change in variable tolerance.
    !!
    !! @param[in,out] this The nonlinear_regression object.
    !! @param[in] x The tolerance value.
    subroutine nr_set_var_tol(this, x)
        class(nonlinear_regression), intent(inout) :: this
        real(real64), intent(in) :: x
        call this%m_solver%set_var_tolerance(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return The tolerance value.
    pure function nr_get_grad_tol(this) result(x)
        class(nonlinear_regression), intent(in) :: this
        real(real64) :: x
        x = this%m_solver%get_gradient_tolerance()
    end function

! --------------------
    !> @brief Sets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return The tolerance value.
    subroutine nr_set_grad_tol(this, x)
        class(nonlinear_regression), intent(inout) :: this
        real(real64), intent(in) :: x
        call this%m_solver%set_gradient_tolerance(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return True if the iteration status should be printed; else, false.
    pure function nr_get_print_status(this) result(x)
        class(nonlinear_regression), intent(in) :: this
        logical :: x
        x = this%m_solver%get_print_status()
    end function

! --------------------
    !> @brief Sets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in,out] this The nonlinear_regression object.
    !! @param[in] x True if the iteration status should be printed; else, false.
    subroutine nr_set_print_status(this, x)
        class(nonlinear_regression), intent(inout) :: this
        logical, intent(in) :: x
        call this%m_solver%set_print_status(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the number of stored data points.
    !!
    !! @param[in] this The nonlinear_regression object.
    !!
    !! @return The number of data points.
    pure function nr_get_num_pts(this) result(n)
        class(nonlinear_regression), intent(in) :: this
        integer(int32) :: n
        if (this%m_init) then
            n = size(this%m_x)
        else
            n = 0
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the x component of the requested data point.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @param[in] ind The one-based index of the data point to retrieve.
    !!
    !! @return The x component of the requested data point.
    pure function nr_get_x(this, ind) result(x)
        class(nonlinear_regression), intent(in) :: this
        integer(int32), intent(in) :: ind
        real(real64) :: x
        if (this%m_init) then
            x = this%m_x(ind)
        else
            x = 0.0d0
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the y component of the requested data point.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @param[in] ind The one-based index of the data point to retrieve.
    !!
    !! @return The y component of the requested data point.
    pure function nr_get_y(this, ind) result(y)
        class(nonlinear_regression), intent(in) :: this
        integer(int32), intent(in) :: ind
        real(real64) :: y
        if (this%m_init) then
            y = this%m_y(ind)
        else
            y = 0.0d0
        end if
    end function

! ******************************************************************************
! LINEAR REGRESSION
! ------------------------------------------------------------------------------
    !> @brief Employs a least squares fit to determine the coefficient A in the
    !! linear system: Y = A * X.
    !!
    !! @param[in] x An N-element array containing the independent variable data.
    !! @param[in] y An N-element array containing the dependent variable
    !!  data corresponding to @p x.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are different sizes.
    !!
    !! @return The scalar coefficient A.
    !!
    !! @par Example
    !! The following example illustrates how to fit a set of data using a
    !! linear least-squares approach.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use curvefit_regression
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: n = 20
    !!     real(real64), parameter :: xMin = 0.0d0
    !!     real(real64), parameter :: xMax = 5.0d0
    !!     real(real64), parameter :: slope = 10.0d0
    !!
    !!     ! Local Variables
    !!     real(real64) :: m, x(n), y(n), yfit(n)
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: d1, d2
    !!     type(legend), pointer :: lgnd
    !!     character(len = 60) :: txt
    !!
    !!     ! Initialization
    !!     x = linspace(xMin, xMax, n)
    !!     call random_number(y)
    !!     y = slope * x + 2.0d0 * (y - 0.5d0)
    !!
    !!     ! Fit the data
    !!     m = linear_least_squares(x, y)
    !!     yfit = m * x
    !!
    !!     ! Plot the data
    !!     call plt%initialize()
    !!     call plt%set_font_size(14)
    !!
    !!     lgnd => plt%get_legend()
    !!     call lgnd%set_horizontal_position(LEGEND_LEFT)
    !!
    !!     write(txt, '(AF7.3AF7.3)') "Actual Slope: ", slope, ", Fitted Slope: ", m
    !!     call plt%set_title(trim(txt))
    !!
    !!     call d1%set_name("Data")
    !!     call d1%set_draw_line(.false.)
    !!     call d1%set_draw_markers(.true.)
    !!     call d1%set_marker_style(MARKER_X)
    !!     call d1%set_line_color(CLR_BLACK)
    !!     call d1%set_line_width(2.0)
    !!     call d1%set_marker_scaling(1.5)
    !!     call d1%define_data(x, y)
    !!
    !!     call d2%set_name("Fit")
    !!     call d2%set_line_color(CLR_BLUE)
    !!     call d2%set_line_width(2.0)
    !!     call d2%define_data(x, yfit)
    !!
    !!     call plt%push(d1)
    !!     call plt%push(d2)
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @image html linear_least_squares_scalar_example.png
    function linear_least_squares_1var(x, y, err) result(a)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x, y
        class(errors), intent(inout), optional, target :: err
        real(real64) :: a

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: n, flag
        type(polynomial) :: poly
        real(real64), allocatable, dimension(:) :: ycopy

        ! Initialization
        a = zero
        n = size(x)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(y) /= n) then
          call errmgr%report_error("linear_least_squares_1var", &
              "Incompatible array dimensions.", CF_ARRAY_SIZE_ERROR)
          return
        end if

        ! Local Memory Allocation
        allocate(ycopy(n), stat = flag)
        if (flag /= 0) then
          call errmgr%report_error("linear_least_squares_1var", &
              "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
          return
        end if

        ! Process
        ycopy = y
        call poly%fit_thru_zero(x, ycopy, 1, err = errmgr)
        if (errmgr%has_error_occurred()) return
        a = poly%get(2)
    end function

! ------------------------------------------------------------------------------
    !> @brief Employs a least squares fit to determine the coefficient A in the
    !! linear system: Y = A * X.
    !!
    !! @param[in] x An M-by-P matrix containing the P data points of the
    !!  M independent variables.
    !! @param[in] y An N-by-P matrix containing the P data points of the N
    !!  dependent variables.
    !! @param[in] thrsh An optional threshold value that defines a lower cutoff
    !!  for singular values.  Any singular values falling below this value will
    !!  have their reciprocal replaced with zero.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if any of the matrix dimensions are not
    !!      compatiable.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!
    !! @return An N-by-M matrix relating Y to X such that: Y = A * X.
    !!
    !! @par Remarks
    !! The algorithm to compute the coefficient matrix A is as follows.
    !! @verbatim
    !! First, multiply through by the pseudo-inverse of A (pinv(A))
    !! pinv(A) * Y = pinv(A) * A * X
    !!
    !! Next, post multiply by the pseudo-inverse of Y
    !! pinv(A) * Y * pinv(Y) = pinv(A) * A * X * pinv(Y)
    !!
    !! Note, Y * pinv(Y) = I
    !! Also, let D = pinv(A) * A, and let B = X * pinv(Y)
    !!
    !! Then,
    !! pinv(A) = D * B
    !!
    !! Solving for A
    !! A = pinv(D * B) = pinv(D) * pinv(B)
    !!
    !! Remembering that D = pinv(A) * A
    !! A = pinv(pinv(A) * A) * pinv(B) = A * pinv(A) * pinv(B)
    !!
    !! Noting A * pinv(A) = I
    !! A = pinv(B)
    !! @endverbatim
    !!
    !! @par Example
    !! The following example illustrates how to apply the linear_least_squares
    !! solver to a system of multiple inputs and multiple outputs.  The data
    !! from this example comes from an axial-torsional load cell calibration.
    !!
    !! @par
    !! The applied load data.
    !! Axial Load Data: | Torsional Load Data:
    !! :--------------: | :--------------:
    !! 0 | 0
    !! 3000 | 0
    !! 6000 | 0
    !! 7500 | 0
    !! 9000 | 0
    !! 12000 | 0
    !! 15000 | 0
    !! 7500 | 0
    !! 0 | 0
    !! 0 | 0
    !! -3000 | 0
    !! -6000 | 0
    !! -7500 | 0
    !! -9000 | 0
    !! -12000 | 0
    !! -15000 | 0
    !! -7500 | 0
    !! 0 | 0
    !! 0 | 0
    !! 0 | 67.7908728
    !! 0 | 135.5817456
    !! 0 | 203.3726184
    !! 0 | 271.1634912
    !! 0 | 338.954364
    !! 0 | 203.3726184
    !! 0 | 0
    !! 0 | 0
    !! 0 | -67.7908728
    !! 0 | -135.5817456
    !! 0 | -203.3726184
    !! 0 | -271.1634912
    !! 0 | -338.954364
    !! 0 | -203.3726184
    !! 0 | 0
    !!
    !! @par
    !! The corresponding bridge data.
    !! Bridge 1 Output: | Bridge 2 Output:
    !! :--------------: | :--------------:
    !! 0 | 0
    !! 0.38905 | 0.00122
    !! 0.77816 | 0.00259
    !! 0.97269 | 0.0029
    !! 1.16714 | 0.00314
    !! 1.556 | 0.00338
    !! 1.94484 | 0.00356
    !! 0.9726 | 0.00477
    !! -0.00001 | -0.00001
    !! 0 | 0
    !! -0.388886 | 0.00021
    !! -0.77775 | 0.00051
    !! -0.97215 | 0.00069
    !! -1.16654 | 0.00088
    !! -1.55533 | 0.0013
    !! -1.9441 | 0.00175
    !! -0.97171 | 0.00058
    !! 0.00004 | 0.00003
    !! 0 | 0
    !! -0.00044 | 0.27156
    !! -0.0013 | 0.54329
    !! -0.0024 | 0.81507
    !! -0.00382 | 1.08682
    !! -0.00528 | 1.35881
    !! -0.00257 | 0.81553
    !! 0.00015 | 0.0001
    !! 0 | 0
    !! 0.00144 | -0.27145
    !! 0.00306 | -0.54312
    !! 0.00446 | -0.81493
    !! 0.00567 | -1.0868
    !! 0.00688 | -1.35879
    !! 0.00451 | -0.81548
    !! -0.00002 | 0
    !! @endverbatim
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use curvefit_regression
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     integer(int32), parameter :: npts = 34
    !!     integer(int32), parameter :: nchan = 2
    !!     real(real64) :: loads(npts, nchan), bridge(npts, nchan), &
    !!         gain(nchan, nchan), calib(npts, nchan), errs(npts, nchan), &
    !!         axialFullScale, torqueFullScale
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: d1, d2
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!
    !!     ! Populate the applied loads matrix
    !!     loads = reshape([0.0d0, 3000.0d0, 6000.0d0, 7500.0d0, 9000.0d0, 12000.0d0, &
    !!         15000.0d0, 7500.0d0, 0.0d0, 0.0d0, -3000.0d0, -6000.0d0, -7500.0d0, &
    !!         -9000.0d0, -12000.0d0, -15000.0d0, -7500.0d0, 0.0d0, 0.0d0, 0.0d0, &
    !!         0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
    !!         0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
    !!         0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
    !!         0.0d0, 0.0d0, 0.0d0, 67.7908728d0, 135.5817456d0, 203.3726184d0, &
    !!         271.1634912d0, 338.954364d0, 203.3726184d0, 0.0d0, 0.0d0, &
    !!         -67.7908728d0, -135.5817456d0, -203.3726184d0, -271.1634912d0, &
    !!         -338.954364d0, -203.3726184d0, 0.0d0], [npts, nchan])
    !!
    !!     ! Populate the bridge output (measured) matrix
    !!     bridge = reshape([0.0d0, 0.38905d0, 0.77816d0, 0.97269d0, 1.16714d0, &
    !!         1.556d0, 1.94484d0, 0.9726d0, -0.00001d0, 0.0d0, -0.388886d0, &
    !!         -0.77775d0, -0.97215d0, -1.16654d0, -1.55533d0, -1.9441d0, -0.97171d0, &
    !!         0.00004d0, 0.0d0, -0.00044d0, -0.0013d0, -0.0024d0, -0.00382d0, &
    !!         -0.00528d0, -0.00257d0, 0.00015d0, 0.0d0, 0.00144d0, 0.00306d0, &
    !!         0.00446d0, 0.00567d0, 0.00688d0, 0.00451d0, -0.00002d0, 0.0d0, &
    !!         0.00122d0, 0.00259d0, 0.0029d0, 0.00314d0, 0.00338d0, 0.00356d0, &
    !!         0.00477d0, -0.00001d0, 0.0d0, 0.00021d0, 0.00051d0, 0.00069d0, &
    !!         0.00088d0, 0.0013d0, 0.00175d0, 0.00058d0, 0.00003d0, 0.0d0, &
    !!         0.27156d0, 0.54329d0, 0.81507d0, 1.08682d0, 1.35881d0, 0.81553d0, &
    !!         0.0001d0, 0.0d0, -0.27145d0, -0.54312d0, -0.81493d0, -1.0868d0, &
    !!         -1.35879d0, -0.81548d0, 0.0d0], [npts, nchan])
    !!
    !!     ! Compute the coefficient matrix.  The transpose operations are necessary as
    !!     ! the data was input in column-major format in the order illustrated in the
    !!     ! comments.  Had the data been input differently, the transpose operations
    !!     ! could have been avoided.
    !!     gain = linear_least_squares(transpose(bridge), transpose(loads))
    !!
    !!     ! Apply the matrix.
    !!     calib = matmul(bridge, transpose(gain))
    !!
    !!     ! Compute errors in the fit
    !!     errs = calib - loads
    !!
    !!     ! Express the errors in terms of percent of full scale
    !!     axialFullScale = maxval(loads(:,1))
    !!     torqueFullScale = maxval(loads(:,2))
    !!     errs(:,1) = 1.0d2 * errs(:,1) / axialFullScale
    !!     errs(:,2) = 1.0d2 * errs(:,2) / torqueFullScale
    !!
    !!     ! Plot the error data vs. index of the applied load
    !!     call plt%initialize()
    !!     call plt%set_font_size(14)
    !!
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!
    !!     call xAxis%set_title("Index")
    !!     call yAxis%set_title("Error [% FS]")
    !!
    !!     call d1%set_name("Axial")
    !!     call d1%set_line_color(CLR_BLUE)
    !!     call d1%set_line_width(2.0)
    !!     call d1%set_draw_markers(.true.)
    !!     call d1%set_marker_style(MARKER_X)
    !!     call d1%define_data(errs(:,1))
    !!
    !!     call d2%set_name("Torsional")
    !!     call d2%set_line_color(CLR_RED)
    !!     call d2%set_line_width(2.0)
    !!     call d2%set_draw_markers(.true.)
    !!     call d2%set_marker_style(MARKER_EMPTY_CIRCLE)
    !!     call d2%define_data(errs(:,2))
    !!
    !!     call plt%push(d1)
    !!     call plt%push(d2)
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! @image html linear_least_squares_mimo_example.png
    function linear_least_squares_nvar(x, y, thrsh, err) result(a)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: x, y
        real(real64), intent(in), optional :: thrsh
        class(errors), intent(inout), optional, target :: err
        real(real64), dimension(size(y,1), size(x,1)) :: a

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: m, n, p, flag
        real(real64), allocatable, dimension(:,:) :: b, yinv, ycopy

        ! Initialization
        m = size(x, 1)
        p = size(x, 2)
        n = size(y, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        ! X is M-by-P
        ! Y is N-by-P
        if (size(y, 2) /= p) then
            call errmgr%report_error("linear_least_squares_nvar", &
                "Incompatible array dimensions.", CF_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(ycopy(n, p), stat = flag)
        if (flag == 0) allocate(yinv(p, n), stat = flag)
        if (flag == 0) allocate(b(m, n), stat = flag)

        if (flag /= 0) then
            call errmgr%report_error("linear_least_squares_nvar", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Compute the pseudo-inverse of Y
        ycopy = y ! Need a copy of Y as mtx_pinverse modifies the matrix
        call mtx_pinverse(ycopy, yinv, tol = thrsh, err = errmgr)
        if (errmgr%has_error_occurred()) return

        ! Compute B = X * pinv(Y)
        call mtx_mult(.false., .false., one, x, yinv, zero, b)

        ! Compute A = pinv(B)
        call mtx_pinverse(b, a, tol = thrsh, err = errmgr)

        ! ----------------------------------------------------------------------
        ! Previous code, replaced with the above code by JAC on 23-April, 2018
        ! ! Local Memory Allocation
        ! allocate(xinv(p, m), stat = flag)
        ! if (flag /= 0) then
        !     call errmgr%report_error("linear_least_squares_nvar", &
        !         "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
        !     return
        ! end if
        !
        ! ! Compute the pseudo-inverse of X
        ! call mtx_pinverse(x, xinv, tol = thrsh, err = errmgr)
        ! if (errmgr%has_error_occurred()) return
        !
        ! ! Compute A = Y * pinv(X)
        ! call mtx_mult(.false., .false., one, y, xinv, zero, a)
        ! ----------------------------------------------------------------------
    end function

! ------------------------------------------------------------------------------
end module
