! curvefit_regression.f90

!> @brief \b curvefit_regression
!!
!! @par Purpose
!! To provide routines for perforing regression operations, and other data
!! smoothing operations on sets of numerical data.
module curvefit_regression
    use curvefit_core
    use linalg_sorting, only : sort
    use ferror, only : errors
    use nonlin_types, only : vecfcn_helper, iteration_behavior
    use nonlin_least_squares, only : least_squares_solver
    use curvefit_statistics, only : mean
    use linalg_solve, only : mtx_pinverse
    use linalg_core, only : mtx_mult
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
    type lowess_smoothing
        private
        !> N-element array of x data points - sorted into ascending order.
        real(dp), allocatable, dimension(:) :: m_x
        !> N-element array of y data points.
        real(dp), allocatable, dimension(:) :: m_y
        !> N-element array containing the robustness weights for each data 
        !! point.
        real(dp), allocatable, dimension(:) :: m_weights
        !> N-element array containing the residuals (Y - YS)
        real(dp), allocatable, dimension(:) :: m_residuals
        !> Scaling parameter used to define the nature of the linear 
        !! interpolations used by the algorithm.
        real(dp) :: m_delta
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
    end type

! ------------------------------------------------------------------------------
    !> @brief A type for supporting nonlinear regression calculations.
    type, extends(vecfcn_helper) :: nonlinear_regression
        private
        !> A pointer to the routine containing the function of interest.
        procedure(reg_fcn), pointer, nopass :: m_rfcn => null()
        !> The x data points.
        real(dp), allocatable, dimension(:) :: m_x
        !> The y data points.
        real(dp), allocatable, dimension(:) :: m_y
        !> The number of coefficients in the function of interest
        integer(i32) :: m_ncoeff = 0
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
        real(dp), intent(inout), dimension(:) :: x
        integer(i32), intent(in) :: npts
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(i32) :: i, n, flag
        real(dp), allocatable, dimension(:) :: buffer

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
        real(dp), intent(in), dimension(:) :: x, y, rw ! N ELEMENT
        real(dp), intent(in) :: xs
        real(dp), intent(out) :: ys
        integer(i32), intent(in) :: nleft, nright
        real(dp), intent(out), dimension(:) :: w ! N ELEMENT
        logical, intent(in) :: userw
        logical, intent(out) :: ok

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: p001 = 1.0d-3
        real(dp), parameter :: p999 = 0.999d0

        ! Local Variables
        integer(i32) :: j, n, nrt
        real(dp) :: range, h, h9, h1, a, b, c, r

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
        real(dp), intent(in), dimension(:) :: x, y
        real(dp), intent(in) :: f
        integer(i32), intent(in) :: nsteps
        real(dp), intent(in) :: delta
        real(dp), intent(out), dimension(:) :: ys, rw, res

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: three = 3.0d0
        real(dp), parameter :: p001 = 1.0d-3
        real(dp), parameter :: p999 = 0.999d0

        ! Local Variables
        logical :: ok
        integer(i32) :: iter, i, j, n, nleft, nright, ns, last, m1, m2
        real(dp) :: d1, d2, denom, alpha, cut, eps, cmad, c1, c9, r

        ! Initialization
        n = size(x)
        ns = max(min(int(f * real(n, dp), i32), n), 2)
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
        real(dp), intent(in), dimension(:) :: x, y
        logical, intent(in), optional :: srt
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: i, n, flag
        integer(i32), allocatable, dimension(:) :: indices
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
    !! @param[in,out] this THe lowess_smoothing object.
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
        real(dp), intent(in) :: f
        class(errors), intent(inout), optional, target :: err
        real(dp), allocatable, dimension(:) :: ys

        ! Local Variables
        integer(i32) :: n, flag
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
        integer(i32) :: n
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
        integer(i32), intent(in) :: ind
        real(dp) :: x
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
        integer(i32), intent(in) :: ind
        real(dp) :: y
        if (this%m_init) then
            y = this%m_y(ind)
        else
            y = 0.0d0
        end if
    end function

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
        real(dp), intent(in), dimension(:) :: x, y
        procedure(reg_fcn), pointer, intent(in) :: fcn
        integer(i32), intent(in) :: ncoeff
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(i32) :: i, n, flag

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
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: f

        ! Local Variables
        integer(i32) :: i, n

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
        integer(i32) :: n
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
        integer(i32) :: n
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
        real(dp), intent(inout), dimension(:) :: c
        real(dp), intent(out), dimension(:), target, optional :: res
        type(iteration_behavior), intent(out), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(i32) :: n, flag
        real(dp), allocatable, target, dimension(:) :: f
        real(dp), pointer, dimension(:) :: fptr

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
        integer(i32) :: n
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
        integer(i32), intent(in) :: n
        call this%m_solver%set_max_fcn_evals(n)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on function value tolerance.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return The tolerance value.
    pure function nr_get_fcn_tol(this) result(x)
        class(nonlinear_regression), intent(in) :: this
        real(dp) :: x
        x = this%m_solver%get_fcn_tolerance()
    end function

! --------------------
    !> @brief Sets the convergence on function value tolerance.
    !!
    !! @param[in,out] this The nonlinear_regression object.
    !! @param[in] x The tolerance value.
    subroutine nr_set_fcn_tol(this, x)
        class(nonlinear_regression), intent(inout) :: this
        real(dp), intent(in) :: x
        call this%m_solver%set_fcn_tolerance(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on change in variable tolerance.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return The tolerance value.
    pure function nr_get_var_tol(this) result(x)
        class(nonlinear_regression), intent(in) :: this
        real(dp) :: x
        x = this%m_solver%get_var_tolerance()
    end function

! --------------------
    !> @brief Sets the convergence on change in variable tolerance.
    !!
    !! @param[in,out] this The nonlinear_regression object.
    !! @param[in] x The tolerance value.
    subroutine nr_set_var_tol(this, x)
        class(nonlinear_regression), intent(inout) :: this
        real(dp), intent(in) :: x
        call this%m_solver%set_var_tolerance(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return The tolerance value.
    pure function nr_get_grad_tol(this) result(x)
        class(nonlinear_regression), intent(in) :: this
        real(dp) :: x
        x = this%m_solver%get_gradient_tolerance()
    end function

! --------------------
    !> @brief Sets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The nonlinear_regression object.
    !! @return The tolerance value.
    subroutine nr_set_grad_tol(this, x)
        class(nonlinear_regression), intent(inout) :: this
        real(dp), intent(in) :: x
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
        integer(i32) :: n
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
        integer(i32), intent(in) :: ind
        real(dp) :: x
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
        integer(i32), intent(in) :: ind
        real(dp) :: y
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
    !! @param[in,out] y An N-element array containing the dependent variable
    !!  data corresponding to @p x.  On output, the contents of this array are
    !!  overwritten as it is used for storage purposes by the algorithm.
    !! @param[out] err
    !!
    !! @return The scalar coefficient A.
    function linear_least_squares_1var(x, y, err) result(a)
        ! Arguments
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(inout), dimension(:) :: y
        class(errors), intent(inout), optional, target :: err
        real(dp) :: a

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(i32) :: n
        type(polynomial) :: poly

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
        end if

        ! Process
        call poly%fit_thru_zero(x, y, 1, err = errmgr)
        if (errmgr%has_error_occurred()) return
        a = poly%get(2)
    end function

! ------------------------------------------------------------------------------
    !> @brief Employs a least squares fit to determine the coefficient A in the
    !! linear system: Y = A * X.
    !!
    !! @param[in,out] x An M-by-P matrix containing the P data points of the
    !!  M independent variables.
    !! @param[in] y An N-by-P matrix containing the P data points of the N
    !!  dependent variables.
    !! @param[in] thrsh An optional threshold value that defines a lower cutoff
    !!  for singular values.  Any singular values falling below this value will
    !!  have their reciprocal replaced with zero.
    !! @param[out] err
    !!
    !! @return An N-by-M matrix relating Y to X such that: Y = A * X.
    !!
    !! @par Remarks
    !! The algorithm attempts to compute the coefficient matrix A as follows.
    !! Y * X**T = A * X * X**T
    !! Y * X**T * INV(X * X**T) = A
    !! This does require that X * X**T does not result in a singular matrix.  To
    !! handle the situation where X * X**T is singular, the Moore-Penrose
    !! pseudo-inverse, computed by means of singular value decomposition, is
    !! utilized to still arrive at a solution that, at minimum, has a minimum
    !! Euclidean norm of its residual.
    !! Let: PINV(X) = X**T * INV(X * X**T),
    !! Then: A = Y * PINV(X)
    function linear_least_squares_nvar(x, y, thrsh, err) result(a)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: x
        real(dp), intent(in), dimension(:,:) :: y
        real(dp), intent(in), optional :: thrsh
        class(errors), intent(inout), optional, target :: err
        real(dp), dimension(size(y,1), size(x,1)) :: a

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(i32) :: m, n, p, flag
        real(dp), allocatable, dimension(:,:) :: xinv

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
        allocate(xinv(p, m), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("linear_least_squares_nvar", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Compute the pseudo-inverse of X
        call mtx_pinverse(x, xinv, tol = thrsh, err = errmgr)
        if (errmgr%has_error_occurred()) return

        ! Compute A = Y * pinv(X)
        call mtx_mult(.false., .false., one, y, xinv, zero, a)
    end function

! ------------------------------------------------------------------------------
end module
