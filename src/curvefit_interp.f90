! curvefit_interp.f90

module curvefit_interp
    use curvefit_core
    use ferror, only : errors
    implicit none
    private
    public :: SPLINE_QUADRATIC_OVER_INTERVAL
    public :: SPLINE_KNOWN_FIRST_DERIVATIVE
    public :: SPLINE_KNOWN_SECOND_DERIVATIVE
    public :: SPLINE_CONTINUOUS_THIRD_DERIVATIVE
    public :: interp_manager
    public :: linear_interp
    public :: polynomial_interp
    public :: spline_interp

! ******************************************************************************
! CONSTANTS
! ------------------------------------------------------------------------------
    !> Indicates that the spline is quadratic over the interval under
    !! consideration (beginning or ending interval).  This is equivalent to
    !! allowing a "free" boundary condition at either the initial or final
    !! point.
    integer(i32), parameter :: SPLINE_QUADRATIC_OVER_INTERVAL = 1000
    !> Indicates a known first derivative at either the beginning or ending 
    !! point.
    integer(i32), parameter :: SPLINE_KNOWN_FIRST_DERIVATIVE = 1001
    !> Indicates a known second derivative at either the beginning or ending 
    !! point.
    integer(i32), parameter :: SPLINE_KNOWN_SECOND_DERIVATIVE = 1002
    !> Indicates a continuous third derivative at either the beginning or ending 
    !! point.
    integer(i32), parameter :: SPLINE_CONTINUOUS_THIRD_DERIVATIVE = 1003

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Describes an abstract base class allowing for interpolation of X-Y
    !! type data sets.
    !!
    !! @par Notes
    !! This interpolation object is based upon the interpolation scheme utilized
    !! by the Numerical Recipes in C++ text.
    type, abstract :: interp_manager
    private
        integer(i32) :: m_order
        integer(i32) :: m_savedIndex
        integer(i32) :: m_indexCheck
        logical :: m_correlated
        real(dp), allocatable, dimension(:) :: m_x
        real(dp), allocatable, dimension(:) :: m_y
    contains
        !> @brief Initializes the interp_manager instance.
        procedure, public :: initialize => im_init
        !> @brief Attempts to locate the index in the array providing a lower 
        !! bounds to the specified interpolation point.
        procedure, non_overridable, public :: locate => im_locate
        !> @brief Attempts to locate the index in the array providing a lower 
        !! bounds to the specified interpolation point.
        procedure, non_overridable, public :: hunt => im_hunt
        !> @brief Interpolates to obtain the function value at the specified
        !!  independent variable.
        generic, public :: interpolate => im_perform, im_perform_array
        !> @brief Performs the actual interpolation.
        procedure(interp_xy), deferred :: raw_interp

        procedure, non_overridable :: im_perform
        procedure, non_overridable :: im_perform_array
    end type

! ------------------------------------------------------------------------------
    !> @brief Extends the interp_manager class allowing for linear, piecewise
    !! interpolation of a data set.
    type, extends(interp_manager) :: linear_interp
    contains
        !> @brief Performs the actual interpolation.
        procedure :: raw_interp => li_raw_interp
    end type

! ------------------------------------------------------------------------------
    !> @brief Extends the interp_manager class allowing for polynomial 
    !! interpolation of a data set.
    type, extends(interp_manager) :: polynomial_interp
    private
        real(dp), allocatable, dimension(:) :: m_c
        real(dp), allocatable, dimension(:) :: m_d
        real(dp) :: m_dy
    contains
        !> @brief Initializes the polynomial_interp instance.
        procedure, public :: initialize => pi_init
        !> @brief Performs the actual interpolation.
        procedure :: raw_interp => pi_raw_interp
    end type

! ------------------------------------------------------------------------------
    !> @brief Extends the interp_manager class allowing for cubic spline 
    !! interpolation of a data set.
    type, extends(interp_manager) :: spline_interp
    private
        real(dp), allocatable, dimension(:) :: m_ypp
        real(dp), allocatable, dimension(:) :: m_a1
        real(dp), allocatable, dimension(:) :: m_a2
        real(dp), allocatable, dimension(:) :: m_a3
        real(dp), allocatable, dimension(:) :: m_a4
        real(dp), allocatable, dimension(:) :: m_a5
        real(dp), allocatable, dimension(:) :: m_b
    contains
        !> @brief Performs the actual interpolation.
        procedure :: raw_interp => si_raw_interp
        !> @brief Computes the second derivative terms for the cubic-spline 
        !! model.
        procedure :: compute_diff2 => si_second_deriv
        !> @brief Initializes the spline_interp instance.
        procedure, public :: initialize => si_init_1
        !> @brief Initializes the spline_interp instance while allowing 
        !! definition of boundary conditions.
        procedure, public :: initialize_spline => si_init_2
    end type


! ******************************************************************************
! ABSTRACT INTERFACES
! ------------------------------------------------------------------------------
interface
    !> @brief Defines the signature of a method used to interpolate a single
    !!  value in an X-Y data set.
    !!
    !! @param[in,out] this The interp_manager based instance.
    !! @param[in] jlo The array index below which @p pt is found in x.
    !! @param[in] pt The independent variable value to interpolate.
    !!
    !! @return The interpolated value.
    function interp_xy(this, jlo, pt) result(yy)
        use curvefit_core, only : dp, i32
        import interp_manager
        class(interp_manager), intent(inout) :: this
        integer(i32), intent(in) :: jlo
        real(dp), intent(in) :: pt
        real(dp) :: yy
    end function
end interface


contains
! ******************************************************************************
! INTERP_MANAGER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Initializes the specified interp_manager instance.
    !!
    !! @param[in,out] this The interp_manager instance.
    !! @param[in] x An N-element array containing the independent variable data.
    !! @param[in] y An N-element array containing the dependent variable data.
    !! @param[in] order The order of the interpolating polynomial.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    subroutine im_init(this, x, y, order, err)
        ! Arguments
        class(interp_manager), intent(inout) :: this
        real(dp), intent(in), dimension(:) :: x, y
        integer(i32), intent(in) :: order
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: i, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        this%m_order = order
        this%m_savedIndex = 1
        this%m_indexCheck = 1
        n = size(x)
        if (size(y) /= n) then
            ! ERROR
            write(errmsg, '(AI0AI0A)') &
                "Expected the dependent variable array to be of length ", &
                size(x), ", but found an array of length ", size(y), "."
            call errmgr%report_error("im_init", trim(errmsg), &
                CF_ARRAY_SIZE_ERROR)
            return
        end if

        if (allocated(this%m_x)) deallocate(this%m_x)
        if (allocated(this%m_y)) deallocate(this%m_y)

        allocate(this%m_x(n), stat = flag)
        if (flag == 0) allocate(this%m_y(n), stat = flag)
        if (flag /= 0) then
            ! ERROR
            call errmgr%report_error("im_init", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Copy the data
        do i = 1, n
            this%m_x(i) = x(i)
            this%m_y(i) = y(i)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Attempts to locate the index in the array providing a lower bounds
    !!  to the specified interpolation point.
    !!
    !! @param[in,out] this The interp_manager instance.
    !! @param[in] pt The interpolation point.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
    !!
    !! @return The array index below @p pt.
    function im_locate(this, pt, err) result(j)
        ! Arguments
        class(interp_manager), intent(inout) :: this
        real(dp), intent(in) :: pt
        class(errors), intent(inout), optional, target :: err
        integer :: j

        ! Local Variables
        integer(i32) :: n, m, jhi, jmid, jlo
        logical :: ascnd
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        j = 0
        n = size(this%m_x)
        m = this%m_order + 1
        ascnd = this%m_x(n) >= this%m_x(1)
        jlo = 1
        jhi = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Ensure data has been defined
        if (.not.allocated(this%m_x) .or. .not.allocated(this%m_y)) then
            call errmgr%report_error("im_locate", "No data has been defined.", &
                CF_NO_DATA_DEFINED_ERROR)
            return
        end if

        ! Process
        do while (jhi - jlo > 1)
            jmid = (jhi + jlo) / 2
            if (pt >= this%m_x(jmid) .eqv. ascnd) then
                jlo = jmid
            else
                jhi = jmid
            end if
        end do

        ! Check to see if we should use a more efficient search approach next
        ! time
        this%m_correlated = abs(jlo - this%m_savedIndex) <= this%m_indexCheck
        this%m_savedIndex = jlo

        ! Output
        j = max(1, min(n + 1 - m, jlo - (m - 1) / 2))
    end function

! ------------------------------------------------------------------------------
    !> @brief Attempts to locate the index in the array providing a lower bounds
    !!  to the specified interpolation point.  This method is typically more
    !!  efficient than locate when the current index does not stray
    !!  too far from the previous.
    !!
    !! @param[in,out] this The interp_manager instance.
    !! @param[in] pt The interpolation point.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
    !!
    !! @return The array index below @p pt.
    function im_hunt(this, pt, err) result(j)
        ! Arguments
        class(interp_manager), intent(inout) :: this
        real(dp), intent(in) :: pt
        class(errors), intent(inout), optional, target :: err
        integer(i32) :: j

        ! Local Variables
        integer(i32) :: jlo, jmid, jhi, inc, n, m
        logical :: ascnd
        class(errors), pointer :: errmgr
        type(errors), target :: deferr


        ! Initialization
        j = 0
        n = size(this%m_x)
        m = this%m_order + 1
        jlo = this%m_savedIndex
        inc = 1
        ascnd = this%m_x(n) > this%m_x(1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Ensure data has been defined
        if (.not.allocated(this%m_x) .or. .not.allocated(this%m_y)) then
            call errmgr%report_error("im_hunt", "No data has been defined.", &
                CF_NO_DATA_DEFINED_ERROR)
            return
        end if

        ! Process
        if (jlo < 1 .or. jlo > n) then
            jlo = 1
            jhi = n
        else
            if (pt >= this%m_x(jlo) .eqv. ascnd) then
                do
                    jhi = jlo + inc
                    if (jhi >= n) then
                        jhi = n
                        exit
                    else if (pt < this%m_x(jhi) .eqv. ascnd) then
                        exit
                    else
                        jlo = jhi
                        inc = inc + inc
                    end if
                end do
            else
                jhi = jlo
                do
                    jlo = jlo - inc
                    if (jlo <= 1) then
                        jlo = 1
                        exit
                    else if (pt >= this%m_x(jlo) .eqv. ascnd) then
                        exit
                    else
                        jhi = jlo
                        inc = inc + inc
                    end if
                end do
            end if
        end if

        ! The hunt is done, so begin the final bisection phase
        do while (jhi - jlo > 1)
            jmid = (jhi + jlo) / 2
            if (pt >= this%m_x(jmid) .eqv. ascnd) then
                jlo = jmid
            else
                jhi = jmid
            end if
        end do

        ! Check to see if we should hunt or locate the next time around
        this%m_correlated = abs(jlo - this%m_savedIndex) <= this%m_indexCheck
        this%m_savedIndex = jlo

        ! Output
        j = max(1, min(n + 1 - m, jlo - (m - 1) / 2))
    end function

! ------------------------------------------------------------------------------
    !> @brief Interpolates to obtain the function value at the specified
    !!  independent variable.
    !!
    !! @param[in,out] this The interp_manager instance.
    !! @param[in] pt The independent variable value to interpolate.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
    !!
    !! @return The interpolated value.
    function im_perform(this, pt, err) result(yy)
        ! Arguments
        class(interp_manager), intent(inout) :: this
        real(dp), intent(in) :: pt
        class(errors), intent(inout), optional, target :: err
        real(dp) :: yy

        ! Local Variables
        integer(i32) :: jlo

        ! Process
        if (this%m_correlated) then
            jlo = this%hunt(pt, err)
        else
            jlo = this%locate(pt, err)
        end if
        yy = this%raw_interp(jlo, pt)
    end function

! ------------------------------------------------------------------------------
    !> @brief Interpolates to obtain the function value at the specified
    !!  independent variables.
    !!
    !! @param[in,out] this The interp_manager instance.
    !! @param[in] pt An M-element array containing the independent variable 
    !!  values to interpolate.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_NO_DATA_DEFINED_ERROR: Occurs if no data has yet been defined.
    !!
    !! @return An M-element array containing the interpolated values.
    function im_perform_array(this, pts, err) result(yy)
        ! Arguments
        class(interp_manager), intent(inout) :: this
        real(dp), intent(in), dimension(:) :: pts
        class(errors), intent(inout), optional, target :: err
        real(dp), dimension(size(pts)) :: yy

        ! Local Variables
        integer(i32) :: i, jlo

        ! Process
        do i = 1, size(pts)
            if (this%m_correlated) then
                jlo = this%hunt(pts(i), err)
            else
                jlo = this%locate(pts(i), err)
            end if
            yy(i) = this%raw_interp(jlo, pts(i))
        end do
    end function

! ******************************************************************************
! LINEAR_INTERP MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Performs the actual linear interpolation.
    !!
    !! @param[in,out] this The linear_interp_mgr instance.
    !! @param[in] jlo The array index below which @p pt is found in x.
    !! @param[in] pt The independent variable value to interpolate.
    !!
    !! @return The interpolated value.
    function li_raw_interp(this, jlo, pt) result(yy)
        ! Arguments
        class(linear_interp), intent(inout) :: this
        integer(i32), intent(in) :: jlo
        real(dp), intent(in) :: pt
        real(dp) :: yy

        ! Process
        if (this%m_x(jlo) == this%m_x(jlo+1)) then
            yy = this%m_y(jlo)
        else
            yy = this%m_y(jlo) + ((pt - this%m_x(jlo)) / &
                (this%m_x(jlo+1) - this%m_x(jlo))) * &
                (this%m_y(jlo+1) - this%m_y(jlo))
        end if
    end function

! ******************************************************************************
! POLYNOMIAL_INTERP MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Initializes the specified polynomial_interp instance.
    !!
    !! @param[in,out] this The polynomial_interp instance.
    !! @param[in] x An N-element array containing the independent variable data.
    !! @param[in] y An N-element array containing the dependent variable data.
    !! @param[in] order The order of the interpolating polynomial.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
    subroutine pi_init(this, x, y, order, err)
        ! Arguments
        class(polynomial_interp), intent(inout) :: this
        real(dp), intent(in), dimension(:) :: x, y
        integer(i32), intent(in) :: order
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        if (order < 1) then
            call errmgr%report_error("pi_init", &
                "A polynomial order greater than or equal to 1 must " // &
                "be specified.", CF_INVALID_INPUT_ERROR)
            return
        end if

        ! Memory Allocation
        call im_init(this, x, y, order, err)
        if (allocated(this%m_c)) deallocate(this%m_c)
        if (allocated(this%m_d)) deallocate(this%m_d)
        m = order + 1
        allocate(this%m_c(m), stat = flag)
        if (flag == 0) allocate(this%m_d(m), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("pi_init", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Performs the actual interpolation.
    !!
    !! @param[in,out] this The polynomial_interp instance.
    !! @param[in] jlo The array index below which @p pt is found in x.
    !! @param[in] pt The independent variable value to interpolate.
    !!
    !! @return The interpolated value.
    function pi_raw_interp(this, jlo, pt) result(yy)
        ! Arguments
        class(polynomial_interp), intent(inout) :: this
        integer(i32), intent(in) :: jlo
        real(dp), intent(in) :: pt
        real(dp) :: yy

        ! Local Variables
        integer(i32) :: i, ind, m, ns, mm, jl
        real(dp) :: den, dif, dift, ho, hp, w

        ! Initialization
        mm = this%m_order + 1
        ns = 1
        jl = jlo - 1
        dif = abs(pt - this%m_x(jl + 1))

        ! Find the index NS of the closest table entry, and then initialize
        ! the C and D arrays.
        do i = 1, mm
            ind = jl + i
            dift = abs(pt - this%m_x(ind))
            if (dift < dif) then
                ns = i
                dif = dift
            end if
            this%m_c(i) = this%m_y(ind)
            this%m_d(i) = this%m_y(ind)
        end do

        ! Define the initial approximation to the interpolated point
        yy = this%m_y(jl + ns)
        ns = ns - 1

        ! Build the tables, and define the interpolated point
        do m = 1, mm-1
            do i = 1, mm - m
                ind = jl + i
                ho = this%m_x(ind) - pt
                hp = this%m_x(ind+m) - pt
                w = this%m_c(i+1) - this%m_d(i)
                den = ho - hp
                den = w / den
                this%m_d(i) = hp * den
                this%m_c(i) = ho * den
            end do
            if (2 * ns < mm - m) then
                this%m_dy = this%m_c(ns + 1)
            else
                this%m_dy = this%m_d(ns)
                ns = ns - 1
            end if
            yy = yy + this%m_dy
        end do
    end function

! ******************************************************************************
! SPLINE_INTERP MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Solves a pentadiagonal system of linear equations.  A
    !!  pentadiagonal matrix is all zeros with the exception of the diagonal,
    !!  and the two immediate sub and super-diagonals.  The entries of row I
    !!  are stored as follows:
    !!      A(I,I-2) -> A1(I)
    !!      A(I,I-1) -> A2(I)
    !!      A(I,I) -> A3(I)
    !!      A(I,I+1) -> A4(I)
    !!      A(I,I+2) -> A5(I)
    !!
    !! @param[in] a1 An N-element array as defined above.
    !! @param[in,out] a2 An N-element array as defined above.  This array is
    !!  overwritten by this routine during the solution process.
    !! @param[in,out] a3 An N-element array as defined above.  This array is
    !!  overwritten by this routine during the solution process.
    !! @param[in,out] a4 An N-element array as defined above.  This array is
    !!  overwritten by this routine during the solution process.
    !! @param[in] a5 An N-element array as defined above.
    !! @param[in,out] b An N-element array containing the right-hand-side.  This
    !!  array is overwritten by this routine during the solution process.
    !! @param[out] x An N-element array that, on output, contains the solution
    !!  to the linear system.
    !!
    !! - [Spline Library](http://people.sc.fsu.edu/~jburkardt/f77_src/spline/spline.html)
    subroutine penta_solve(a1, a2, a3, a4, a5, b, x)
        ! Arguments
        real(dp), intent(in), dimension(:) :: a1, a5
        real(dp), intent(inout), dimension(:) :: a2, a3, a4, b
        real(dp), intent(out), dimension(:) :: x

        ! Local Variables
        integer(i32) :: i, n
        real(dp) :: xmult

        ! Initialization
        n = size(a1)

        ! Process
        do i = 2, n - 1
            xmult = a2(i) / a3(i - 1)
            a3(i) = a3(i) - xmult * a4(i - 1)
            a4(i) = a4(i) - xmult * a5(i - 1)
            b(i) = b(i) - xmult * b(i - 1)
            xmult = a1(i + 1) - xmult * a4(i - 1)
            a2(i + 1) = a2(i + 1) - xmult * a4(i - 1)
            a3(i + 1) = a3(i + 1) - xmult * a5(i - 1)
            b(i + 1) = b(i + 1) - xmult * b(i - 1)
        end do

        xmult = a2(n) / a3(n - 1)
        a3(n) = a3(n) - xmult * a4(n - 1)
        x(n) = (b(n) - xmult * b(n - 1)) / a3(n)
        x(n - 1) = (b(n - 1) - a4(n - 1) * x(n)) / a3(n - 1)
        do i = n - 2, 1, -1
            x(i) = (b(i) - a4(i) * x(i + 1) - a5(i) * x(i + 2)) / a3(i)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Performs the actual interpolation.
    !!
    !! @param[in,out] this The spline_interp instance.
    !! @param[in] jlo The array index below which @p pt is found in x.
    !! @param[in] pt The independent variable value to interpolate.
    !!
    !! @return The interpolated value.
    function si_raw_interp(this, jlo, pt) result(yy)
        ! Arguments
        class(spline_interp), intent(inout) :: this
        integer(i32), intent(in) :: jlo
        real(dp), intent(in) :: pt
        real(dp) :: yy

        ! Parameters
        real(dp), parameter :: half = 0.5d0
        real(dp), parameter :: three = 3.0d0
        real(dp), parameter :: six = 6.0d0

        ! Local Variables
        integer(i32) :: right
        real(dp) :: dt, h

        ! Initialization
        right = jlo + 1
        dt = pt - this%m_x(jlo)
        h = this%m_x(right) - this%m_x(jlo)

        ! Process
        yy = this%m_y(jlo) + dt * ((this%m_y(right) - this%m_y(jlo)) / h - &
            (this%m_ypp(right) / six + this%m_ypp(jlo) / three) * h + &
            dt * (half * this%m_ypp(jlo) + &
            dt * ((this%m_ypp(right) - this%m_ypp(jlo)) / (six * h))))
    end function

! ------------------------------------------------------------------------------
    !> @brief Computes the second derivative terms for the cubic-spline model.
    !!
    !! @param[in,out] this The spline_interp_mgr instance.
    !! @param[in] ibcbeg Defines the nature of the boundary condition at the
    !!  beginning of the spline.
    !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
    !!      initial interval.
    !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
    !!      initial point is provided in @p ybcbeg.
    !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
    !!      initial point is provided in @p ybcbeg.
    !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
    !!      at x(2).
    !! @param[in] ybcbeg If needed, the value of the initial point boundary
    !!  condition.
    !! @param[in] ibcend Defines the nature of the boundary condition at the
    !!  end of the spline.
    !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
    !!      final interval.
    !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
    !!      initial point is provided in @p ybcend.
    !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
    !!      initial point is provided in @p ybcend.
    !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
    !!      at x(n-1).
    !! @param[in] ybcend If needed, the value of the final point boundary
    !!  condition.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!
    !! @par Remarks
    !! This code is a slight modification of the SPLINE_CUBIC_SET routine from
    !! the [SPLINE]
    !! (http://people.sc.fsu.edu/~jburkardt/f77_src/spline/spline.html) library.
    subroutine si_second_deriv(this, ibcbeg, ybcbeg, ibcend, ybcend, err)
        ! Arguments
        class(spline_interp), intent(inout) :: this
        integer(i32), intent(in) :: ibcbeg, ibcend
        real(dp), intent(in) :: ybcbeg, ybcend
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: three = 3.0d0
        real(dp), parameter :: six = 6.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(i32) :: i, n, flag

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        n = size(this%m_x)

        ! Allocate Memory
        if (allocated(this%m_ypp)) deallocate(this%m_ypp)
        if (allocated(this%m_a1)) deallocate(this%m_a1)
        if (allocated(this%m_a2)) deallocate(this%m_a2)
        if (allocated(this%m_a3)) deallocate(this%m_a3)
        if (allocated(this%m_a4)) deallocate(this%m_a4)
        if (allocated(this%m_a5)) deallocate(this%m_a5)
        if (allocated(this%m_b)) deallocate(this%m_b)
        allocate(this%m_ypp(n), stat = flag)
        if (flag == 0) allocate(this%m_a1(n), stat = flag)
        if (flag == 0) allocate(this%m_a2(n), stat = flag)
        if (flag == 0) allocate(this%m_a3(n), stat = flag)
        if (flag == 0) allocate(this%m_a4(n), stat = flag)
        if (flag == 0) allocate(this%m_a5(n), stat = flag)
        if (flag == 0) allocate(this%m_b(n), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("si_second_deriv", &
                "Insufficient memory available.", CF_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Zero out the matrix
        do i = 1, n
            this%m_a1(i) = zero
            this%m_a2(i) = zero
            this%m_a3(i) = zero
            this%m_a4(i) = zero
            this%m_a5(i) = zero
        end do

        ! Set the first equation
        select case (ibcbeg)
        case (SPLINE_QUADRATIC_OVER_INTERVAL)
            this%m_b(1) = zero
            this%m_a3(1) = one
            this%m_a4(1) = one
        case (SPLINE_KNOWN_FIRST_DERIVATIVE)
            this%m_b(1) = (this%m_y(2) - this%m_y(1)) / &
                (this%m_x(2) - this%m_x(1)) - ybcbeg
            this%m_a3(1) = (this%m_x(2) - this%m_x(1)) / three
            this%m_a4(1) = (this%m_x(2) - this%m_x(1)) / six
        case (SPLINE_KNOWN_SECOND_DERIVATIVE)
            this%m_b(1) = ybcbeg
            this%m_a3(1) = one
            this%m_a4(1) = zero
        case (SPLINE_CONTINUOUS_THIRD_DERIVATIVE)
            this%m_b(1) = zero
            this%m_a3(1) = this%m_x(2) - this%m_x(3)
            this%m_a4(1) = this%m_x(3) - this%m_x(1)
            this%m_a5(1) = this%m_x(1) - this%m_x(2)
        case default
            this%m_b(1) = zero
            this%m_a3(1) = one
            this%m_a4(1) = one
        end select

        ! Set the intermediate equations
        do i = 2, n - 1
            this%m_b(i) = (this%m_y(i+1) - this%m_y(i)) / &
                (this%m_x(i+1) - this%m_x(i)) - &
                (this%m_y(i) - this%m_y(i-1)) / (this%m_x(i) - this%m_x(i-1))
            this%m_a2(i) = (this%m_x(i+1) - this%m_x(i)) / six
            this%m_a3(i) = (this%m_x(i+1) - this%m_x(i-1)) / three
            this%m_a4(i) = (this%m_x(i) - this%m_x(i-1)) / six
        end do

        ! Set the last equation
        select case (ibcend)
        case (SPLINE_QUADRATIC_OVER_INTERVAL)
            this%m_b(n) = zero
            this%m_a2(n) = -one
            this%m_a3(n) = one
        case (SPLINE_KNOWN_FIRST_DERIVATIVE)
            this%m_b(n) = ybcend - (this%m_y(n) - this%m_y(n-1)) / &
                (this%m_x(n) - this%m_x(n-1))
            this%m_a2(n) = (this%m_x(n) - this%m_x(n-1)) / six
            this%m_a3(n) = (this%m_x(n) - this%m_x(n-1)) / three
        case (SPLINE_KNOWN_SECOND_DERIVATIVE)
            this%m_b(n) = ybcend
            this%m_a2(n) = zero
            this%m_a3(n) = one
        case (SPLINE_CONTINUOUS_THIRD_DERIVATIVE)
            this%m_b(n) = zero
            this%m_a1(n) = this%m_x(n-1) - this%m_x(n)
            this%m_a2(n) = this%m_x(n) - this%m_x(n-2)
            this%m_a3(n) = this%m_x(n-2) - this%m_x(n-1)
        case default
            this%m_b(n) = zero
            this%m_a2(n) = -one
            this%m_a3(n) = one
        end select

        ! Define the second derivative
        if (n == 2 .and. ibcbeg == 0 .and. ibcend == 0) then
            ! Deal with the special case of N = 2, and IBCBEG = IBCEND = 0
            this%m_ypp(1) = zero
            this%m_ypp(2) = zero
        else
            ! Solve the linear system
            call penta_solve(this%m_a1, this%m_a2, this%m_a3, this%m_a4, &
                this%m_a5, this%m_b, this%m_ypp)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Initializes the specified spline_interp instance.  The end points
    !! are considered free such that the interpolant is quadratic over both the
    !! initial and final intervals.
    !!
    !! @param[in,out] this The spline_interp instance.
    !! @param[in] x An N-element array containing the independent variable data.
    !! @param[in] y An N-element array containing the dependent variable data.
    !! @param[in] order The order of the interpolating polynomial.  This 
    !!  parameter is ignored as the spline is a cubic approximation.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
    subroutine si_init_1(this, x, y, order, err)
        ! Arguments
        class(spline_interp), intent(inout) :: this
        real(dp), intent(in), dimension(:) :: x, y
        integer(i32), intent(in) :: order
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(i32) :: dummy

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        dummy = order ! Avoids complaining by the compiler

        ! Initialize the base object
        call im_init(this, x, y, 3, err)

        ! Evaluate the second derivatives
        call this%compute_diff2(SPLINE_QUADRATIC_OVER_INTERVAL, zero, &
            SPLINE_QUADRATIC_OVER_INTERVAL, zero, errmgr)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Initializes the specified spline_interp instance.
    !!
    !! @param[in,out] this The spline_interp instance.
    !! @param[in] x An N-element array containing the independent variable data.
    !! @param[in] y An N-element array containing the dependent variable data.
    !! @param[in] ibcbeg An optional input that defines the nature of the 
    !!  boundary condition at the beginning of the spline.  If no parameter, or 
    !!  an invalid parameter, is specified, the default free condition 
    !!  (SPLINE_QUADRATIC_OVER_INTERVAL) is used.
    !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
    !!      initial interval.  No value is required for @p ybcbeg.
    !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
    !!      initial point is provided in @p ybcbeg.
    !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
    !!      initial point is provided in @p ybcbeg.
    !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
    !!      at x(2).  No value is required for @p ybcbeg.
    !! @param[in] ybcbeg If needed, the value of the initial point boundary
    !!  condition.  If needed, but not supplied, a default value of zero will
    !!  be used.
    !! @param[in] ibcend An optional input that defines the nature of the 
    !!  boundary condition at the end of the spline.  If no parameter, or an 
    !!  invalid parameter, is specified, the default free condition 
    !!  (SPLINE_QUADRATIC_OVER_INTERVAL) is used.
    !!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
    !!      final interval.  No value is required for @p ybcend.
    !!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
    !!      initial point is provided in @p ybcend.
    !!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
    !!      initial point is provided in @p ybcend.
    !!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
    !!      at x(n-1).  No value is required for @p ybcend.
    !! @param[in] ybcend If needed, the value of the final point boundary
    !!  condition.  If needed, but not supplied, a default value of zero will
    !!  be used.
    !!
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
    !!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - CF_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
    subroutine si_init_2(this, x, y, ibcbeg, ybcbeg, ibcend, ybcend, err)
        ! Arguments
        class(spline_interp), intent(inout) :: this
        real(dp), intent(in), dimension(:) :: x, y
        integer(i32), intent(in), optional :: ibcbeg, ibcend
        real(dp), intent(in), optional :: ybcbeg, ybcend
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: ibeg, iend
        real(dp) :: ybeg, yend
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        ibeg = SPLINE_QUADRATIC_OVER_INTERVAL
        iend = SPLINE_QUADRATIC_OVER_INTERVAL
        ybeg = zero
        yend = zero
        if (present(ibcbeg)) ibeg = ibcbeg
        if (present(ybcbeg)) ybeg = ybcbeg
        if (present(ibcend)) iend = ibcend
        if (present(ybcend)) yend = ybcend

        ! Input Check
        if (ibeg /= SPLINE_CONTINUOUS_THIRD_DERIVATIVE .or. &
            ibeg /= SPLINE_KNOWN_SECOND_DERIVATIVE .or. &
            ibeg /= SPLINE_KNOWN_FIRST_DERIVATIVE .or. &
            ibeg /= SPLINE_QUADRATIC_OVER_INTERVAL) &
                ibeg = SPLINE_QUADRATIC_OVER_INTERVAL
        if (iend /= SPLINE_CONTINUOUS_THIRD_DERIVATIVE .or. &
            iend /= SPLINE_KNOWN_SECOND_DERIVATIVE .or. &
            iend /= SPLINE_KNOWN_FIRST_DERIVATIVE .or. &
            iend /= SPLINE_QUADRATIC_OVER_INTERVAL) &
                iend = SPLINE_QUADRATIC_OVER_INTERVAL

        ! Initialize the base object
        call im_init(this, x, y, 3, err)

        ! Evaluate the second derivatives
        call this%compute_diff2(ibeg, ybeg, iend, yend, errmgr)
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
