! curvefit_interp.f90

module curvefit_interp
    use curvefit_core
    use ferror, only : errors
    implicit none
    private
    public :: interp_manager
    public :: linear_interp
    public :: polynomial_interp

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
        !> @brief Initializes the specified interp_manager instance.
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
        !> @brief Initializes the specified polynomial_interp instance.
        procedure, public :: initialize => pi_init
        !> @brief Performs the actual interpolation.
        procedure :: raw_interp => pi_raw_interp
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
    !> @brief Performs the actual linear interpolation.
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

        ! Process
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

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
