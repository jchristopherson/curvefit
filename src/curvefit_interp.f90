! curvefit_interp.f90

module curvefit_interp
    use curvefit_core
    implicit none

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
    integer(i32) :: m_npts
    integer(i32) :: m_order
    integer(i32) :: m_savedIndex
    integer(i32) :: m_indexCheck
    logical :: m_correlated
contains
    !> @brief Initializes the specified interp_manager instance.
    procedure, public :: initialize => im_init
    !> @brief Attempts to locate the index in the array providing a lower bounds
    !!  to the specified interpolation point.
    procedure, non_overridable, public :: locate => im_locate
    !> @brief Attempts to locate the index in the array providing a lower bounds
    !!  to the specified interpolation point.
    procedure, non_overridable, public :: hunt => interp_hunt
    !> @brief Interpolates to obtain the function value at the specified
    !!  independent variable.
    generic, public :: interpolate => im_perform, im_perform_array
    !> @brief Performs the actual interpolation.
    procedure(interp_1d), deferred :: raw_interp

    procedure, non_overridable :: im_perform
    procedure, non_overridable :: im_perform_array
end type

! ******************************************************************************
! ABSTRACT INTERFACES
! ------------------------------------------------------------------------------
interface
    !> @brief Defines the signature of a method used to interpolate a single
    !!  value in a 1D data set.
    !!
    !! @param[in,out] this The InterpolateManager instance.
    !! @param[in] jlo The array index below which @p pt is found in @p x.
    !! @param[in] x An N-element array of the independent data points.
    !! @param[in] y An N-element array of the corresponding dependent data
    !!  points.
    !! @param[in] pt The independent variable value to interpolate.
    !!
    !! @return The interpolated value.
    function interp_1d(this, jlo, x, y, pt) result(yy)
        use curvefit_core, only : dp, i32
        import interp_manager
        class(interp_manager), intent(inout) :: this
        integer(i32), intent(in) :: jlo
        real(dp), intent(in), dimension(:) :: x, y
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
    !! @param[in] npts The number of data points.
    !! @param[in] order The order of the interpolating polynomial.
    subroutine im_init(this, npts, order)
        class(interp_manager), intent(inout) :: this
        integer(i32), intent(in) :: npts, order
        this%m_npts = npts
        this%m_order = order
        this%m_savedIndex = 1
        this%m_indexCheck = 1
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Attempts to locate the index in the array providing a lower bounds
    !!  to the specified interpolation point.
    !!
    !! @param[in,out] this The interp_manager instance.
    !! @param[in] x The array of independent variable data.
    !! @param[in] pt The interpolation point.
    !!
    !! @return The array index in @p x below @p pt.
    function im_locate(this, x, pt) result(j)
        ! Arguments
        class(interp_manager), intent(inout) :: this
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(in) :: pt
        integer :: j

        ! Local Variables
        integer(i32) :: n, m, jhi, jmid, jlo
        logical :: ascnd

        ! Initialization
        n = this%m_npts
        m = this%m_order + 1
        ascnd = x(n) >= x(1)
        jlo = 1
        jhi = n

        ! Process
        do while (jhi - jlo > 1)
            jmid = (jhi + jlo) / 2
            if (pt >= x(jmid) .eqv. ascnd) then
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
    !! @param[in] x The array of independent variable data.
    !! @param[in] pt The interpolation point.
    !!
    !! @return The array index in @p x below @p pt.
    function im_hunt(this, x, pt) result(j)
        ! Arguments
        class(interp_manager), intent(inout) :: this
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(in) :: pt
        integer(i32) :: j

        ! Local Variables
        integer(i32) :: jlo, jmid, jhi, inc, n, m
        logical :: ascnd

        ! Initialization
        n = this%m_npts
        m = this%m_order + 1
        jlo = this%m_savedIndex
        inc = 1
        ascnd = x(n) > x(1)

        ! Process
        if (jlo < 1 .or. jlo > n) then
            jlo = 1
            jhi = n
        else
            if (pt >= x(jlo) .eqv. ascnd) then
                do
                    jhi = jlo + inc
                    if (jhi >= n) then
                        jhi = n
                        exit
                    else if (pt < x(jhi) .eqv. ascnd) then
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
                    else if (pt >= x(jlo) .eqv. ascnd) then
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
            if (pt >= x(jmid) .eqv. ascnd) then
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
    !! @param[in] x An N-element array of the independent data points.
    !! @param[in] y An N-element array of the corresponding dependent data
    !!  points.
    !! @param[in] pt The independent variable value to interpolate.
    !!
    !! @return The interpolated value.
    function im_perform(this, x, y, pt) result(yy)
        ! Arguments
        class(interp_manager), intent(inout) :: this
        real(dp), intent(in), dimension(:) :: x, y
        real(dp), intent(in) :: pt
        real(dp) :: yy

        ! Local Variables
        integer(i32) :: jlo

        ! Process
        if (this%m_correlated) then
            jlo = this%hunt(x, pt)
        else
            jlo = this%locate(x, pt)
        end if
        yy = this%raw_interp(jlo, x, y, pt)
    end function

! ------------------------------------------------------------------------------
    !> @brief Interpolates to obtain the function value at the specified
    !!  independent variables.
    !!
    !! @param[in,out] this The interp_manager instance.
    !! @param[in] x An N-element array of the independent data points.
    !! @param[in] y An N-element array of the corresponding dependent data
    !!  points.
    !! @param[in] pt An M-element array containing the independent variable 
    !!  values to interpolate.
    !!
    !! @return An M-element array containing the interpolated values.
    function im_perform_array(this, x, y, pts) result(yy)
        ! Arguments
        class(interp_manager), intent(inout) :: this
        real(dp), intent(in), dimension(:) :: x, y
        real(dp), intent(in), dimension(:) :: pts
        real(dp), dimension(size(pts)) :: yy

        ! Local Variables
        integer(i32) :: i, jlo

        ! Process
        do i = 1, size(pts)
            if (this%m_correlated) then
                jlo = this%hunt(x, pts(i))
            else
                jlo = this%locate(x, pts(i))
            end if
            yy(i) = this%raw_interp(jlo, x, y, pts(i))
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
