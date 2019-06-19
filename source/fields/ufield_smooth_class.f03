module ufield_smooth_class

use ufield_class
use param
use sys
use mpi

implicit none

private

character(len=20), parameter :: cls_name = "ufield_smooth"
integer, parameter :: cls_level = 4

public :: ufield_smooth

type :: ufield_smooth

  private

  integer :: type
  integer :: order
  real, dimension(:), pointer :: scoef => null()

  contains

  generic :: new => init_ufield_smooth
  generic :: del => end_ufield_smooth
  procedure, private :: init_ufield_smooth
  procedure, private :: end_ufield_smooth
  procedure :: smooth_f1
  procedure :: if_smooth
  procedure :: get_scoef, get_order

end type ufield_smooth

contains

subroutine init_ufield_smooth( this, type, order )
!-----------------------------------------------------------------------------------------

  implicit none

  class( ufield_smooth ), intent(inout) :: this
  integer, intent(in) :: type
  integer, intent(in) :: order

  real, dimension(:,:), allocatable :: kern3
  real :: comp, total

  integer :: i, j

  character(len=20), save :: sname = "init_ufield_smooth"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%type = type
  this%order = order

  select case ( this%type )

  case ( p_smooth_binomial )

    allocate( kern3( 3, this%order ) )
    allocate( this%scoef( -this%order:this%order ) )
    do j = 1, this%order
      kern3(:,j) = (/ 1.0d0, 2.0d0, 1.0d0 /)
    enddo
    this%scoef = combine_kernel3( kern3, this%order )

  case ( p_smooth_compensated )

    allocate( kern3( 3, this%order ) )
    allocate( this%scoef( -this%order:this%order ) )
    do j = 1, this%order - 1
      kern3(:,j) = (/  1.0d0, 2.0d0, 1.0d0 /)
    enddo

    comp = ( 4.0 + 2.0 * ( this%order - 1 ) ) / ( this%order - 1 )

    kern3(:,this%order) = (/ -1.0d0, comp, -1.0d0 /)

    this%scoef = combine_kernel3( kern3, this%order )

  case ( p_smooth_none )

    this%order = 0

  case default

    call write_err( 'Invalid smooth type!' )

  end select

  ! normalize coefficients
  if ( this%order > 0 ) then
    total = 0.0
    do j = -this%order, this%order
     total = total + this%scoef(j)
    enddo
    this%scoef = this%scoef / total
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_ufield_smooth

!-----------------------------------------------------------------------------------------
subroutine end_ufield_smooth( this )

  implicit none

  class( ufield_smooth ), intent(inout) :: this

  character(len=20), save :: sname = "end_ufield_smooth"
  
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( associated(this%scoef) ) deallocate( this%scoef )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_ufield_smooth

!-----------------------------------------------------------------------------------------
subroutine smooth_f1( this, uf, mode )
!-----------------------------------------------------------------------------------------

  implicit none

  class( ufield_smooth ), intent(inout) :: this
  class( ufield ), intent(inout) :: uf
  integer, intent(in) :: mode
  ! logical, intent(in) :: q_cons
  ! real, intent(inout), optional :: q_ax

  real, dimension(:,:), pointer :: pf => null()

  integer :: ext_lb, ext_ub, int_lb, int_ub
  integer :: nsm, dim, i, j, k, idproc, nvp, nrp, noff
  real :: pha

  real, dimension(:,:), allocatable, save :: o
  real, dimension(:), allocatable, save :: f

  character(len=20), save :: sname = "smooth_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call start_tprof( 'smooth' )

  pf => uf%get_f1()

  ! exterior boundaries of grid
  ext_lb = lbound( pf, 2 )
  ext_ub = ubound( pf, 2 )

  nsm    = this%order
  dim    = uf%get_dim()
  idproc = uf%pp%getlidproc()
  nvp    = uf%get_nvp(1)
  nrp    = uf%get_ndp(1)
  noff   = uf%get_noff(1)

  ! interior boundaries of grid
  int_lb = 1
  int_ub = nrp
  if ( idproc == 0 )   int_lb = 1 + nsm
  if ( idproc == nvp-1 ) int_ub = nrp - nsm

  if ( .not. allocated(o) ) then
    allocate( o(p_max_xdim,-nsm:nsm), f(p_max_xdim) )
  endif
  o = 0.0
  f = 0.0

  if ( nsm > 0 ) then

    do k = -nsm, nsm-1
    ! do k = 0, nsm-1
      do i = 1, dim
        o(i,k) = pf(i, 1+k)
      enddo
    enddo

    ! smooth the axis
    if ( idproc == 0 ) then

      if ( mod(mode,2) == 0 ) then
        pha = 1.0
      else
        pha = -1.0
      endif

      do k = -nsm, -1
        do i = 1, dim
          o(i,k) = pha * pf(i,-k)
        enddo
      enddo

      do j = 1, int_lb-1
        do i = 1, dim

          o(i, nsm) = pf(i, j+nsm)
          
          f(i) = 0.0
          ! do k = 1-j, nsm
          do k = -nsm, nsm
            f(i) = f(i) + o(i,k) * this%scoef(k)
          enddo

          ! do k = -j, nsm-1
          do k = -nsm, nsm-1
            o(i,k) = o(i,k+1)
          enddo

          pf(i,j) = f(i)

        enddo
      enddo

      ! if ( q_cons ) then

      !   if ( .not. present(q_ax) ) call write_err( 'Argument q_ax must be provided.' )
      !   if ( dim /= 1 ) call write_err( 'Charge conserved smooth only applies for scalar field.' )

      !   q_ax = 0.0

      !   ! add the charge smoothed out
      !   do j = 1, nsm
      !     do i = -nsm, -j
      !       q_ax = q_ax + pf(1,j) * this%scoef(i) * real(i+j-0.5)
      !     enddo
      !   enddo

      !   ! deduct the charge smoothed in
      !   do j = 1, nsm
      !     do i = 1, nsm-j+1
      !       q_ax = q_ax - real(j-0.5) * pf(1,i) * this%scoef(i)
      !     enddo
      !   enddo

      ! endif

    endif

    ! smooth interior region
    do j = int_lb, int_ub
      do i = 1, dim

        o(i, nsm) = pf(i, j+nsm)

        f(i) = 0.0
        do k = -nsm, nsm
          f(i) = f(i) + o(i,k) * this%scoef(k)
        enddo

        do k = -nsm, nsm-1
          o(i,k) = o(i,k+1)
        enddo

        pf(i,j) = f(i)

      enddo
    enddo

    ! smooth outer boundary
    if ( idproc == nvp-1 ) then
      do j = int_ub+1, nrp
        do i = 1, dim

          k = j - int_ub
          o(i, nsm) = pf(i, nrp-k+1)

          f(i) = 0.0
          ! do k = -nsm, nsm-j+int_ub
          do k = -nsm, nsm
            f(i) = f(i) + o(i,k) * this%scoef(k)
          enddo

          do k = -nsm, nsm-1
            o(i,k) = o(i,k+1)
          enddo

          pf(i,j) = f(i)

        enddo
      enddo
    endif

  endif

  call stop_tprof( 'smooth' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine smooth_f1
!-----------------------------------------------------------------------------------------

! !-----------------------------------------------------------------------------------------
! subroutine smooth_f1( this, uf, dr )
! !-----------------------------------------------------------------------------------------

!   implicit none

!   class( ufield_smooth ), intent(inout) :: this
!   class( ufield ), intent(inout) :: uf
!   real, intent(in) :: dr

!   real, dimension(:,:), pointer :: pf => null()

!   integer :: ext_lb, ext_ub, int_lb, int_ub
!   integer :: nsm, dim, i, j, k, idproc, nvp, nrp, noff
!   real :: r

!   real, dimension(:,:), allocatable, save :: o
!   real, dimension(:), allocatable, save :: f

!   character(len=20), save :: sname = "smooth_f1"

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   call start_tprof( 'smooth' )

!   pf => uf%get_f1()

!   ! exterior boundaries of grid
!   ext_lb = lbound( pf, 2 )
!   ext_ub = ubound( pf, 2 )

!   nsm    = this%order
!   dim    = uf%get_dim()
!   idproc = uf%pp%getlidproc()
!   nvp    = uf%get_nvp(1)
!   nrp    = uf%get_ndp(1)
!   noff   = uf%get_noff(1)

!   ! interior boundaries of grid
!   int_lb = 1
!   int_ub = nrp
!   if ( idproc == 0 )   int_lb = 1 + nsm
!   if ( idproc == nvp-1 ) int_ub = nrp - nsm

!   if ( .not. allocated(o) ) then
!     allocate( o(p_max_xdim,-nsm:nsm), f(p_max_xdim) )
!   endif
!   o = 0.0
!   f = 0.0

!   if ( nsm > 0 ) then

!     do k = -nsm, nsm-1
!     ! do k = 0, nsm-1
!       do i = 1, dim
!         r = ( noff + 0.5 + k ) * dr
!         o(i,k) = pf(i, 1+k) * r
!       enddo
!     enddo

!     ! smooth the axis
!     if ( idproc == 0 ) then

!       do k = -nsm, -1
!         do i = 1, dim
!           o(i,k) = o(i,-k-1)
!         enddo
!       enddo

!       do j = 1, int_lb-1
!         do i = 1, dim

!           r = (noff+j+nsm-0.5) * dr
!           o(i, nsm) = pf(i, j+nsm) * r
          
!           f(i) = 0.0
!           ! do k = 1-j, nsm
!           do k = -nsm, nsm
!             f(i) = f(i) + o(i,k) * this%scoef(k)
!           enddo

!           ! do k = -j, nsm-1
!           do k = -nsm, nsm-1
!             o(i,k) = o(i,k+1)
!           enddo

!           r = (noff+j-0.5) * dr
!           pf(i,j) = f(i) / r

!         enddo
!       enddo
!     endif

!     ! smooth interior region
!     do j = int_lb, int_ub
!       do i = 1, dim

!         r = (noff+j+nsm-0.5) * dr
!         o(i, nsm) = pf(i, j+nsm) * r

!         f(i) = 0.0
!         do k = -nsm, nsm
!           f(i) = f(i) + o(i,k) * this%scoef(k)
!         enddo

!         do k = -nsm, nsm-1
!           o(i,k) = o(i,k+1)
!         enddo

!         r = (noff+j-0.5) * dr
!         pf(i,j) = f(i) / r

!       enddo
!     enddo

!     ! smooth outer boundary
!     if ( idproc == nvp-1 ) then
!       do j = int_ub+1, nrp
!         do i = 1, dim

!           k = j - int_ub
!           r = (noff+nrp-k+0.5) * dr
!           o(i, nsm) = pf(i, nrp-k+1) * r

!           f(i) = 0.0
!           ! do k = -nsm, nsm-j+int_ub
!           do k = -nsm, nsm
!             f(i) = f(i) + o(i,k) * this%scoef(k)
!           enddo

!           do k = -nsm, nsm-1
!             o(i,k) = o(i,k+1)
!           enddo

!           r = (noff+j-0.5) * dr
!           pf(i,j) = f(i) / r

!         enddo
!       enddo
!     endif

!   endif

!   call stop_tprof( 'smooth' )
!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine smooth_f1
! !-----------------------------------------------------------------------------------------

! !-----------------------------------------------------------------------------------------
! subroutine smooth_f1( this, uf, dr )
! !-----------------------------------------------------------------------------------------

!   implicit none

!   class( ufield_smooth ), intent(inout) :: this
!   class( ufield ), intent(inout) :: uf
!   real, intent(in) :: dr

!   real, dimension(:,:), pointer :: pf => null()

!   integer :: ext_lb, ext_ub, int_lb, int_ub
!   integer :: nsm, dim, i, j, k, idproc, nvp, nrp, noff
!   real :: r

!   real, dimension(:,:), allocatable, save :: o
!   real, dimension(:), allocatable, save :: f, qloss, qcomp

!   character(len=20), save :: sname = "smooth_f1"

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   pf => uf%get_f1()

!   ! exterior boundaries of grid
!   ext_lb = lbound( pf, 2 )
!   ext_ub = ubound( pf, 2 )

!   nsm    = this%order
!   dim    = uf%get_dim()
!   idproc = uf%pp%getlidproc()
!   nvp    = uf%get_nvp(1)
!   nrp    = uf%get_ndp(1)
!   noff   = uf%get_noff(1)

!   ! interior boundaries of grid
!   int_lb = 1
!   int_ub = nrp
!   if ( idproc == 0 )   int_lb = 1 + nsm
!   if ( idproc == nvp-1 ) int_ub = nrp - nsm

!   if ( .not. allocated(o) ) then
!     allocate( o(p_max_xdim,-nsm:nsm), f(p_max_xdim), &
!       qloss(p_max_xdim), qcomp(p_max_xdim) )
!   endif
!   o = 0.0
!   f = 0.0

!   if ( nsm > 0 ) then

!     do k = -nsm, nsm-1
!     ! do k = 0, nsm-1
!       do i = 1, dim
!         o(i,k) = pf(i, 1+k)
!       enddo
!     enddo

!     ! smooth the axis
!     if ( idproc == 0 ) then

!       do k = -nsm, -1
!         do i = 1, dim
!           o(i,k) = o(i,-k-1)
!         enddo
!       enddo

!       qloss = 0.0
!       qcomp = 0.0

!       do j = 1, int_lb-1 ! 1, nsm
!         do i = 1, dim

!           o(i, nsm) = pf(i, j+nsm)
          
!           f(i) = 0.0
!           ! do k = 1-j, nsm
!           do k = -nsm, nsm
!             f(i) = f(i) + o(i,k) * this%scoef(k)
!           enddo

!           ! calculate charge loss and compensation
!           r = (j-0.5) * dr
!           do k = 1, nsm-j+1
!             qloss(i) = qloss(i) + o(i,0) * (r-dr) * (2-k-j) * this%scoef(k+j-1)
!             qcomp(i) = qcomp(i) + o(i,k-j) * this%scoef(k) * r
!           enddo

!           ! do k = -j, nsm-1
!           do k = -nsm, nsm-1
!             o(i,k) = o(i,k+1)
!           enddo

!           pf(i,j) = f(i)

!         enddo
!       enddo

!       pf(1:dim,0) = pf(1:dim,0) + ( qloss(1:dim) - qcomp(1:dim) ) / dr

!     endif

!     ! smooth interior region
!     do j = int_lb, int_ub
!       do i = 1, dim

!         o(i, nsm) = pf(i, j+nsm)

!         f(i) = 0.0
!         do k = -nsm, nsm
!           f(i) = f(i) + o(i,k) * this%scoef(k)
!         enddo

!         do k = -nsm, nsm-1
!           o(i,k) = o(i,k+1)
!         enddo

!         pf(i,j) = f(i)

!       enddo
!     enddo

!     ! smooth outer boundary
!     if ( idproc == nvp-1 ) then
!       do j = int_ub+1, nrp
!         do i = 1, dim

!           k = j - int_ub
!           o(i, nsm) = pf(i, nrp-k+1)

!           f(i) = 0.0
!           ! do k = -nsm, nsm-j+int_ub
!           do k = -nsm, nsm
!             f(i) = f(i) + o(i,k) * this%scoef(k)
!           enddo

!           do k = -nsm, nsm-1
!             o(i,k) = o(i,k+1)
!           enddo

!           pf(i,j) = f(i)

!         enddo
!       enddo
!     endif

!   endif

!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine smooth_f1
! !-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Combines n 3 point kernels into a single 2n+1 point kernel
!-----------------------------------------------------------------------------------------
function combine_kernel3( kern3, n )

  implicit none

  real, dimension(:,:), intent(in) :: kern3
  integer, intent(in) :: n

  real, dimension( 2*n+1 ) :: combine_kernel3
  real :: f1, f2
  integer :: k, i

  combine_kernel3 = 0.0
  combine_kernel3(n+1) = 1.0

  do k = 1, n
    f1 = 0.
    do i = (n+1) - (k-1), (n+1) + k
      combine_kernel3(i-1) = kern3(3,k) * combine_kernel3(i) + combine_kernel3(i-1)
      f2                   = kern3(1,k) * combine_kernel3(i)
      combine_kernel3(i)   = kern3(2,k) * combine_kernel3(i) + f1
      f1 = f2
    enddo
  enddo

end function combine_kernel3
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function if_smooth( this )
!-----------------------------------------------------------------------------------------

  implicit none

  logical :: if_smooth

  class( ufield_smooth ), intent(in) :: this

  if_smooth = ( this%order > 0 )

end function if_smooth
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_order( this )
!-----------------------------------------------------------------------------------------

  implicit none

  class( ufield_smooth ), intent(in) :: this
  integer :: get_order

  get_order = this%order

end function get_order
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_scoef( this )
!-----------------------------------------------------------------------------------------

  implicit none

  class( ufield_smooth ), intent(in) :: this
  real, dimension(:), pointer :: get_scoef

  get_scoef => this%scoef

end function get_scoef
!-----------------------------------------------------------------------------------------

end module ufield_smooth_class