module ufield_smooth_class

use ufield_class
use param
use system
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
  real, dimension(:), pointer :: norm_bnd => null()

  contains

  generic :: new => init_ufield_smooth
  generic :: del => end_ufield_smooth
  procedure, private :: init_ufield_smooth
  procedure, private :: end_ufield_smooth
  procedure :: smooth_f1
  procedure :: if_smooth

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

    ! normalization factor at boundary
    total = 0.0
    allocate( this%norm_bnd( 1:this%order ) )
    do j = -this%order, 0
      total = total + this%scoef(j)
    enddo
    this%norm_bnd(1) = total
    do j = 1, this%order-1
      total = total + this%scoef(j)
      this%norm_bnd(j+1) = total
    enddo
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
  if ( associated(this%norm_bnd) ) deallocate( this%norm_bnd )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_ufield_smooth

!-----------------------------------------------------------------------------------------
subroutine smooth_f1( this, uf )
!-----------------------------------------------------------------------------------------

  implicit none

  class( ufield_smooth ), intent(inout) :: this
  class( ufield ), intent(inout) :: uf

  real, dimension(:,:), pointer :: pf => null()

  integer :: ext_lb, ext_ub, int_lb, int_ub
  integer :: nsm, dim, j, k, idproc, nvp, nrp

  real, dimension(:,:), allocatable, save :: o
  real, dimension(:), allocatable, save :: f

  character(len=20), save :: sname = "smooth_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  pf => uf%get_f1()

  ! exterior boundaries of grid
  ext_lb = lbound( pf, 2 )
  ext_ub = ubound( pf, 2 )

  nsm    = this%order
  dim    = uf%get_dim()
  idproc = uf%pp%getlidproc()
  nvp    = uf%get_nvp(1)
  nrp    = uf%get_ndp(1)

  ! interior boundaries of grid
  int_lb = 1
  int_ub = nrp
  if ( idproc == 0 )   int_lb = 1 + nsm
  if ( idproc == nvp-1 ) int_ub = nrp - nsm

  if ( .not. allocated(o) ) then
    allocate( o(dim,-nsm:nsm), f(dim) )
  endif
  o = 0.0
  f = 0.0

  if ( nsm > 0 ) then

    do k = -nsm, nsm-1
      o(:,k) = pf(:, 1+k)
    enddo

    ! smooth the axis
    if ( idproc == 0 ) then

      do j = 1, int_lb-1

        o(:, nsm) = pf(:, j+nsm)
        f(:) = o(:,1-j) * this%scoef(1-j)

        do k = 2-j, nsm
          f(:) = f(:) + o(:,k) * this%scoef(k)
        enddo

        do k = -j, nsm-1
          o(:,k) = o(:,k+1)
        enddo

        pf(:,j) = f(:) / this%norm_bnd(j)

      enddo
    endif

    ! smooth interior region
    do j = int_lb, int_ub

      o(:, nsm) = pf(:, j+nsm)

      f(:) = o(:, -nsm) * this%scoef(-nsm)

      do k = -nsm + 1, nsm
        f(:) = f(:) + o(:,k) * this%scoef(k)
      enddo

      do k = -nsm, nsm-1
        o(:,k) = o(:,k+1)
      enddo

      pf(:,j) = f(:)

    enddo

    ! smooth outer boundary
    if ( idproc == nvp-1 ) then
      do j = int_ub+1, nrp

        ! o(:, nsm) = pf(:, j+nsm) ! no need for this
        f(:) = o(:,-nsm) * this%scoef(-nsm)

        do k = -nsm+1, nsm-j+int_ub
          f(:) = f(:) + o(:,k) * this%scoef(k)
        enddo

        do k = -nsm, nsm-j+int_ub-1
          o(:,k) = o(:,k+1)
        enddo

        pf(:,j) = f(:) / this%norm_bnd(nrp-j+1)

      enddo
    endif

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine smooth_f1
!-----------------------------------------------------------------------------------------

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

end module ufield_smooth_class