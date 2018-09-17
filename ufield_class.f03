module ufield_class

use param

implicit none

private

public :: ufield

type :: ufield

  private

  real, dimension(:,:), pointer, public :: f1 => null()
  real, dimension(:,:,:), pointer :: f2 => null()
  integer, dimension(2) :: nd ! number of global grid points
  integer, dimension(2) :: ndp ! number of local grid points
  integer, dimension(2) :: nvp
  integer :: dim ! dimension of array (guard cell excluded)
  integer, dimension(2) :: noff
  integer, dimension(2,2) :: gc_num ! number of guard cells
  logical :: has_2d

  real, dimension(:), pointer :: buf => null() ! data buffer used for MPI

  contains

  generic :: new => init_ufield
  generic :: del => end_ufield
  ! generic :: write_hdf5
  generic :: get_nd => get_nd_all, get_nd_dim
  generic :: get_ndp => get_ndp_all, get_ndp_dim
  generic :: get_gc_num => get_gc_num_all, get_gc_num_dim
  generic :: get_nvp => get_nvp_all, get_nvp_dim
  procedure :: get_dim
  procedure :: copy_slice
  procedure :: get_f1
  procedure :: get_f2

  generic :: assignment(=) => assign_array
  generic :: operator(+) => add_array, add_scalar1, add_scalar2
  generic :: operator(-) => sub_array, sub_scalar1, sub_scalar2
  generic :: operator(*) => dot_array, dot_scalar1, dot_scalar2

  procedure, private :: init_ufield, end_ufield
  procedure, private :: get_nd_all, get_nd_dim
  procedure, private :: get_ndp_all, get_ndp_dim
  procedure, private :: get_gc_num_all, get_gc_num_dim
  procedure, private :: get_nvp_all, get_nvp_dim

  procedure, private, pass(a1) :: add_array, add_scalar1, add_scalar2
  procedure, private, pass(a1) :: dot_array, dot_scalar1, dot_scalar2
  procedure, private, pass(a1) :: sub_array, sub_scalar1, sub_scalar2
  procedure, private :: assign_array

end type ufield

contains

subroutine init_ufield( this, dim, nd, nvp, gc_num, has_2d )

  implicit none

  class( ufield ), intent(inout) :: this
  integer, intent(in) :: dim
  integer, intent(in), dimension(2) :: nd, nvp
  integer, intent(in), dimension(2,2) :: gc_num
  logical, intent(in), optional :: has_2d

  character(len=32), save :: sname = "init_ufield:"

  this%nd = nd
  this%nvp = nvp
  this%dim = dim
  this%ndp = nd / this%nvp
  this%noff = (/0,0/) ! single core
  this%gc_num = gc_num

  allocate( this%f1( dim, 1-this%gc_num(p_lower,1):this%ndp(1)+this%gc_num(p_upper,1) ) )
  this%f1 = 0.0

  if ( present(has_2d) ) then
    this%has_2d = has_2d
    allocate( this%f2( dim, &
              1-this%gc_num(p_lower,1):this%ndp(1)+this%gc_num(p_upper,1), &
              1-this%gc_num(p_lower,2):this%ndp(2)+this%gc_num(p_upper,2) ) )
    this%f2 = 0.0
  else
    this%has_2d = .false.
  endif

end subroutine

subroutine end_ufield( this )

  implicit none

  class( ufield ), intent(inout) :: this

  character(len=32), save :: sname = "end_ufield:"

  if ( associated( this%f1 ) ) deallocate( this%f1 )
  if ( associated( this%f2 ) ) deallocate( this%f2 )
  
end subroutine end_ufield

subroutine copy_slice( this, idx, dir )

  implicit none

  class( ufield ), intent(inout) :: this
  integer, intent(in) :: idx, dir

  integer :: i, j, lb, ub

  lb = 1-this%gc_num(p_lower,1)
  ub = this%ndp(1)+this%gc_num(p_upper,1)

  if ( .not. this%has_2d ) then
    print *, "The field has no 2D layout."
    stop
  endif

  select case ( dir )

  case ( p_copy_1to2 )

    do j = lb, ub
      do i = 1, this%dim
        this%f2(i,j,idx) = this%f1(i,j)
      enddo
    enddo

  case ( p_copy_2to1 )

    do j = lb, ub
      do i = 1, this%dim
        this%f1(i,j) = this%f2(i,j,idx)
      enddo
    enddo

  end select
  
end subroutine copy_slice

subroutine assign_array( this, that )

  implicit none

  class( ufield ), intent(inout) :: this
  class(*), intent(in) :: that

  integer :: i, j

  select type (that)
    
    type is (real)

      do j = 1, this%ndp(1)
        do i = 1, this%dim
          this%f1(i,j) = that
        enddo
      enddo

    class is (ufield)

      do j = 1, this%ndp(1)
        do i = 1, this%dim
          this%f1(i,j) = that%f1(i,j)
        enddo
      enddo

    class default

      print *, "invalid assignment type!"
      stop

  end select

end subroutine assign_array

function add_array( a1, a2 ) result( a3 )

  implicit none

  class( ufield ), intent(in) :: a1
  class( ufield ), intent(in) :: a2
  class( ufield ), allocatable :: a3

  integer :: i, j

  allocate(a3)
  call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )

  do j = 1, a1%ndp(1)
    do i = 1, a1%dim
      a3%f1(i,j) = a1%f1(i,j) + a2%f1(i,j)
    enddo
  enddo

end function add_array

function add_scalar1( a1, a2 ) result( a3 )

  implicit none

  class( ufield ), intent(in) :: a1
  real, intent(in) :: a2
  class( ufield ), allocatable :: a3

  integer :: i, j

  allocate( a3 )
  ! note that buffer has no 2D array
  call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )

  do j = 1, a1%ndp(1)
    do i = 1, a1%dim
      a3%f1(i,j) = a1%f1(i,j) + a2
    enddo
  enddo

end function add_scalar1

function add_scalar2( a2, a1 ) result( a3 )

  implicit none

  real, intent(in) :: a2
  class( ufield ), intent(in) :: a1
  class( ufield ), allocatable :: a3

  integer :: i, j

  allocate( a3 )
    ! note that buffer has no 2D array
  call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )

  do j = 1, a1%ndp(1)
    do i = 1, a1%dim
      a3%f1(i,j) = a1%f1(i,j) + a2
    enddo
  enddo
    
end function add_scalar2

function sub_array( a1, a2 ) result( a3 )

  implicit none

  class( ufield ), intent(in) :: a1
  class( ufield ), intent(in) :: a2
  class( ufield ), allocatable :: a3

  integer :: i, j

  allocate(a3)
  call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )

  do j = 1, a1%ndp(1)
    do i = 1, a1%dim
      a3%f1(i,j) = a1%f1(i,j) - a2%f1(i,j)
    enddo
  enddo

end function sub_array

function sub_scalar1( a1, a2 ) result( a3 )

  implicit none

  class( ufield ), intent(in) :: a1
  real, intent(in) :: a2
  class( ufield ), allocatable :: a3

  integer :: i, j

  allocate( a3 )
  ! note that buffer has no 2D array
  call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )

  do j = 1, a1%ndp(1)
    do i = 1, a1%dim
      a3%f1(i,j) = a1%f1(i,j) - a2
    enddo
  enddo

end function sub_scalar1

function sub_scalar2( a2, a1 ) result( a3 )

  implicit none

  real, intent(in) :: a2
  class( ufield ), intent(in) :: a1
  class( ufield ), allocatable :: a3

  integer :: i, j

  allocate( a3 )
    ! note that buffer has no 2D array
  call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )

  do j = 1, a1%ndp(1)
    do i = 1, a1%dim
      a3%f1(i,j) = a2 - a1%f1(i,j)
    enddo
  enddo
    
end function sub_scalar2

function dot_array( a1, a2 ) result( a3 )

  implicit none

  class( ufield ), intent(in) :: a1
  class( ufield ), intent(in) :: a2
  class( ufield ), allocatable :: a3

  integer :: i, j

  allocate(a3)
  call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )

  do j = 1, a1%ndp(1)
    do i = 1, a1%dim
      a3%f1(i,j) = a1%f1(i,j) * a2%f1(i,j)
    enddo
  enddo

end function dot_array

function dot_scalar1( a1, a2 ) result( a3 )

  implicit none

  class( ufield ), intent(in) :: a1
  real, intent(in) :: a2
  class( ufield ), allocatable :: a3

  integer :: i, j

  allocate( a3 )
  ! note that buffer has no 2D array
  call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )

  do j = 1, a1%ndp(1)
    do i = 1, a1%dim
      a3%f1(i,j) = a1%f1(i,j) * a2
    enddo
  enddo

end function dot_scalar1

function dot_scalar2( a2, a1 ) result( a3 )

  implicit none

  real, intent(in) :: a2
  class( ufield ), intent(in) :: a1
  class( ufield ), allocatable :: a3

  integer :: i, j

  allocate( a3 )
    ! note that buffer has no 2D array
  call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )

  do j = 1, a1%ndp(1)
    do i = 1, a1%dim
      a3%f1(i,j) = a1%f1(i,j) * a2
    enddo
  enddo
    
end function dot_scalar2

function get_dim( this )

  implicit none

  class( ufield ), intent(in) :: this
  integer :: get_dim

  get_dim = this%dim

end function get_dim

function get_nd_all( this )

  implicit none

  class( ufield ), intent(in) :: this
  integer, dimension(2) :: get_nd_all

  get_nd_all = this%nd

end function get_nd_all

function get_nd_dim( this, dim )

  implicit none

  class( ufield ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_nd_dim

  get_nd_dim = this%nd(dim)
  
end function get_nd_dim

function get_ndp_all( this )

  implicit none

  class( ufield ), intent(in) :: this
  integer, dimension(2) :: get_ndp_all

  get_ndp_all = this%ndp

end function get_ndp_all

function get_ndp_dim( this, dim )

  implicit none

  class( ufield ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_ndp_dim

  get_ndp_dim = this%ndp(dim)
  
end function get_ndp_dim

function get_gc_num_all( this )

  implicit none

  class( ufield ), intent(in) :: this
  integer, dimension(2,2) :: get_gc_num_all

  get_gc_num_all = this%gc_num

end function get_gc_num_all

function get_gc_num_dim( this, dim )

  implicit none

  class( ufield ), intent(in) :: this
  integer, intent(in) :: dim
  integer, dimension(2) :: get_gc_num_dim

  get_gc_num_dim = this%gc_num(:,dim)
  
end function get_gc_num_dim

function get_nvp_all( this )

  implicit none

  class( ufield ), intent(in) :: this
  integer, dimension(2) :: get_nvp_all

  get_nvp_all = this%nvp

end function get_nvp_all

function get_nvp_dim( this, dim )

  implicit none

  class( ufield ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_nvp_dim

  get_nvp_dim = this%nvp(dim)
  
end function get_nvp_dim

function get_f1( this )

  implicit none

  class( ufield ), intent(in) :: this
  real, dimension(:,:), pointer :: get_f1

  get_f1 => this%f1

end function get_f1

function get_f2( this )

  implicit none

  class( ufield ), intent(in) :: this
  real, dimension(:,:,:), pointer :: get_f2

  get_f2 => this%f2
  
end function get_f2

end module ufield_class