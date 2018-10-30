module ufield_class

use parallel_pipe_class
use grid_class
use hdf5io_class
use param
use system
use mpi

implicit none

private

character(len=20), parameter :: cls_name = "ufield"
integer, parameter :: cls_level = 4

public :: ufield

type :: ufield

  private

  class( parallel_pipe ), pointer, public :: pp => null()
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

  generic :: new        => init_ufield, init_ufield_cp
  generic :: del        => end_ufield
  generic :: get_nd     => get_nd_all, get_nd_dim
  generic :: get_ndp    => get_ndp_all, get_ndp_dim
  generic :: get_gc_num => get_gc_num_all, get_gc_num_dim
  generic :: get_nvp    => get_nvp_all, get_nvp_dim
  generic :: get_noff   => get_noff_all, get_noff_dim
  generic :: copy_gc    => copy_gc_f1, copy_gc_f2
  generic :: write_hdf5 => write_hdf5_single, write_hdf5_pipe
  procedure :: get_dim
  procedure :: copy_slice
  procedure :: get_f1
  procedure :: get_f2

  generic :: assignment(=) => assign_array
  generic :: operator(+) => add_array, add_scalar1, add_scalar2
  generic :: operator(-) => sub_array, sub_scalar1, sub_scalar2
  generic :: operator(*) => dot_array, dot_scalar1, dot_scalar2

  procedure, private :: init_ufield, init_ufield_cp, end_ufield
  procedure, private :: get_nd_all, get_nd_dim
  procedure, private :: get_ndp_all, get_ndp_dim
  procedure, private :: get_gc_num_all, get_gc_num_dim
  procedure, private :: get_nvp_all, get_nvp_dim
  procedure, private :: get_noff_all, get_noff_dim
  procedure, private :: copy_gc_f1, copy_gc_f2
  procedure, private :: write_hdf5_single, write_hdf5_pipe

  procedure, private, pass(a1) :: add_array, add_scalar1, add_scalar2
  procedure, private, pass(a1) :: dot_array, dot_scalar1, dot_scalar2
  procedure, private, pass(a1) :: sub_array, sub_scalar1, sub_scalar2
  procedure, private :: assign_array

end type ufield

contains

subroutine init_ufield( this, pp, gp, dim, gc_num, has_2d )

  implicit none

  class( ufield ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: dim
  integer, intent(in), dimension(2,2) :: gc_num
  logical, intent(in), optional :: has_2d

  character(len=20), save :: sname = "init_ufield"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%pp => pp
  this%nd = gp%get_nd()
  this%nvp = gp%get_nvp()
  this%dim = dim
  this%ndp = gp%get_ndp()
  this%noff = gp%get_noff()
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

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_ufield

subroutine init_ufield_cp( this, that, has_2d )

  implicit none

  class( ufield ), intent(inout) :: this
  class( ufield ), intent(in) :: that
  logical, intent(in), optional :: has_2d

  character(len=20), save :: sname = "init_ufield_cp"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%pp => that%pp
  this%nd = that%get_nd()
  this%nvp = that%get_nvp()
  this%dim = that%dim
  this%ndp = that%get_ndp()
  this%noff = that%get_noff()
  this%gc_num = that%get_gc_num()

  allocate( this%f1( this%dim, 1-this%gc_num(p_lower,1):this%ndp(1)+this%gc_num(p_upper,1) ) )
  this%f1 = 0.0

  if ( present(has_2d) ) then
    this%has_2d = has_2d
    allocate( this%f2( this%dim, &
              1-this%gc_num(p_lower,1):this%ndp(1)+this%gc_num(p_upper,1), &
              1-this%gc_num(p_lower,2):this%ndp(2)+this%gc_num(p_upper,2) ) )
    this%f2 = 0.0
  else
    this%has_2d = .false.
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_ufield_cp

subroutine end_ufield( this )

  implicit none

  class( ufield ), intent(inout) :: this

  character(len=20), save :: sname = "end_ufield"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( associated( this%f1 ) ) deallocate( this%f1 )
  if ( associated( this%f2 ) ) deallocate( this%f2 )

  call write_dbg( cls_name, sname, cls_level, 'ends' )
  
end subroutine end_ufield

subroutine write_hdf5_single( this, file, dim )

  implicit none

  class( ufield ), intent(inout) :: this
  class( hdf5file ), intent(in) :: file
  integer, intent(in) :: dim

  ! local data
  integer :: ierr
  integer, dimension(2) :: lsize, gsize
  integer :: noff
  character(len=32), save :: sname = 'write_hdf5_single'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  gsize = this%nd
  lsize = this%ndp
  noff = this%noff(1)

  call pwfield( this%pp, file, this%f2(dim,1:,1:), gsize, lsize, noff, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_single

subroutine write_hdf5_pipe( this, file, dim, rtag, stag, id )

  implicit none

  class( ufield ), intent(inout) :: this
  class( hdf5file ), intent(in) :: file
  integer, intent(in) :: dim, rtag, stag
  integer, intent(inout) :: id

  ! local data
  integer :: nstage, stageid, ierr, start_pos
  integer, dimension(2) :: lsize, gsize, noff
  character(len=32), save :: sname = 'write_hdf5_pipe'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nstage = this%pp%getnstage()
  stageid = this%pp%getstageid()

  gsize = this%nd
  lsize = this%ndp
  noff = this%noff

  if ( nstage == 1 ) then
    call pwfield_pipe( this%pp, file, this%f2(dim,1:,1:), gsize, lsize, noff, &
      rtag, stag, id, ierr )
  else
    if ( stageid == 0 ) then
      lsize(2) = lsize(2) + 1
      start_pos = 1
    elseif ( stageid == nstage-1 ) then
      lsize(2) = lsize(2) - 1
      noff(2) = noff(2) + 1
      start_pos = 2
    else
      noff(2) = noff(2) + 1
      start_pos = 2
    endif
    call pwfield_pipe( this%pp, file, this%f2(dim,1:,start_pos:), gsize, lsize, noff, &
      rtag, stag, id, ierr )
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_pipe

subroutine copy_slice( this, idx, dir )

  implicit none

  class( ufield ), intent(inout) :: this
  integer, intent(in) :: idx, dir

  integer :: i, j, lb, ub
  character(len=20), save :: sname = "copy_slice"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  lb = 1-this%gc_num(p_lower,1)
  ub = this%ndp(1)+this%gc_num(p_upper,1)

  if ( .not. this%has_2d ) then
    call write_err( "The field has no 2D layout." )
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

  call write_dbg( cls_name, sname, cls_level, 'ends' )
  
end subroutine copy_slice

subroutine copy_gc_f1( this )

  implicit none

  class( ufield ), intent(inout) :: this

  integer :: idproc, idproc_left, idproc_right, nvp, comm
  integer :: nrp, count, dtype
  integer :: tag = 1, msgid, ierr
  integer, dimension(MPI_STATUS_SIZE) :: stat


  idproc = this%pp%getlidproc()
  idproc_left =  idproc - 1
  idproc_right = idproc + 1
  nrp = this%ndp(1)
  nvp = this%pp%getlnvp()
  comm = this%pp%getlgrp()
  dtype = this%pp%getmreal()

  ! forward message passing
  if ( this%gc_num(p_lower,1) > 0 ) then
    count = this%dim * this%gc_num(p_lower,1)
    ! receiver
    if ( idproc > 0 ) then
      call MPI_IRECV( this%f1(1,1-this%gc_num(p_lower,1)), count, dtype, &
        idproc_left, tag, comm, msgid, ierr )
    endif
    ! sender
    if ( idproc < nvp-1 ) then
      call MPI_SEND( this%f1(1,nrp+1-this%gc_num(p_lower,1)), count, dtype, &
        idproc_right, tag, comm, ierr )
    endif
    ! wait receiving finish
    if ( idproc > 0 ) then
      call MPI_WAIT( msgid, stat, ierr )
    endif
  endif

  ! backward message passing
  if ( this%gc_num(p_upper,1) > 0 ) then
    count = this%dim * this%gc_num(p_upper,1)
    ! receiver
    if ( idproc < nvp-1 ) then
      call MPI_IRECV( this%f1(1,nrp+1), count, dtype, &
        idproc_right, tag, comm, msgid, ierr )
    endif
    ! sender
    if ( idproc > 0 ) then
      call MPI_SEND( this%f1(1,1), count, dtype, &
        idproc_left, tag, comm, ierr )
    endif
    ! wait receiving finish
    if ( idproc < nvp-1 ) then
      call MPI_WAIT( msgid, stat, ierr )
    endif
  endif

end subroutine copy_gc_f1

subroutine copy_gc_f2( this, dir )

  implicit none

  class( ufield ), intent(inout) :: this
  integer, intent(in) :: dir

  integer :: idproc, idproc_next, idproc_last, nstage, stageid, nvp, comm
  integer :: nzp, count, dtype
  integer :: tag = 1, msgid, ierr
  integer, dimension(MPI_STATUS_SIZE) :: stat

  nvp = this%nvp(1)
  nstage = this%pp%getnstage()
  stageid = this%pp%getstageid()
  idproc = this%pp%getidproc()
  idproc_last = idproc - nvp
  idproc_next = idproc + nvp
  nzp = this%ndp(2)
  comm = this%pp%getlworld()
  dtype = this%pp%getmreal()

  select case (dir)

  ! forward message passing
  case ( p_mpi_forward )
  
    if ( this%gc_num(p_lower,2) > 0 ) then
      count = this%dim * size(this%f2,2) * this%gc_num(p_lower,2)
      ! receiver
      if ( stageid > 0 ) then
        call MPI_IRECV( this%f2(1,1-this%gc_num(p_lower,1),1-this%gc_num(p_lower,2)), &
          count, dtype, idproc_last, tag, comm, msgid, ierr )
      endif
      ! sender
      if ( stageid < nstage-1 ) then
        call MPI_SEND( this%f2(1,1-this%gc_num(p_lower,1),nzp+1-this%gc_num(p_lower,2)), &
          count, dtype, idproc_next, tag, comm, ierr )
      endif
      ! wait receiving finish
      if ( stageid > 0 ) then
        call MPI_WAIT( msgid, stat, ierr )
      endif
    endif

  ! backward message passing
  case ( p_mpi_backward )

    if ( this%gc_num(p_upper,2) > 0 ) then
      count = this%dim * size(this%f2,2) * this%gc_num(p_upper,2)
      ! receiver
      if ( stageid < nstage-1 ) then
        call MPI_IRECV( this%f2(1,1-this%gc_num(p_lower,1),nzp+1), &
          count, dtype, idproc_next, tag, comm, msgid, ierr )
      endif
      ! sender
      if ( stageid > 0 ) then
        call MPI_SEND( this%f2(1,1-this%gc_num(p_lower,1),1), &
          count, dtype, idproc_last, tag, comm, ierr )
      endif
      ! wait receiving finish
      if ( stageid < nstage-1 ) then
        call MPI_WAIT( msgid, stat, ierr )
      endif
    endif

  end select

end subroutine copy_gc_f2

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

      call write_err( "invalid assignment type!" )
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
  ! call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )
  call a3%new(a1)

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
  ! call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )
  call a3%new(a1)

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
  ! call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )
  call a3%new(a1)

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
  ! call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )
  call a3%new(a1)

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
  ! call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )
  call a3%new(a1)

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
  ! call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )
  call a3%new(a1)

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
  ! call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )
  call a3%new(a1)

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
  ! call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )
  call a3%new(a1)

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
  ! call a3%new( a1%dim, a1%get_nd(), a1%get_nvp(), a1%get_gc_num() )
  call a3%new(a1)

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

function get_noff_all( this )

  implicit none

  class( ufield ), intent(in) :: this
  integer, dimension(2) :: get_noff_all

  get_noff_all = this%noff

end function get_noff_all

function get_noff_dim( this, dim )

  implicit none

  class( ufield ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_noff_dim

  get_noff_dim = this%noff(dim)
  
end function get_noff_dim

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