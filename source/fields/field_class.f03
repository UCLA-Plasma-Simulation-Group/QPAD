module field_class

use parallel_module
use options_class
use ufield_class
use ufield_smooth_class
use hdf5io_class
use param
use sysutil_module
use mpi

implicit none

private

interface add_f1
  module procedure add_f1_binary
  module procedure add_f1_binary_dim
  module procedure add_f1_unitary
  module procedure add_f1_unitary_dim
end interface

interface add_f2
  module procedure add_f2_binary
  module procedure add_f2_unitary
end interface

interface sub_f1
  module procedure sub_f1_binary
  module procedure sub_f1_binary_dim
  module procedure sub_f1_unitary
  module procedure sub_f1_unitary_dim
end interface

interface sub_f2
  module procedure sub_f2_binary
  module procedure sub_f2_unitary
end interface

interface dot_f1
  module procedure dot_f1_unitary
  module procedure dot_f1_unitary_dim
end interface

public :: field
public :: add_f1, add_f2, sub_f1, sub_f2, dot_f1

character(len=20), parameter :: cls_name = "field"
integer, parameter :: cls_level = 3

type :: field

  ! private

  class( ufield ), dimension(:), pointer :: rf_re => null()
  class( ufield ), dimension(:), pointer :: rf_im => null()
  type( ufield_smooth ) :: smooth

  real :: dr, dxi
  integer :: max_mode, dim
  integer :: entity
  real, dimension(:,:,:), allocatable :: psend_buf, precv_buf

  contains

  procedure :: alloc => alloc_field
  generic :: new => init_field, init_field_cp
  procedure :: del => end_field
  generic :: get_rf_re => get_rf_re_all, get_rf_re_mode
  generic :: get_rf_im => get_rf_im_all, get_rf_im_mode
  generic :: write_hdf5 => write_hdf5_single, write_hdf5_pipe
  procedure :: smooth_f1
  procedure :: copy_slice
  procedure :: get_dr, get_dxi, get_max_mode, get_dim
  procedure :: copy_gc_f1, copy_gc_f2
  procedure :: acopy_gc_f1, acopy_gc_f2
  generic :: pipe_send => pipe_send_f1, pipe_send_f2
  generic :: pipe_recv => pipe_recv_f1, pipe_recv_f2

  procedure, private :: init_field, init_field_cp, end_field
  procedure, private :: get_rf_re_all, get_rf_re_mode, get_rf_im_all, get_rf_im_mode
  procedure, private :: write_hdf5_single, write_hdf5_pipe
  procedure, private :: pipe_send_f1, pipe_recv_f1
  procedure, private :: pipe_send_f2, pipe_recv_f2

  generic :: assignment(=) => assign_f1
  generic :: as            => assign_f2

  procedure, private :: assign_f1, assign_f2

end type field


contains

! =====================================================================
! Class field implementation
! =====================================================================
subroutine alloc_field( this, max_mode )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: max_mode

  ! placeholder, do nothing

end subroutine alloc_field

subroutine init_field( this, opts, dim, max_mode, gc_num, &
  entity, smooth_type, smooth_order, has_2d )

  implicit none

  class( field ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: max_mode, dim
  integer, intent(in), dimension(2,2) :: gc_num
  integer, intent(in), optional :: entity, smooth_type, smooth_order
  logical, intent(in), optional :: has_2d

  integer :: i, entity_, smooth_type_, smooth_order_
  logical :: has_2d_
  integer, dimension(2,2) :: gc_num_new
  character(len=20), save :: sname = "init_field"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ! set values of optional arguments
  entity_       = p_entity_none
  smooth_type_  = p_smooth_none
  smooth_order_ = 0
  has_2d_       = .true.
  if ( present(entity) ) entity_ = entity
  if ( present(smooth_type) ) smooth_type_ = smooth_type
  if ( present(smooth_order) ) smooth_order_ = smooth_order
  if ( present(has_2d) ) has_2d_ = has_2d

  this%dim      = dim
  this%max_mode = max_mode
  this%dr       = opts%get_dr()
  this%dxi      = opts%get_dxi()
  this%entity   = entity_

  ! gc_num_new(:,1) = gc_num(:,1)
  ! gc_num_new(:,2) = gc_num(:,2)

  call this%smooth%new( smooth_type_, smooth_order_ )
  gc_num_new(1,1) = max( gc_num(1,1), smooth_order_ )
  gc_num_new(2,1) = max( gc_num(2,1), smooth_order_ )
  gc_num_new(:,2) = gc_num(:,2)

  ! if ( present(smooth_type) .and. present(smooth_order) ) then
  !   call this%smooth%new( smooth_type, smooth_order )
  !   gc_num_new(1,1) = max( gc_num(1,1), smooth_order )
  !   gc_num_new(2,1) = max( gc_num(2,1), smooth_order )
  !   gc_num_new(:,2) = gc_num(:,2)
  ! else
  !   call this%smooth%new( p_smooth_none, 0 )
  ! endif

  allocate( this%rf_re(0:max_mode) )
  allocate( this%rf_im(max_mode) )
  do i = 0, this%max_mode
    call this%rf_re(i)%new( opts, dim, i, gc_num_new, has_2d=has_2d_ )
    if (i==0) cycle
    call this%rf_im(i)%new( opts, dim, i, gc_num_new, has_2d=has_2d_ )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field

subroutine init_field_cp( this, that )

  implicit none

  class( field ), intent(inout) :: this
  class( field ), intent(in) :: that

  integer :: i

  this%dim       = that%get_dim()
  this%max_mode = that%get_max_mode()
  this%dr        = that%get_dr()
  this%dxi       = that%get_dxi()
  this%entity    = that%entity

  allocate( this%rf_re(0:this%max_mode) )
  allocate( this%rf_im(this%max_mode) )
  do i = 0, this%max_mode
    call this%rf_re(i)%new( that%get_rf_re(i) )
    if (i==0) cycle
    call this%rf_im(i)%new( that%get_rf_im(i) )
  enddo

end subroutine init_field_cp

subroutine end_field( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "end_field"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%rf_re(0)%del()
  do i = 1, this%max_mode
    call this%rf_re(i)%del()
    call this%rf_im(i)%del()
  enddo
  deallocate( this%rf_re, this%rf_im )

  call this%smooth%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field

subroutine copy_slice( this, idx, dir )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: idx, dir

  integer :: i
  character(len=20), save :: sname = "copy_slice"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'copy slices' )

  do i = 0, this%max_mode
    call this%rf_re(i)%copy_slice( idx, dir )
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_slice( idx, dir )
  enddo

  call stop_tprof( 'copy slices' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_slice

subroutine copy_gc_f1( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "copy_gc_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%rf_re(i)%copy_gc_f1()
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_gc_f1()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_f1

subroutine copy_gc_f2( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "copy_gc_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%rf_re(i)%copy_gc_f2()
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_gc_f2()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_f2

subroutine acopy_gc_f1( this, dir, ncell )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: dir
  integer, intent(in), optional :: ncell

  integer :: i, nc
  character(len=20), save :: sname = "acopy_gc_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( present(ncell) ) then
    nc = ncell
  else
    nc = 1
  endif

  do i = 0, this%max_mode
    call this%rf_re(i)%acopy_gc_f1( dir, nc )
    if ( i == 0 ) cycle
    call this%rf_im(i)%acopy_gc_f1( dir, nc )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine acopy_gc_f1

subroutine acopy_gc_f2( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "acopy_gc_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%rf_re(i)%acopy_gc_f2()
    if ( i == 0 ) cycle
    call this%rf_im(i)%acopy_gc_f2()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine acopy_gc_f2

! subroutine pipe_gc_send( this, tag, sid )

!   implicit none

!   class( field ), intent(inout) :: this
!   integer, intent(in) :: tag
!   integer, intent(inout) :: sid

!   integer :: i, j, k, m
!   integer :: idproc, idproc_des, stageid, nstage, lnvp, comm
!   integer :: n1p, nzp, count, ierr
!   integer, dimension(2) :: gc1, gc2
!   real, dimension(:,:,:,:), allocatable, save :: buf
!   character(len=20), save :: sname = "pipe_gc_send"

!   call write_dbg( cls_name, sname, cls_level, 'starts' )
!   call start_tprof( 'pipeline' )

!   lnvp       = num_procs_loc()
!   stageid    = id_stage()
!   nstage     = num_stages()
!   idproc     = id_proc()
!   idproc_des = idproc + lnvp
!   n1p        = size(this%rf_re(0)%f1, 2)
!   nzp        = this%rf_re(0)%get_ndp(2)
!   comm       = comm_world()

!   gc1 = this%rf_re(0)%get_gc_num(1)
!   gc2 = this%rf_re(0)%get_gc_num(2)

!   if ( gc2(p_upper) <= 0 ) call write_err( 'No guard cells for pipeline guard cell copy!' )

!   if ( stageid == nstage-1 ) then
!     sid = MPI_REQUEST_NULL
!     call stop_tprof( 'pipeline' )
!     call write_dbg( cls_name, sname, cls_level, 'ends' )
!     return
!   endif

!   if ( .not. allocated(buf) ) then
!     allocate( buf(this%dim, n1p, gc2(p_upper), 0:2*this%max_mode) )
!   endif

!   ! pack m=0 mode
!   do k = 1, gc2(p_upper)
!     do j = 1, n1p
!       do i = 1, this%dim
!         buf(i,j,k,0) = this%rf_re(0)%f2(i, j-gc1(p_lower), nzp+k)
!       enddo
!     enddo
!   enddo

!   ! pack m>0 modes
!   if (this%max_mode > 0) then
!     do m = 1, this%max_mode
!       do k = 1, gc2(p_upper)
!         do j = 1, n1p
!           do i = 1, this%dim
!             buf(i,j,k,2*m-1) = this%rf_re(m)%f2(i, j-gc1(p_lower), nzp+k)
!             buf(i,j,k,2*m  ) = this%rf_im(m)%f2(i, j-gc1(p_lower), nzp+k)
!           enddo
!         enddo
!       enddo
!     enddo
!   endif

!   count = size(buf)
!   call MPI_ISEND( buf, count, p_dtype_real, idproc_des, tag, comm, sid, ierr )
!   ! check for error
!   if (ierr /= 0) call write_err('MPI_ISEND failed.')

!   call stop_tprof( 'pipeline' )
!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine pipe_gc_send

! subroutine pipe_gc_recv( this, tag )

!   implicit none

!   class( field ), intent(inout) :: this
!   integer, intent(in) :: tag
!   ! integer, intent(inout) :: rid

!   integer :: i, j, k, m
!   integer :: idproc, idproc_src, stageid, lnvp, comm
!   integer :: n1p, count, ierr
!   integer, dimension(2) :: gc1, gc2
!   integer, dimension(MPI_STATUS_SIZE) :: stat
!   real, dimension(:,:,:,:), allocatable, save :: buf
!   character(len=20), save :: sname = "pipe_gc_recv"

!   call write_dbg( cls_name, sname, cls_level, 'starts' )
!   call start_tprof( 'pipeline' )

!   lnvp       = num_procs_loc()
!   stageid    = id_stage()
!   idproc     = id_proc()
!   idproc_src = idproc - lnvp
!   n1p        = size(this%rf_re(0)%f1, 2)
!   comm       = comm_world()

!   gc1 = this%rf_re(0)%get_gc_num(1)
!   gc2 = this%rf_re(0)%get_gc_num(2)

!   if ( gc2(p_upper) <= 0 ) call write_err( 'No guard cells for pipeline guard cell copy!' )

!   if ( stageid == 0 ) then
!     call stop_tprof( 'pipeline' )
!     call write_dbg( cls_name, sname, cls_level, 'ends' )
!     return
!   endif

!   if ( .not. allocated(buf) ) then
!     allocate( buf(this%dim, n1p, gc2(p_upper), 0:2*this%max_mode) )
!   endif

!   count = size(buf)
!   call MPI_RECV( buf, count, p_dtype_real, idproc_src, tag, comm, stat, ierr )
!   ! check for error
!   if (ierr /= 0) call write_err('MPI_RECV failed.')

!   ! unpack m=0 mode
!   do k = 1, gc2(p_upper)
!     do j = 1, n1p
!       do i = 1, this%dim
!         this%rf_re(0)%f2(i, j-gc1(p_lower), k) = &
!           this%rf_re(0)%f2(i, j-gc1(p_lower), k) + buf(i,j,k,0)
!       enddo
!     enddo
!   enddo

!   ! unpack m>0 modes
!   if (this%max_mode > 0) then
!     do m = 1, this%max_mode
!       do k = 1, gc2(p_upper)
!         do j = 1, n1p
!           do i = 1, this%dim
!             this%rf_re(m)%f2(i, j-gc1(p_lower), k) = &
!               this%rf_re(m)%f2(i, j-gc1(p_lower), k) + buf(i,j,k,2*m-1)
!             this%rf_im(m)%f2(i, j-gc1(p_lower), k) = &
!               this%rf_im(m)%f2(i, j-gc1(p_lower), k) + buf(i,j,k,2*m  )
!           enddo
!         enddo
!       enddo
!     enddo
!   endif

!   call stop_tprof( 'pipeline' )
!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine pipe_gc_recv

subroutine pipe_send_f2( this, stag, id, dir, pos )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: stag
  integer, intent(inout) :: id
  character(len=*), intent(in) :: dir, pos

  integer :: i, j, m, idproc_des, idx_send
  integer :: nzp, n1p, count, ierr
  integer, dimension(2) :: gc
  character(len=20), save :: sname = "pipe_send_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )
  
  nzp = this%rf_re(0)%get_ndp(2)
  n1p = size( this%rf_re(0)%f1, 2 )
  gc  = this%rf_re(0)%get_gc_num(1)

  select case ( trim(dir) )
  case ( 'forward' )

    if ( id_stage() == num_stages() - 1 ) then
      id = MPI_REQUEST_NULL
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    select case ( trim(pos) )
    case ( 'inner' )
      idx_send = nzp
    case ( 'guard' )
      idx_send = nzp + 1
    case default
      call write_err( 'Invalid cell type! Valid options are "inner" and "guard".' )
    end select

    idproc_des = id_proc() + num_procs_loc()

  case ( 'backward' )

    if ( id_stage() == 0 ) then
      id = MPI_REQUEST_NULL
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    select case ( trim(pos) )
    case ( 'inner' )
      idx_send = 1
    case ( 'guard' )
      idx_send = 0
    case default
      call write_err( 'Invalid cell type! Valid options are "inner" and "guard".' )
    end select

    idproc_des = id_proc() - num_procs_loc()

  case default
    call write_err( 'Invalid data transfer direction! Valid options are "forward"&
      & and "backward".' )
  end select

  if ( .not. allocated( this%psend_buf ) ) then
    allocate( this%psend_buf( this%dim, n1p, 2*this%max_mode+1 ) )
  endif

  ! copy m = 0 mode
  do j = 1, n1p
    do i = 1, this%dim
      this%psend_buf(i,j,1) = this%rf_re(0)%f2( i, j-gc(1), idx_send )
    enddo
  enddo

  ! copy m > 0 mode
  if ( this%max_mode > 0 ) then
    do m = 1, this%max_mode
      do j = 1, n1p
        do i = 1, this%dim
          this%psend_buf(i,j,2*m) = this%rf_re(m)%f2( i, j-gc(1), idx_send )
        enddo
      enddo
      do j = 1, n1p
        do i = 1, this%dim
          this%psend_buf(i,j,2*m+1) = this%rf_im(m)%f2( i, j-gc(1), idx_send )
        enddo
      enddo
    enddo
  endif

  count = size( this%psend_buf )
  call mpi_isend( this%psend_buf, count, p_dtype_real, idproc_des, stag, &
    comm_world(), id, ierr )

  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_ISEND failed.' )
  endif

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_send_f2

subroutine pipe_recv_f2( this, rtag, dir, pos, mode )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: rtag
  character(len=*), intent(in) :: dir, pos, mode

  integer :: i, j, m, idx_recv, idproc_src, n1p, nzp
  integer :: count, ierr
  integer, dimension(2) :: gc
  integer, dimension(MPI_STATUS_SIZE) :: stat
  character(len=20), save :: sname = "pipe_precv_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  nzp = this%rf_re(0)%get_ndp(2)
  n1p = size( this%rf_re(0)%f1, 2 )
  gc  = this%rf_re(0)%get_gc_num(1)

  select case ( trim(dir) )
  case ( 'forward' )
  
    if ( id_stage() == 0 ) then
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    select case ( trim(pos) )
    case ( 'inner' )
      idx_recv = 1
    case ( 'guard' )
      idx_recv = 0
    case default
      call write_err( 'Invalid cell type! Valid options are "inner" and "guard".' )
    end select

    idproc_src = id_proc() - num_procs_loc()

  case ( 'backward' )

    if ( id_stage() == num_stages() - 1 ) then
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    select case ( trim(pos) )
    case ( 'inner' )
      idx_recv = nzp
    case ( 'guard' )
      idx_recv = nzp + 1
    case default
      call write_err( 'Invalid cell type! Valid options are "inner" and "guard".' )
    end select

    idproc_src = id_proc() + num_procs_loc()

  case default
    call write_err( 'Invalid data transfer direction! Valid options are "forward"&
      & and "backward".' )
  end select

  if ( .not. allocated( this%precv_buf ) ) then
    allocate( this%precv_buf( this%dim, n1p, 2*this%max_mode+1 ) )
  endif

  count = size( this%precv_buf )
  call mpi_recv( this%precv_buf, count, p_dtype_real, idproc_src, rtag, comm_world(), stat, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_RECV failed.' )
  endif

  select case ( trim(mode) )
  case ( 'replace' )

    ! copy m=0 mode
    do j = 1, n1p
      do i = 1, this%dim
        this%rf_re(0)%f2( i, j-gc(1), idx_recv ) = this%precv_buf(i,j,1)
      enddo
    enddo

    ! copy m>0 mode
    if ( this%max_mode > 0 ) then
      do m = 1, this%max_mode
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_re(m)%f2( i, j-gc(1), idx_recv ) = this%precv_buf(i,j,2*m)
          enddo
        enddo
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_im(m)%f2( i, j-gc(1), idx_recv ) = this%precv_buf(i,j,2*m+1)
          enddo
        enddo
      enddo
    endif

  case ( 'add' )

    ! copy m=0 mode
    do j = 1, n1p
      do i = 1, this%dim
        this%rf_re(0)%f2( i, j-gc(1), idx_recv ) = &
        this%rf_re(0)%f2( i, j-gc(1), idx_recv ) + this%precv_buf(i,j,1)
      enddo
    enddo

    ! copy m>0 mode
    if ( this%max_mode > 0 ) then
      do m = 1, this%max_mode
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_re(m)%f2( i, j-gc(1), idx_recv ) = &
            this%rf_re(m)%f2( i, j-gc(1), idx_recv ) + this%precv_buf(i,j,2*m)
          enddo
        enddo
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_im(m)%f2( i, j-gc(1), idx_recv ) = &
            this%rf_im(m)%f2( i, j-gc(1), idx_recv ) + this%precv_buf(i,j,2*m+1)
          enddo
        enddo
      enddo
    endif

  case default
    call write_err( 'Invalid data transfer mode! Valid options are "replace" &
      &and "add".' )
  end select

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_recv_f2

subroutine pipe_send_f1( this, stag, id, dir )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: stag
  integer, intent(inout) :: id
  character(len=*), intent(in) :: dir

  integer :: i, j, m, idproc_des
  integer :: n1p, count, ierr
  integer, dimension(2) :: gc
  character(len=20), save :: sname = "pipe_send_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )
  
  n1p = size( this%rf_re(0)%f1, 2 )
  gc  = this%rf_re(0)%get_gc_num(1)

  select case ( trim(dir) )
  case ( 'forward' )

    if ( id_stage() == num_stages() - 1 ) then
      id = MPI_REQUEST_NULL
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    idproc_des = id_proc() + num_procs_loc()

  case ( 'backward' )

    if ( id_stage() == 0 ) then
      id = MPI_REQUEST_NULL
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    idproc_des = id_proc() - num_procs_loc()

  case default
    call write_err( 'Invalid data transfer direction! Valid options are "forward"&
      & and "backward".' )
  end select

  if ( .not. allocated( this%psend_buf ) ) then
    allocate( this%psend_buf( this%dim, n1p, 2*this%max_mode+1 ) )
  endif

  ! copy m = 0 mode
  do j = 1, n1p
    do i = 1, this%dim
      this%psend_buf(i,j,1) = this%rf_re(0)%f1( i, j-gc(1) )
    enddo
  enddo

  ! copy m > 0 mode
  if ( this%max_mode > 0 ) then
    do m = 1, this%max_mode
      do j = 1, n1p
        do i = 1, this%dim
          this%psend_buf(i,j,2*m) = this%rf_re(m)%f1( i, j-gc(1) )
        enddo
      enddo
      do j = 1, n1p
        do i = 1, this%dim
          this%psend_buf(i,j,2*m+1) = this%rf_im(m)%f1( i, j-gc(1) )
        enddo
      enddo
    enddo
  endif

  count = size( this%psend_buf )
  call mpi_isend( this%psend_buf, count, p_dtype_real, idproc_des, stag, comm_world(), id, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_ISEND failed.' )
  endif

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_send_f1

subroutine pipe_recv_f1( this, rtag, dir, mode )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: rtag
  character(len=*), intent(in) :: dir, mode

  integer :: i, j, m, idproc_src, n1p
  integer :: count, ierr
  integer, dimension(2) :: gc
  integer, dimension(MPI_STATUS_SIZE) :: stat
  character(len=20), save :: sname = "pipe_precv_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  n1p = size( this%rf_re(0)%f1, 2 )
  gc  = this%rf_re(0)%get_gc_num(1)

  select case ( trim(dir) )
  case ( 'forward' )
  
    if ( id_stage() == 0 ) then
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    idproc_src = id_proc() - num_procs_loc()

  case ( 'backward' )

    if ( id_stage() == num_stages() - 1 ) then
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    idproc_src = id_proc() + num_procs_loc()

  case default
    call write_err( 'Invalid data transfer direction! Valid options are "forward"&
      & and "backward".' )
  end select

  if ( .not. allocated( this%precv_buf ) ) then
    allocate( this%precv_buf( this%dim, n1p, 2*this%max_mode+1 ) )
  endif

  count = size( this%precv_buf )
  call mpi_recv( this%precv_buf, count, p_dtype_real, idproc_src, rtag, comm_world(), stat, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_RECV failed.' )
  endif

  select case ( trim(mode) )
  case ( 'replace' )

    ! copy m=0 mode
    do j = 1, n1p
      do i = 1, this%dim
        this%rf_re(0)%f1( i, j-gc(1) ) = this%precv_buf(i,j,1)
      enddo
    enddo

    ! copy m>0 mode
    if ( this%max_mode > 0 ) then
      do m = 1, this%max_mode
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_re(m)%f1( i, j-gc(1) ) = this%precv_buf(i,j,2*m)
          enddo
        enddo
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_im(m)%f1( i, j-gc(1) ) = this%precv_buf(i,j,2*m+1)
          enddo
        enddo
      enddo
    endif

  case ( 'add' )

    ! copy m=0 mode
    do j = 1, n1p
      do i = 1, this%dim
        this%rf_re(0)%f1( i, j-gc(1) ) = &
        this%rf_re(0)%f1( i, j-gc(1) ) + this%precv_buf(i,j,1)
      enddo
    enddo

    ! copy m>0 mode
    if ( this%max_mode > 0 ) then
      do m = 1, this%max_mode
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_re(m)%f1( i, j-gc(1) ) = &
            this%rf_re(m)%f1( i, j-gc(1) ) + this%precv_buf(i,j,2*m)
          enddo
        enddo
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_im(m)%f1( i, j-gc(1) ) = &
            this%rf_im(m)%f1( i, j-gc(1) ) + this%precv_buf(i,j,2*m+1)
          enddo
        enddo
      enddo
    endif

  case default
    call write_err( 'Invalid data transfer mode! Valid options are "replace" &
      &and "add".' )
  end select

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_recv_f1

subroutine write_hdf5_single( this, files, dim )

  implicit none

  class( field ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim

  integer :: i
  character(len=32), save :: sname = 'write_hdf5_single'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode

    if ( i == 0 ) then
      call this%rf_re(i)%write_hdf5( files(1), dim )
      cycle
    endif

    call this%rf_re(i)%write_hdf5( files(2*i), dim )
    call this%rf_im(i)%write_hdf5( files(2*i+1), dim )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_single

subroutine write_hdf5_pipe( this, files, dim, rtag, stag, id )

  implicit none

  class( field ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim, rtag, stag
  integer, intent(inout) :: id

  integer :: i
  character(len=32), save :: sname = 'write_hdf5_pipe'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode

    if ( i == 0 ) then
      call this%rf_re(i)%write_hdf5( files(1), dim, rtag, stag, id )
      cycle
    endif

    call this%rf_re(i)%write_hdf5( files(2*i), dim, rtag, stag, id )
    call this%rf_im(i)%write_hdf5( files(2*i+1), dim, rtag, stag, id )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_pipe

subroutine smooth_f1( this )

  implicit none

  class( field ), intent(inout) :: this
  ! logical, intent(in) :: q_cons

  integer :: i
  character(len=32), save :: sname = 'smooth_f1'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. this%smooth%if_smooth() ) return

  do i = 0, this%max_mode

    if ( i == 0 ) then
      call this%smooth%smooth_f1( this%rf_re(i) )
      cycle
    endif

    call this%smooth%smooth_f1( this%rf_re(i) )
    call this%smooth%smooth_f1( this%rf_im(i) )

  enddo

  call this%acopy_gc_f1( dir=p_mpi_bothway, ncell=this%smooth%get_order() )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine smooth_f1

function get_dr( this )

  implicit none

  class( field ), intent(in) :: this
  real :: get_dr

  get_dr = this%dr

end function get_dr

function get_dxi( this )

  implicit none

  class( field ), intent(in) :: this
  real :: get_dxi

  get_dxi = this%dxi

end function get_dxi

function get_max_mode( this )

  implicit none

  class( field ), intent(in) :: this
  integer :: get_max_mode

  get_max_mode = this%max_mode

end function get_max_mode

function get_dim( this )

  implicit none

  class( field ), intent(in) :: this
  integer :: get_dim

  get_dim = this%dim

end function get_dim

function get_rf_re_all( this )

  implicit none

  class( field ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_rf_re_all

  get_rf_re_all => this%rf_re

end function get_rf_re_all

function get_rf_re_mode( this, mode )

  implicit none

  class( field ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_rf_re_mode

  get_rf_re_mode => this%rf_re(mode)

end function get_rf_re_mode

function get_rf_im_all( this )

  implicit none

  class( field ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_rf_im_all

  get_rf_im_all => this%rf_im

end function get_rf_im_all

function get_rf_im_mode( this, mode )

  implicit none

  class( field ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_rf_im_mode

  get_rf_im_mode => this%rf_im(mode)

end function get_rf_im_mode

subroutine assign_f1( this, that )

  implicit none

  class( field ), intent(inout) :: this
  class(*), intent(in) :: that

  integer :: i

  select type (that)

  type is (real)

    do i = 0, this%max_mode
      this%rf_re(i) = that
      if (i == 0) cycle
      this%rf_im(i) = that
    enddo

  class is (field)

    do i = 0, this%max_mode
      this%rf_re(i) = that%rf_re(i)
      if (i == 0) cycle
      this%rf_im(i) = that%rf_im(i)
    enddo

  class default

    call write_err( "invalid assignment type!" )

  end select

end subroutine assign_f1

subroutine assign_f2( this, that )

  implicit none

  class( field ), intent(inout) :: this
  class(*), intent(in) :: that

  integer :: i

  select type (that)

  type is (real)

    do i = 0, this%max_mode
      call this%rf_re(i)%as( that )
      if (i == 0) cycle
      call this%rf_im(i)%as( that )
    enddo

  class is (field)

    do i = 0, this%max_mode
      call this%rf_re(i)%as( that%get_rf_re(i) )
      if (i == 0) cycle
      call this%rf_im(i)%as( that%get_rf_im(i) )
    enddo

  class default

    call write_err( "invalid assignment type!" )

  end select

end subroutine assign_f2

subroutine add_f1_binary_dim( a1, a2, a3, dim1, dim2, dim3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3
  integer, intent(in), dimension(:) :: dim1, dim2, dim3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%max_mode
    call add_f1( a1%rf_re(i), a2%rf_re(i), ua3_re(i), dim1, dim2, dim3 )
    if (i==0) cycle
    call add_f1( a1%rf_im(i), a2%rf_im(i), ua3_im(i), dim1, dim2, dim3 )
  enddo

end subroutine add_f1_binary_dim

subroutine add_f1_unitary_dim( a1, a2, dim1, dim2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) :: a2
  integer, intent(in), dimension(:) :: dim1, dim2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%max_mode
    call add_f1( a1%rf_re(i), ua2_re(i), dim1, dim2 )
    if (i==0) cycle
    call add_f1( a1%rf_im(i), ua2_im(i), dim1, dim2 )
  enddo

end subroutine add_f1_unitary_dim

subroutine add_f1_binary( a1, a2, a3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%max_mode
    call add_f1( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
    if (i==0) cycle
    call add_f1( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
  enddo

end subroutine add_f1_binary

subroutine add_f1_unitary( a1, a2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) :: a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%max_mode
    call add_f1( a1%rf_re(i), ua2_re(i) )
    if (i==0) cycle
    call add_f1( a1%rf_im(i), ua2_im(i) )
  enddo

end subroutine add_f1_unitary

subroutine sub_f1_binary_dim( a1, a2, a3, dim1, dim2, dim3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3
  integer, intent(in), dimension(:) :: dim1, dim2, dim3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%max_mode
    call sub_f1( a1%rf_re(i), a2%rf_re(i), ua3_re(i), dim1, dim2, dim3 )
    if (i==0) cycle
    call sub_f1( a1%rf_im(i), a2%rf_im(i), ua3_im(i), dim1, dim2, dim3 )
  enddo

end subroutine sub_f1_binary_dim

subroutine sub_f1_unitary_dim( a1, a2, dim1, dim2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) :: a2
  integer, intent(in), dimension(:) :: dim1, dim2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%max_mode
    call sub_f1( a1%rf_re(i), ua2_re(i), dim1, dim2 )
    if (i==0) cycle
    call sub_f1( a1%rf_im(i), ua2_im(i), dim1, dim2 )
  enddo

end subroutine sub_f1_unitary_dim

subroutine sub_f1_binary( a1, a2, a3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%max_mode
    call sub_f1( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
    if (i==0) cycle
    call sub_f1( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
  enddo

end subroutine sub_f1_binary

subroutine sub_f1_unitary( a1, a2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) :: a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%max_mode
    call sub_f1( a1%rf_re(i), ua2_re(i) )
    if (i==0) cycle
    call sub_f1( a1%rf_im(i), ua2_im(i) )
  enddo

end subroutine sub_f1_unitary

subroutine dot_f1_unitary( a1, a2 )

  implicit none

  real, intent(in) :: a1
  class( field ), intent(inout) :: a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a2%max_mode
    call dot_f1( a1, ua2_re(i) )
    if (i==0) cycle
    call dot_f1( a1, ua2_im(i) )
  enddo

end subroutine dot_f1_unitary

subroutine dot_f1_unitary_dim( a1, a2, dim )

  implicit none

  real, intent(in) :: a1
  class( field ), intent(inout) :: a2
  integer, intent(in), dimension(:) :: dim

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a2%max_mode
    call dot_f1( a1, ua2_re(i), dim )
    if (i==0) cycle
    call dot_f1( a1, ua2_im(i), dim )
  enddo

end subroutine dot_f1_unitary_dim

subroutine add_f2_binary( a1, a2, a3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%max_mode
    call add_f2( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
    if (i==0) cycle
    call add_f2( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
  enddo

end subroutine add_f2_binary

subroutine add_f2_unitary( a1, a2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) ::a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%max_mode
    call add_f2( a1%rf_re(i), ua2_re(i) )
    if (i==0) cycle
    call add_f2( a1%rf_im(i), ua2_im(i) )
  enddo

end subroutine add_f2_unitary

subroutine sub_f2_binary( a1, a2, a3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%max_mode
    call sub_f2( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
    if (i==0) cycle
    call sub_f2( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
  enddo

end subroutine sub_f2_binary

subroutine sub_f2_unitary( a1, a2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) ::a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%max_mode
    call sub_f2( a1%rf_re(i), ua2_re(i) )
    if (i==0) cycle
    call sub_f2( a1%rf_im(i), ua2_im(i) )
  enddo

end subroutine sub_f2_unitary

end module field_class