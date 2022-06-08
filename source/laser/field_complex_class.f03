module field_complex_class

use parallel_module
use input_class
use options_class
use ufield_class
use ufield_smooth_class
use hdf5io_class
use param
use sysutil_module
use mpi
use kwargs_class
use ppmsg_class

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

public :: field_complex
public :: add_f1, add_f2, sub_f1, sub_f2, dot_f1

character(len=20), parameter :: cls_name = "field_complex"
integer, parameter :: cls_level = 3

type :: field_complex

  private

  class( ufield ), dimension(:), pointer, public :: cfr_re => null()
  class( ufield ), dimension(:), pointer, public :: cfr_im => null()
  class( ufield ), dimension(:), pointer, public :: cfi_re => null()
  class( ufield ), dimension(:), pointer, public :: cfi_im => null()

  real, public :: dr, dz
  integer, public :: max_mode, dim
  integer, dimension(2,2), public :: gc_num
  real, dimension(:), allocatable :: psend_buf, precv_buf

  contains

  procedure :: alloc => alloc_field_complex
  procedure :: new => init_field_complex
  procedure :: del => end_field_complex
  generic :: get_cfr_re => get_cfr_re_all, get_cfr_re_mode
  generic :: get_cfr_im => get_cfr_im_all, get_cfr_im_mode
  generic :: get_cfi_re => get_cfi_re_all, get_cfi_re_mode
  generic :: get_cfi_im => get_cfi_im_all, get_cfi_im_mode
  generic :: write_hdf5 => write_hdf5_single, write_hdf5_pipe
  procedure :: read_hdf5
  procedure :: copy_slice
  procedure :: copy_gc_f1
  procedure :: copy_gc_f2
  procedure :: acopy_gc_f1
  procedure :: acopy_gc_f2
  generic :: pipe_send => pipe_send_f1, pipe_send_f2
  generic :: pipe_recv => pipe_recv_f1, pipe_recv_f2

  procedure, private :: init_field_complex
  procedure, private :: end_field_complex
  procedure, private :: get_cfr_re_all, get_cfr_re_mode
  procedure, private :: get_cfr_im_all, get_cfr_im_mode
  procedure, private :: get_cfi_re_all, get_cfi_re_mode
  procedure, private :: get_cfi_im_all, get_cfi_im_mode
  procedure, private :: write_hdf5_single, write_hdf5_pipe
  procedure, private :: pipe_send_f1, pipe_recv_f1
  procedure, private :: pipe_send_f2, pipe_recv_f2
  procedure, private :: check_send_buf
  procedure, private :: check_recv_buf

  generic :: assignment(=) => assign_f1
  generic :: as            => assign_f2

  procedure, private :: assign_f1, assign_f2

end type field_complex

contains

! =====================================================================
! Class field implementation
! =====================================================================
subroutine alloc_field_complex( this, input, opts, id )

  implicit none

  class( field_complex ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  integer, intent(in) :: id

  ! placeholder, do nothing

end subroutine alloc_field_complex

subroutine init_field_complex( this, opts, dim, max_mode, gc_num, only_f1, kwargs )

  implicit none
  class( field_complex ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: max_mode, dim
  integer, intent(in), dimension(2,2) :: gc_num
  logical, intent(in), optional :: only_f1
  type( kw_list ), intent(in), optional :: kwargs

  integer :: i
  logical :: only_f1_tmp
  character(len=32), save :: sname = "init_field_complex"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ! set values of optional arguments
  only_f1_tmp = .false.
  if ( present(only_f1) ) only_f1_tmp = only_f1

  this%dim      = dim
  this%max_mode = max_mode
  this%dr       = opts%get_dr()
  this%dz       = opts%get_dxi()
  this%gc_num   = gc_num

  allocate( this%cfr_re(0:max_mode) )
  allocate( this%cfi_re(0:max_mode) )
  allocate( this%cfr_im(max_mode) )
  allocate( this%cfi_im(max_mode) )
  do i = 0, this%max_mode
    call this%cfr_re(i)%new( opts, dim, i, gc_num, has_2d=.not. only_f1_tmp )
    call this%cfi_re(i)%new( opts, dim, i, gc_num, has_2d=.not. only_f1_tmp )
    if (i==0) cycle
    call this%cfr_im(i)%new( opts, dim, i, gc_num, has_2d=.not. only_f1_tmp )
    call this%cfi_im(i)%new( opts, dim, i, gc_num, has_2d=.not. only_f1_tmp )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_complex

subroutine end_field_complex( this )

  implicit none
  class( field_complex ), intent(inout) :: this

  integer :: i
  character(len=32), save :: sname = "end_field_complex"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%cfr_re(i)%del()
    call this%cfi_re(i)%del()
    if ( i == 0 ) cycle
    call this%cfr_im(i)%del()
    call this%cfi_im(i)%del()
  enddo
  deallocate( this%cfr_re, this%cfr_im, this%cfi_re, this%cfi_im )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_complex

subroutine copy_slice( this, idx, dir )

  implicit none

  class( field_complex ), intent(inout) :: this
  integer, intent(in) :: idx, dir

  integer :: i
  character(len=32), save :: sname = "copy_slice"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'copy slices' )

  do i = 0, this%max_mode
    call this%cfr_re(i)%copy_slice( idx, dir )
    call this%cfi_re(i)%copy_slice( idx, dir )
    if ( i == 0 ) cycle
    call this%cfr_im(i)%copy_slice( idx, dir )
    call this%cfi_im(i)%copy_slice( idx, dir )
  enddo

  call stop_tprof( 'copy slices' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_slice

subroutine copy_gc_f1( this )

  implicit none

  class( field_complex ), intent(inout) :: this

  integer :: i
  character(len=32), save :: sname = "copy_gc_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%cfr_re(i)%copy_gc_f1()
    call this%cfi_re(i)%copy_gc_f1()
    if ( i == 0 ) cycle
    call this%cfr_im(i)%copy_gc_f1()
    call this%cfi_im(i)%copy_gc_f1()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_f1

subroutine copy_gc_f2( this )

  implicit none

  class( field_complex ), intent(inout) :: this

  integer :: i
  character(len=32), save :: sname = "copy_gc_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%cfr_re(i)%copy_gc_f2()
    call this%cfi_re(i)%copy_gc_f2()
    if ( i == 0 ) cycle
    call this%cfr_im(i)%copy_gc_f2()
    call this%cfi_im(i)%copy_gc_f2()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_f2

subroutine acopy_gc_f1( this, dir, ncell )

  implicit none

  class( field_complex ), intent(inout) :: this
  integer, intent(in) :: dir
  integer, intent(in), optional :: ncell

  integer :: i, nc
  character(len=32), save :: sname = "acopy_gc_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nc = 1; if ( present(ncell) ) nc = ncell

  do i = 0, this%max_mode
    call this%cfr_re(i)%acopy_gc_f1( dir, nc )
    call this%cfi_re(i)%acopy_gc_f1( dir, nc )
    if ( i == 0 ) cycle
    call this%cfr_im(i)%acopy_gc_f1( dir, nc )
    call this%cfi_im(i)%acopy_gc_f1( dir, nc )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine acopy_gc_f1

subroutine acopy_gc_f2( this )

  implicit none

  class( field_complex ), intent(inout) :: this

  integer :: i
  character(len=32), save :: sname = "acopy_gc_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%cfr_re(i)%acopy_gc_f2()
    call this%cfi_re(i)%acopy_gc_f2()
    if ( i == 0 ) cycle
    call this%cfr_im(i)%acopy_gc_f2()
    call this%cfi_im(i)%acopy_gc_f2()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine acopy_gc_f2

subroutine pipe_send_f2( this, pp_msg, dir, pos_type, num_slices )

  implicit none

  class( field_complex ), intent(inout) :: this
  type( ppmsg ), intent(inout) :: pp_msg
  character(len=*), intent(in) :: dir, pos_type
  integer, intent(in) :: num_slices

  integer :: i, j, k, m, idproc_des, idx_send, stride1, stride2, stride3
  integer :: nzp, n1p, count, ierr, buf_idx, offset
  integer, dimension(2) :: gc
  character(len=20), save :: sname = "pipe_send_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )
  
  nzp = this%cfr_re(0)%get_ndp(2)
  n1p = size( this%cfr_re(0)%f1, 2 )
  gc  = this%gc_num(:,1)

  select case ( trim(dir) )
    case ( 'forward' )

      if ( id_stage() == num_stages() - 1 ) then
        call pp_msg%set_idle()
        call stop_tprof( 'pipeline' )
        call write_dbg( cls_name, sname, cls_level, 'ends' )
        return
      endif

      select case ( trim(pos_type) )
      case ( 'inner' )
        idx_send = nzp - num_slices + 1
      case ( 'guard' )
        idx_send = nzp + 1
      case default
        call write_err( 'Invalid cell type! Valid options are "inner" and "guard".' )
      end select

      idproc_des = id_proc() + num_procs_loc()

    case ( 'backward' )

      if ( id_stage() == 0 ) then
        call pp_msg%set_idle()
        call stop_tprof( 'pipeline' )
        call write_dbg( cls_name, sname, cls_level, 'ends' )
        return
      endif

      select case ( trim(pos_type) )
      case ( 'inner' )
        idx_send = 1
      case ( 'guard' )
        idx_send = 1 - num_slices
      case default
        call write_err( 'Invalid cell type! Valid options are "inner" and "guard".' )
      end select

      idproc_des = id_proc() - num_procs_loc()

    case default
      call write_err( 'Invalid data transfer direction! Valid options are "forward"&
        & and "backward".' )
  end select

  count = this%dim * n1p * num_slices * ( 4 * this%max_mode + 2 )
  call this%check_send_buf( count )

  ! copy m = 0 mode
  do k = 1, num_slices
    stride1 = 2 * this%dim * n1p * (k-1)
    do j = 1, n1p
      stride2 = 2 * this%dim * (j-1)
      do i = 1, this%dim
        buf_idx = (2*i-1) + stride1 + stride2
        this%psend_buf( buf_idx   ) = this%cfr_re(0)%f2( i, j-gc(1), idx_send+k-1 )
        this%psend_buf( buf_idx+1 ) = this%cfi_re(0)%f2( i, j-gc(1), idx_send+k-1 )
      enddo
    enddo
  enddo

  ! copy m > 0 mode
  offset = 2 * this%dim * n1p * num_slices
  do m = 1, this%max_mode
    stride1 = 4 * this%dim * n1p * num_slices * (m-1)
    do k = 1, num_slices
      stride2 = 4 * this%dim * n1p * (k-1)
      do j = 1, n1p
        stride3 = 4 * this%dim * (j-1)
        do i = 1, this%dim
          buf_idx = (4*i-3) + stride1 + stride2 + stride3 + offset
          this%psend_buf( buf_idx   ) = this%cfr_re(m)%f2( i, j-gc(1), idx_send+k-1 )
          this%psend_buf( buf_idx+1 ) = this%cfr_im(m)%f2( i, j-gc(1), idx_send+k-1 )
          this%psend_buf( buf_idx+2 ) = this%cfi_re(m)%f2( i, j-gc(1), idx_send+k-1 )
          this%psend_buf( buf_idx+3 ) = this%cfi_im(m)%f2( i, j-gc(1), idx_send+k-1 )
        enddo
      enddo
    enddo
  enddo

  call pp_msg%wait_task()
  call mpi_isend( this%psend_buf, count, p_dtype_real, idproc_des, pp_msg%tag, &
    comm_world(), pp_msg%id, ierr )

  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_ISEND failed.' )
  endif

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_send_f2

subroutine pipe_recv_f2( this, pp_msg, dir, pos_type, mode, num_slices )

  implicit none

  class( field_complex ), intent(inout) :: this
  type( ppmsg ), intent(inout) :: pp_msg
  character(len=*), intent(in) :: dir, pos_type, mode
  integer, intent(in) :: num_slices

  integer :: i, j, k, m, idx_recv, idproc_src, n1p, nzp, stride1, stride2, stride3
  integer :: count, ierr, buf_idx, offset
  integer, dimension(2) :: gc
  integer, dimension(MPI_STATUS_SIZE) :: stat
  character(len=32), save :: sname = "pipe_precv_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  nzp = this%cfr_re(0)%get_ndp(2)
  n1p = size( this%cfr_re(0)%f1, 2 )
  gc  = this%gc_num(:,1)

  select case ( trim(dir) )
  case ( 'forward' )
  
    if ( id_stage() == 0 ) then
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    select case ( trim(pos_type) )
    case ( 'inner' )
      idx_recv = 1
    case ( 'guard' )
      idx_recv = 1 - num_slices
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

    select case ( trim(pos_type) )
    case ( 'inner' )
      idx_recv = nzp - num_slices + 1
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

  count = this%dim * n1p * num_slices * ( 4 * this%max_mode + 2 )
  call this%check_recv_buf( count )

  call mpi_recv( this%precv_buf, count, p_dtype_real, idproc_src, pp_msg%tag, comm_world(), stat, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_RECV failed.' )
  endif

  select case ( trim(mode) )
  case ( 'replace' )

    ! copy m = 0 mode
    do k = 1, num_slices
      stride1 = 2 * this%dim * n1p * (k-1)
      do j = 1, n1p
        stride2 = 2 * this%dim * (j-1)
        do i = 1, this%dim
          buf_idx = (2*i-1) + stride1 + stride2
          this%cfr_re(0)%f2( i, j-gc(1), idx_recv+k-1 ) = this%precv_buf( buf_idx   )
          this%cfi_re(0)%f2( i, j-gc(1), idx_recv+k-1 ) = this%precv_buf( buf_idx+1 )
        enddo
      enddo
    enddo

    ! copy m > 0 modes
    offset = 2 * this%dim * n1p * num_slices
    do m = 1, this%max_mode
      stride1 = 4 * this%dim * n1p * num_slices * (m-1)
      do k = 1, num_slices
        stride2 = 4 * this%dim * n1p * (k-1)
        do j = 1, n1p
          stride3 = 4 * this%dim * (j-1)
          do i = 1, this%dim
            buf_idx = (4*i-3) + stride1 + stride2 + stride3 + offset
            this%cfr_re(m)%f2( i, j-gc(1), idx_recv+k-1 ) = this%precv_buf( buf_idx   )
            this%cfr_im(m)%f2( i, j-gc(1), idx_recv+k-1 ) = this%precv_buf( buf_idx+1 )
            this%cfi_re(m)%f2( i, j-gc(1), idx_recv+k-1 ) = this%precv_buf( buf_idx+2 )
            this%cfi_im(m)%f2( i, j-gc(1), idx_recv+k-1 ) = this%precv_buf( buf_idx+3 )
          enddo
        enddo
      enddo
    enddo

  case ( 'add' )

    ! copy m = 0 mode
    do k = 1, num_slices
      stride1 = 2 * this%dim * n1p * (k-1)
      do j = 1, n1p
        stride2 = 2 * this%dim * (j-1)
        do i = 1, this%dim
          buf_idx = (2*i-1) + stride1 + stride2
          this%cfr_re(0)%f2( i, j-gc(1), idx_recv+k-1 ) = &
          this%cfr_re(0)%f2( i, j-gc(1), idx_recv+k-1 ) + this%precv_buf( buf_idx   )
          this%cfi_re(0)%f2( i, j-gc(1), idx_recv+k-1 ) = &
          this%cfi_re(0)%f2( i, j-gc(1), idx_recv+k-1 ) + this%precv_buf( buf_idx+1 )
        enddo
      enddo
    enddo

    ! copy m > 0 modes
    offset = 2 * this%dim * n1p * num_slices
    do m = 1, this%max_mode
      stride1 = 4 * this%dim * n1p * num_slices * (m-1)
      do k = 1, num_slices
        stride2 = 4 * this%dim * n1p * (k-1)
        do j = 1, n1p
          stride3 = 4 * this%dim * (j-1)
          do i = 1, this%dim
            buf_idx = (4*i-3) + stride1 + stride2 + stride3 + offset
            this%cfr_re(m)%f2( i, j-gc(1), idx_recv+k-1 ) = &
            this%cfr_re(m)%f2( i, j-gc(1), idx_recv+k-1 ) + this%precv_buf( buf_idx   )
            this%cfr_im(m)%f2( i, j-gc(1), idx_recv+k-1 ) = &
            this%cfr_im(m)%f2( i, j-gc(1), idx_recv+k-1 ) + this%precv_buf( buf_idx+1 )
            this%cfi_re(m)%f2( i, j-gc(1), idx_recv+k-1 ) = &
            this%cfi_re(m)%f2( i, j-gc(1), idx_recv+k-1 ) + this%precv_buf( buf_idx+2 )
            this%cfi_im(m)%f2( i, j-gc(1), idx_recv+k-1 ) = &
            this%cfi_im(m)%f2( i, j-gc(1), idx_recv+k-1 ) + this%precv_buf( buf_idx+3 )
          enddo
        enddo
      enddo
    enddo

  case default
    call write_err( 'Invalid data transfer mode! Valid options are "replace" &
      &and "add".' )
  end select

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_recv_f2

subroutine pipe_send_f1( this, pp_msg, dir )

  implicit none
  class( field_complex ), intent(inout) :: this
  type( ppmsg ), intent(inout) :: pp_msg
  character(len=*), intent(in) :: dir

  integer :: i, j, m, idproc_des, buf_idx, stride1, stride2, offset
  integer :: n1p, count, ierr
  integer, dimension(2) :: gc
  character(len=32), save :: sname = "pipe_send_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )
  
  n1p = size( this%cfr_re(0)%f1, 2 )
  gc  = this%gc_num(:,1)

  select case ( trim(dir) )
  case ( 'forward' )

    if ( id_stage() == num_stages() - 1 ) then
      call pp_msg%set_idle()
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    idproc_des = id_proc() + num_procs_loc()

  case ( 'backward' )

    if ( id_stage() == 0 ) then
      call pp_msg%set_idle
      call stop_tprof( 'pipeline' )
      call write_dbg( cls_name, sname, cls_level, 'ends' )
      return
    endif

    idproc_des = id_proc() - num_procs_loc()

  case default
    call write_err( 'Invalid data transfer direction! Valid options are "forward"&
      & and "backward".' )
  end select

  count = this%dim * n1p * ( 4 * this%max_mode + 2 )
  call this%check_send_buf( count )

  ! copy m = 0 mode
  do j = 1, n1p
    stride1 = 2 * this%dim * (j-1)
    do i = 1, this%dim
      buf_idx = (2*i-1) + stride1
      this%psend_buf( buf_idx   ) = this%cfr_re(0)%f1( i, j-gc(1) )
      this%psend_buf( buf_idx+1 ) = this%cfi_re(0)%f1( i, j-gc(1) )
    enddo
  enddo

  ! copy m > 0 modes
  offset = 2 * this%dim * n1p
  do m = 1, this%max_mode
    stride1 = 4 * this%dim * n1p * (m-1)
    do j = 1, n1p
      stride2 = 4 * this%dim * (j-1)
      do i = 1, this%dim
        buf_idx = (4*i-3) + stride1 + stride2 + offset
        this%psend_buf( buf_idx   ) = this%cfr_re(m)%f1( i, j-gc(1) )
        this%psend_buf( buf_idx+1 ) = this%cfr_im(m)%f1( i, j-gc(1) )
        this%psend_buf( buf_idx+2 ) = this%cfi_re(m)%f1( i, j-gc(1) )
        this%psend_buf( buf_idx+3 ) = this%cfi_im(m)%f1( i, j-gc(1) )
      enddo
    enddo
  enddo

  call mpi_isend( this%psend_buf, count, p_dtype_real, idproc_des, pp_msg%tag, comm_world(), pp_msg%id, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_ISEND failed.' )
  endif

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_send_f1

subroutine pipe_recv_f1( this, pp_msg, dir, mode )

  implicit none
  class( field_complex ), intent(inout) :: this
  type( ppmsg ), intent(inout) :: pp_msg
  character(len=*), intent(in) :: dir, mode

  integer :: i, j, m, idproc_src, n1p, stride1, stride2, offset, buf_idx
  integer :: count, ierr
  integer, dimension(2) :: gc
  integer, dimension(MPI_STATUS_SIZE) :: stat
  character(len=32), save :: sname = "pipe_precv_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  n1p = size( this%cfr_re(0)%f1, 2 )
  gc  = this%gc_num(:,1)

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

  count = this%dim * n1p * ( 4 * this%max_mode + 2 )
  call this%check_recv_buf( count )

  call mpi_recv( this%precv_buf, count, p_dtype_real, idproc_src, pp_msg%tag, comm_world(), stat, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_RECV failed.' )
  endif

  select case ( trim(mode) )
  case ( 'replace' )

    ! copy m = 0 mode
    do j = 1, n1p
      stride1 = 2 * this%dim * (j-1)
      do i = 1, this%dim
        buf_idx = (2*i-1) + stride1
        this%cfr_re(0)%f1( i, j-gc(1) ) = this%precv_buf( buf_idx   )
        this%cfi_re(0)%f1( i, j-gc(1) ) = this%precv_buf( buf_idx+1 )
      enddo
    enddo

    ! copy m > 0 modes
    offset = 2 * this%dim * n1p
    do m = 1, this%max_mode
      stride1 = 4 * this%dim * n1p * (m-1)
      do j = 1, n1p
        stride2 = 4 * this%dim * (j-1)
        do i = 1, this%dim
          buf_idx = (4*i-3) + stride1 + stride2 + offset
          this%cfr_re(m)%f1( i, j-gc(1) ) = this%precv_buf( buf_idx   )
          this%cfr_im(m)%f1( i, j-gc(1) ) = this%precv_buf( buf_idx+1 )
          this%cfi_re(m)%f1( i, j-gc(1) ) = this%precv_buf( buf_idx+2 )
          this%cfi_im(m)%f1( i, j-gc(1) ) = this%precv_buf( buf_idx+3 )
        enddo
      enddo
    enddo

  case ( 'add' )

    ! copy m=0 mode
    do j = 1, n1p
      stride1 = 2 * this%dim * (j-1)
      do i = 1, this%dim
        buf_idx = (2*i-1) + stride1
        this%cfr_re(0)%f1( i, j-gc(1) ) = &
        this%cfr_re(0)%f1( i, j-gc(1) ) + this%precv_buf( buf_idx   )
        this%cfi_re(0)%f1( i, j-gc(1) ) = &
        this%cfi_re(0)%f1( i, j-gc(1) ) + this%precv_buf( buf_idx+1 )
      enddo
    enddo

    ! copy m > 0 mode
    offset = 2 * this%dim * n1p
    do m = 1, this%max_mode
      stride1 = 4 * this%dim * n1p * (m-1)
      do j = 1, n1p
        stride2 = 4 * this%dim * (j-1)
        do i = 1, this%dim
          buf_idx = (4*i-3) + stride1 + stride2 + offset
          this%cfr_re(m)%f1( i, j-gc(1) ) = &
          this%cfr_re(m)%f1( i, j-gc(1) ) + this%precv_buf( buf_idx   )
          this%cfr_im(m)%f1( i, j-gc(1) ) = &
          this%cfr_im(m)%f1( i, j-gc(1) ) + this%precv_buf( buf_idx+1 )
          this%cfi_re(m)%f1( i, j-gc(1) ) = &
          this%cfi_re(m)%f1( i, j-gc(1) ) + this%precv_buf( buf_idx+2 )
          this%cfi_im(m)%f1( i, j-gc(1) ) = &
          this%cfi_im(m)%f1( i, j-gc(1) ) + this%precv_buf( buf_idx+3 )
        enddo
      enddo
    enddo

  case default
    call write_err( 'Invalid data transfer mode! Valid options are "replace" &
      &and "add".' )
  end select

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_recv_f1

subroutine check_send_buf( this, size_min )
  implicit none
  class( field_complex ), intent(inout) :: this
  integer, intent(in) :: size_min
  character(len=32), save :: sname = 'check_send_buf'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  if ( .not. allocated( this%psend_buf ) ) then
    allocate( this%psend_buf( size_min ) )
  elseif ( size( this%psend_buf ) < size_min ) then
    deallocate( this%psend_buf )
    allocate( this%psend_buf( size_min ) )
  endif
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine check_send_buf

subroutine check_recv_buf( this, size_min )
  implicit none
  class( field_complex ), intent(inout) :: this
  integer, intent(in) :: size_min
  character(len=32), save :: sname = 'check_recv_buf'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  if ( .not. allocated( this%precv_buf ) ) then
    allocate( this%precv_buf( size_min ) )
  elseif ( size( this%precv_buf ) < size_min ) then
    deallocate( this%precv_buf )
    allocate( this%precv_buf( size_min ) )
  endif
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine check_recv_buf

subroutine write_hdf5_single( this, files, dim )

  implicit none

  class( field_complex ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim

  integer :: i
  character(len=32), save :: sname = 'write_hdf5_single'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode

    if ( i == 0 ) then
      call this%cfr_re(i)%write_hdf5( files(1), dim )
      call this%cfi_re(i)%write_hdf5( files(2), dim )
      cycle
    endif

    call this%cfr_re(i)%write_hdf5( files(4*i-1), dim )
    call this%cfr_im(i)%write_hdf5( files(4*i  ), dim )
    call this%cfi_re(i)%write_hdf5( files(4*i+1), dim )
    call this%cfi_im(i)%write_hdf5( files(4*i+2), dim )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_single

subroutine write_hdf5_pipe( this, files, dim, rtag, stag, id )

  implicit none

  class( field_complex ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim, rtag, stag
  integer, intent(inout) :: id

  integer :: i
  character(len=32), save :: sname = 'write_hdf5_pipe'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode

    if ( i == 0 ) then
      call this%cfr_re(i)%write_hdf5( files(1), dim, rtag, stag, id )
      call this%cfi_re(i)%write_hdf5( files(2), dim, rtag, stag, id )
      cycle
    endif

    call this%cfr_re(i)%write_hdf5( files(4*i-1), dim, rtag, stag, id )
    call this%cfr_im(i)%write_hdf5( files(4*i  ), dim, rtag, stag, id )
    call this%cfi_re(i)%write_hdf5( files(4*i+1), dim, rtag, stag, id )
    call this%cfi_im(i)%write_hdf5( files(4*i+2), dim, rtag, stag, id )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_pipe

subroutine read_hdf5( this, files, dim )

  implicit none

  class( field_complex ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim

  integer :: i
  character(len=32), save :: sname = 'read_hdf5'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode

    if ( i == 0 ) then
      call this%cfr_re(i)%read_hdf5( files(1), dim )
      call this%cfi_re(i)%read_hdf5( files(2), dim )
      cycle
    endif

    call this%cfr_re(i)%read_hdf5( files(4*i-1), dim )
    call this%cfr_im(i)%read_hdf5( files(4*i  ), dim )
    call this%cfi_re(i)%read_hdf5( files(4*i+1), dim )
    call this%cfi_im(i)%read_hdf5( files(4*i+2), dim )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine read_hdf5

function get_cfr_re_all( this )
  implicit none
  class( field_complex ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_cfr_re_all
  get_cfr_re_all => this%cfr_re
end function get_cfr_re_all

function get_cfr_re_mode( this, mode )
  implicit none
  class( field_complex ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_cfr_re_mode
  get_cfr_re_mode => this%cfr_re(mode)
end function get_cfr_re_mode

function get_cfr_im_all( this )
  implicit none
  class( field_complex ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_cfr_im_all
  get_cfr_im_all => this%cfr_im
end function get_cfr_im_all

function get_cfr_im_mode( this, mode )
  implicit none
  class( field_complex ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_cfr_im_mode
  get_cfr_im_mode => this%cfr_im(mode)
end function get_cfr_im_mode

function get_cfi_re_all( this )
  implicit none
  class( field_complex ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_cfi_re_all
  get_cfi_re_all => this%cfi_re
end function get_cfi_re_all

function get_cfi_re_mode( this, mode )
  implicit none
  class( field_complex ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_cfi_re_mode
  get_cfi_re_mode => this%cfi_re(mode)
end function get_cfi_re_mode

function get_cfi_im_all( this )
  implicit none
  class( field_complex ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_cfi_im_all
  get_cfi_im_all => this%cfi_im
end function get_cfi_im_all

function get_cfi_im_mode( this, mode )
  implicit none
  class( field_complex ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_cfi_im_mode
  get_cfi_im_mode => this%cfi_im(mode)
end function get_cfi_im_mode

subroutine assign_f1( this, that )

  implicit none

  class( field_complex ), intent(inout) :: this
  class(*), intent(in) :: that

  integer :: i

  select type (that)

  type is (real)

    do i = 0, this%max_mode
      this%cfr_re(i) = that
      this%cfi_re(i) = that
      if (i == 0) cycle
      this%cfr_im(i) = that
      this%cfi_im(i) = that
    enddo

  class is ( field_complex )

    do i = 0, this%max_mode
      this%cfr_re(i) = that%cfr_re(i)
      this%cfi_re(i) = that%cfi_re(i)
      if (i == 0) cycle
      this%cfr_im(i) = that%cfr_im(i)
      this%cfi_im(i) = that%cfi_im(i)
    enddo

  class default
    call write_err( "invalid assignment type!" )
  end select

end subroutine assign_f1

subroutine assign_f2( this, that )

  implicit none

  class( field_complex ), intent(inout) :: this
  class(*), intent(in) :: that

  integer :: i

  select type (that)

  type is (real)

    do i = 0, this%max_mode
      call this%cfr_re(i)%as( that )
      call this%cfi_re(i)%as( that )
      if (i == 0) cycle
      call this%cfr_im(i)%as( that )
      call this%cfi_im(i)%as( that )
    enddo

  class is ( field_complex )

    do i = 0, this%max_mode
      call this%cfr_re(i)%as( that%get_cfr_re(i) )
      call this%cfi_re(i)%as( that%get_cfi_re(i) )
      if (i == 0) cycle
      call this%cfr_im(i)%as( that%get_cfr_im(i) )
      call this%cfi_im(i)%as( that%get_cfi_im(i) )
    enddo

  class default
    call write_err( "invalid assignment type!" )
  end select

end subroutine assign_f2

subroutine add_f1_binary_dim( a1, a2, a3, dim1, dim2, dim3 )

  implicit none

  class( field_complex ), intent(in) :: a1, a2
  class( field_complex ), intent(inout) :: a3
  integer, intent(in), dimension(:) :: dim1, dim2, dim3

  class( ufield ), dimension(:), pointer :: ua3r_re => null(), ua3r_im => null()
  class( ufield ), dimension(:), pointer :: ua3i_re => null(), ua3i_im => null()
  integer :: i

  ua3r_re => a3%get_cfr_re()
  ua3r_im => a3%get_cfr_im()
  ua3i_re => a3%get_cfi_re()
  ua3i_im => a3%get_cfi_im()

  do i = 0, a1%max_mode
    call add_f1( a1%cfr_re(i), a2%cfr_re(i), ua3r_re(i), dim1, dim2, dim3 )
    call add_f1( a1%cfi_re(i), a2%cfi_re(i), ua3i_re(i), dim1, dim2, dim3 )
    if (i==0) cycle
    call add_f1( a1%cfr_im(i), a2%cfr_im(i), ua3r_im(i), dim1, dim2, dim3 )
    call add_f1( a1%cfi_im(i), a2%cfi_im(i), ua3i_im(i), dim1, dim2, dim3 )
  enddo

end subroutine add_f1_binary_dim

subroutine add_f1_unitary_dim( a1, a2, dim1, dim2 )

  implicit none

  class( field_complex ), intent(in) :: a1
  class( field_complex ), intent(inout) :: a2
  integer, intent(in), dimension(:) :: dim1, dim2

  class( ufield ), dimension(:), pointer :: ua2r_re => null(), ua2r_im => null()
  class( ufield ), dimension(:), pointer :: ua2i_re => null(), ua2i_im => null()
  integer :: i

  ua2r_re => a2%get_cfr_re()
  ua2r_im => a2%get_cfr_im()
  ua2i_re => a2%get_cfi_re()
  ua2i_im => a2%get_cfi_im()

  do i = 0, a1%max_mode
    call add_f1( a1%cfr_re(i), ua2r_re(i), dim1, dim2 )
    call add_f1( a1%cfi_re(i), ua2i_re(i), dim1, dim2 )
    if (i==0) cycle
    call add_f1( a1%cfr_im(i), ua2r_im(i), dim1, dim2 )
    call add_f1( a1%cfi_im(i), ua2i_im(i), dim1, dim2 )
  enddo

end subroutine add_f1_unitary_dim

subroutine add_f1_binary( a1, a2, a3 )

  implicit none

  class( field_complex ), intent(in) :: a1, a2
  class( field_complex ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3r_re => null(), ua3r_im => null()
  class( ufield ), dimension(:), pointer :: ua3i_re => null(), ua3i_im => null()
  integer :: i

  ua3r_re => a3%get_cfr_re()
  ua3r_im => a3%get_cfr_im()
  ua3i_re => a3%get_cfi_re()
  ua3i_im => a3%get_cfi_im()

  do i = 0, a1%max_mode
    call add_f1( a1%cfr_re(i), a2%cfr_re(i), ua3r_re(i) )
    call add_f1( a1%cfi_re(i), a2%cfi_re(i), ua3i_re(i) )
    if (i==0) cycle
    call add_f1( a1%cfr_im(i), a2%cfr_im(i), ua3r_im(i) )
    call add_f1( a1%cfi_im(i), a2%cfi_im(i), ua3i_im(i) )
  enddo

end subroutine add_f1_binary

subroutine add_f1_unitary( a1, a2 )

  implicit none

  class( field_complex ), intent(in) :: a1
  class( field_complex ), intent(inout) :: a2

  class( ufield ), dimension(:), pointer :: ua2r_re => null(), ua2r_im => null()
  class( ufield ), dimension(:), pointer :: ua2i_re => null(), ua2i_im => null()
  integer :: i

  ua2r_re => a2%get_cfr_re()
  ua2r_im => a2%get_cfr_im()
  ua2i_re => a2%get_cfi_re()
  ua2i_im => a2%get_cfi_im()

  do i = 0, a1%max_mode
    call add_f1( a1%cfr_re(i), ua2r_re(i) )
    call add_f1( a1%cfi_re(i), ua2i_re(i) )
    if (i==0) cycle
    call add_f1( a1%cfr_im(i), ua2r_im(i) )
    call add_f1( a1%cfi_im(i), ua2i_im(i) )
  enddo

end subroutine add_f1_unitary

subroutine sub_f1_binary_dim( a1, a2, a3, dim1, dim2, dim3 )

  implicit none

  class( field_complex ), intent(in) :: a1, a2
  class( field_complex ), intent(inout) :: a3
  integer, intent(in), dimension(:) :: dim1, dim2, dim3

  class( ufield ), dimension(:), pointer :: ua3r_re => null(), ua3r_im => null()
  class( ufield ), dimension(:), pointer :: ua3i_re => null(), ua3i_im => null()
  integer :: i

  ua3r_re => a3%get_cfr_re()
  ua3r_im => a3%get_cfr_im()
  ua3i_re => a3%get_cfi_re()
  ua3i_im => a3%get_cfi_im()

  do i = 0, a1%max_mode
    call sub_f1( a1%cfr_re(i), a2%cfr_re(i), ua3r_re(i), dim1, dim2, dim3 )
    call sub_f1( a1%cfi_re(i), a2%cfi_re(i), ua3i_re(i), dim1, dim2, dim3 )
    if (i==0) cycle
    call sub_f1( a1%cfr_im(i), a2%cfr_im(i), ua3r_im(i), dim1, dim2, dim3 )
    call sub_f1( a1%cfi_im(i), a2%cfi_im(i), ua3i_im(i), dim1, dim2, dim3 )
  enddo

end subroutine sub_f1_binary_dim

subroutine sub_f1_unitary_dim( a1, a2, dim1, dim2 )

  implicit none

  class( field_complex ), intent(in) :: a1
  class( field_complex ), intent(inout) :: a2
  integer, intent(in), dimension(:) :: dim1, dim2

  class( ufield ), dimension(:), pointer :: ua2r_re => null(), ua2r_im => null()
  class( ufield ), dimension(:), pointer :: ua2i_re => null(), ua2i_im => null()
  integer :: i

  ua2r_re => a2%get_cfr_re()
  ua2r_im => a2%get_cfr_im()
  ua2i_re => a2%get_cfi_re()
  ua2i_im => a2%get_cfi_im()

  do i = 0, a1%max_mode
    call sub_f1( a1%cfr_re(i), ua2r_re(i), dim1, dim2 )
    call sub_f1( a1%cfi_re(i), ua2i_re(i), dim1, dim2 )
    if (i==0) cycle
    call sub_f1( a1%cfr_im(i), ua2r_im(i), dim1, dim2 )
    call sub_f1( a1%cfi_im(i), ua2i_im(i), dim1, dim2 )
  enddo

end subroutine sub_f1_unitary_dim

subroutine sub_f1_binary( a1, a2, a3 )

  implicit none

  class( field_complex ), intent(in) :: a1, a2
  class( field_complex ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3r_re => null(), ua3r_im => null()
  class( ufield ), dimension(:), pointer :: ua3i_re => null(), ua3i_im => null()
  integer :: i

  ua3r_re => a3%get_cfr_re()
  ua3r_im => a3%get_cfr_im()
  ua3i_re => a3%get_cfi_re()
  ua3i_im => a3%get_cfi_im()

  do i = 0, a1%max_mode
    call sub_f1( a1%cfr_re(i), a2%cfr_re(i), ua3r_re(i) )
    call sub_f1( a1%cfi_re(i), a2%cfi_re(i), ua3i_re(i) )
    if (i==0) cycle
    call sub_f1( a1%cfr_im(i), a2%cfr_im(i), ua3r_im(i) )
    call sub_f1( a1%cfi_im(i), a2%cfi_im(i), ua3i_im(i) )
  enddo

end subroutine sub_f1_binary

subroutine sub_f1_unitary( a1, a2 )

  implicit none

  class( field_complex ), intent(in) :: a1
  class( field_complex ), intent(inout) :: a2

  class( ufield ), dimension(:), pointer :: ua2r_re => null(), ua2r_im => null()
  class( ufield ), dimension(:), pointer :: ua2i_re => null(), ua2i_im => null()
  integer :: i

  ua2r_re => a2%get_cfr_re()
  ua2r_im => a2%get_cfr_im()
  ua2i_re => a2%get_cfi_re()
  ua2i_im => a2%get_cfi_im()

  do i = 0, a1%max_mode
    call sub_f1( a1%cfr_re(i), ua2r_re(i) )
    call sub_f1( a1%cfi_re(i), ua2i_re(i) )
    if (i==0) cycle
    call sub_f1( a1%cfr_im(i), ua2r_im(i) )
    call sub_f1( a1%cfi_im(i), ua2i_im(i) )
  enddo

end subroutine sub_f1_unitary

subroutine dot_f1_unitary( a1, a2 )

  implicit none

  real, intent(in) :: a1
  class( field_complex ), intent(inout) :: a2

  class( ufield ), dimension(:), pointer :: ua2r_re => null(), ua2r_im => null()
  class( ufield ), dimension(:), pointer :: ua2i_re => null(), ua2i_im => null()
  integer :: i

  ua2r_re => a2%get_cfr_re()
  ua2i_re => a2%get_cfi_re()
  ua2r_im => a2%get_cfr_im()
  ua2i_im => a2%get_cfi_im()

  do i = 0, a2%max_mode
    call dot_f1( a1, ua2r_re(i) )
    call dot_f1( a1, ua2i_re(i) )
    if (i==0) cycle
    call dot_f1( a1, ua2r_im(i) )
    call dot_f1( a1, ua2i_im(i) )
  enddo

end subroutine dot_f1_unitary

subroutine dot_f1_unitary_dim( a1, a2, dim )

  implicit none

  real, intent(in) :: a1
  class( field_complex ), intent(inout) :: a2
  integer, intent(in), dimension(:) :: dim

  class( ufield ), dimension(:), pointer :: ua2r_re => null(), ua2r_im => null()
  class( ufield ), dimension(:), pointer :: ua2i_re => null(), ua2i_im => null()
  integer :: i

  ua2r_re => a2%get_cfr_re()
  ua2i_re => a2%get_cfi_re()
  ua2r_im => a2%get_cfr_im()
  ua2i_im => a2%get_cfi_im()

  do i = 0, a2%max_mode
    call dot_f1( a1, ua2r_re(i), dim )
    call dot_f1( a1, ua2i_re(i), dim )
    if (i==0) cycle
    call dot_f1( a1, ua2r_im(i), dim )
    call dot_f1( a1, ua2i_im(i), dim )
  enddo

end subroutine dot_f1_unitary_dim

subroutine add_f2_binary( a1, a2, a3 )

  implicit none

  class( field_complex ), intent(in) :: a1, a2
  class( field_complex ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3r_re => null(), ua3r_im => null()
  class( ufield ), dimension(:), pointer :: ua3i_re => null(), ua3i_im => null()
  integer :: i

  ua3r_re => a3%get_cfr_re()
  ua3i_re => a3%get_cfi_re()
  ua3r_im => a3%get_cfr_im()
  ua3i_im => a3%get_cfi_im()

  do i = 0, a1%max_mode
    call add_f2( a1%cfr_re(i), a2%cfr_re(i), ua3r_re(i) )
    call add_f2( a1%cfi_re(i), a2%cfi_re(i), ua3i_re(i) )
    if (i==0) cycle
    call add_f2( a1%cfr_im(i), a2%cfr_im(i), ua3r_im(i) )
    call add_f2( a1%cfi_im(i), a2%cfi_im(i), ua3i_im(i) )
  enddo

end subroutine add_f2_binary

subroutine add_f2_unitary( a1, a2 )

  implicit none

  class( field_complex ), intent(in) :: a1
  class( field_complex ), intent(inout) ::a2

  class( ufield ), dimension(:), pointer :: ua2r_re => null(), ua2r_im => null()
  class( ufield ), dimension(:), pointer :: ua2i_re => null(), ua2i_im => null()
  integer :: i

  ua2r_re => a2%get_cfr_re()
  ua2i_re => a2%get_cfi_re()
  ua2r_im => a2%get_cfr_im()
  ua2i_im => a2%get_cfi_im()

  do i = 0, a1%max_mode
    call add_f2( a1%cfr_re(i), ua2r_re(i) )
    call add_f2( a1%cfi_re(i), ua2i_re(i) )
    if (i==0) cycle
    call add_f2( a1%cfr_im(i), ua2r_im(i) )
    call add_f2( a1%cfi_im(i), ua2i_im(i) )
  enddo

end subroutine add_f2_unitary

subroutine sub_f2_binary( a1, a2, a3 )

  implicit none

  class( field_complex ), intent(in) :: a1, a2
  class( field_complex ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3r_re => null(), ua3r_im => null()
  class( ufield ), dimension(:), pointer :: ua3i_re => null(), ua3i_im => null()
  integer :: i

  ua3r_re => a3%get_cfr_re()
  ua3i_re => a3%get_cfi_re()
  ua3r_im => a3%get_cfr_im()
  ua3i_im => a3%get_cfi_im()

  do i = 0, a1%max_mode
    call sub_f2( a1%cfr_re(i), a2%cfr_re(i), ua3r_re(i) )
    call sub_f2( a1%cfi_re(i), a2%cfi_re(i), ua3i_re(i) )
    if (i==0) cycle
    call sub_f2( a1%cfr_im(i), a2%cfr_im(i), ua3r_im(i) )
    call sub_f2( a1%cfi_im(i), a2%cfi_im(i), ua3i_im(i) )
  enddo

end subroutine sub_f2_binary

subroutine sub_f2_unitary( a1, a2 )

  implicit none

  class( field_complex ), intent(in) :: a1
  class( field_complex ), intent(inout) ::a2

  class( ufield ), dimension(:), pointer :: ua2r_re => null(), ua2r_im => null()
  class( ufield ), dimension(:), pointer :: ua2i_re => null(), ua2i_im => null()
  integer :: i

  ua2r_re => a2%get_cfr_re()
  ua2i_re => a2%get_cfi_re()
  ua2r_im => a2%get_cfr_im()
  ua2i_im => a2%get_cfi_im()

  do i = 0, a1%max_mode
    call sub_f2( a1%cfr_re(i), ua2r_re(i) )
    call sub_f2( a1%cfi_re(i), ua2i_re(i) )
    if (i==0) cycle
    call sub_f2( a1%cfr_im(i), ua2r_im(i) )
    call sub_f2( a1%cfi_im(i), ua2i_im(i) )
  enddo

end subroutine sub_f2_unitary

end module field_complex_class