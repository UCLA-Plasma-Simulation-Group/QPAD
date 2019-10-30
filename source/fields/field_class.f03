module field_class

use parallel_pipe_class
use grid_class
use ufield_class
use ufield_smooth_class
use hdf5io_class
use param
use sys
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

  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()
  class( ufield ), dimension(:), pointer :: rf_re => null()
  class( ufield ), dimension(:), pointer :: rf_im => null()
  type( ufield_smooth ) :: smooth

  real :: dr, dxi
  integer :: num_modes, dim
  integer :: entity
  real, dimension(:,:,:,:), allocatable :: psend_buf, precv_buf

  contains

  generic :: new => init_field, init_field_cp
  procedure :: del => end_field
  generic :: get_rf_re => get_rf_re_all, get_rf_re_mode
  generic :: get_rf_im => get_rf_im_all, get_rf_im_mode
  generic :: write_hdf5 => write_hdf5_single, write_hdf5_pipe
  procedure :: smooth_f1
  procedure :: copy_slice
  procedure :: get_dr, get_dxi, get_num_modes, get_dim
  procedure :: copy_gc_f1, copy_gc_f2
  procedure :: acopy_gc_f1, acopy_gc_f2
  procedure :: pipe_send, pipe_recv
  procedure :: pipe_gc_send, pipe_gc_recv

  procedure, private :: init_field, init_field_cp, end_field
  procedure, private :: get_rf_re_all, get_rf_re_mode, get_rf_im_all, get_rf_im_mode
  procedure, private :: write_hdf5_single, write_hdf5_pipe

  generic :: assignment(=)   => assign_f1
  generic :: as              => assign_f2

  procedure, private :: assign_f1, assign_f2

end type field


contains

! =====================================================================
! Class field implementation
! =====================================================================
subroutine init_field( this, pp, gp, dim, num_modes, gc_num, &
  entity, smooth_type, smooth_order )

  implicit none

  class( field ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, dim
  integer, intent(in), dimension(2,2) :: gc_num
  integer, intent(in), optional :: entity, smooth_type, smooth_order

  integer :: i
  integer, dimension(2,2) :: gc_num_new
  character(len=20), save :: sname = "init_field"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%pp        => pp
  this%gp        => gp
  this%dim       = dim
  this%num_modes = num_modes
  this%dr        = gp%get_dr()
  this%dxi       = gp%get_dxi()

  if ( present(entity) ) then
    this%entity = entity
  else
    this%entity = p_entity_none
  endif

  if ( present(smooth_type) .and. present(smooth_order) ) then
    call this%smooth%new( smooth_type, smooth_order )
    gc_num_new(1,1) = max( gc_num(1,1), smooth_order )
    gc_num_new(2,1) = max( gc_num(2,1), smooth_order )
    gc_num_new(:,2) = gc_num(:,2)
  else
    call this%smooth%new( p_smooth_none, 0 )
  endif

  allocate( this%rf_re(0:num_modes) )
  allocate( this%rf_im(num_modes) )
  do i = 0, this%num_modes
    call this%rf_re(i)%new( pp, gp, dim, i, gc_num, has_2d=.true. )
    if (i==0) cycle
    call this%rf_im(i)%new( pp, gp, dim, i, gc_num, has_2d=.true. )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field

subroutine init_field_cp( this, that )

  implicit none

  class( field ), intent(inout) :: this
  class( field ), intent(in) :: that

  integer :: i

  this%pp        => that%pp
  this%gp        => that%gp
  this%dim       = that%get_dim()
  this%num_modes = that%get_num_modes()
  this%dr        = that%get_dr()
  this%dxi       = that%get_dxi()

  allocate( this%rf_re(0:this%num_modes) )
  allocate( this%rf_im(this%num_modes) )
  do i = 0, this%num_modes
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
  do i = 1, this%num_modes
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

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_slice( idx, dir )
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_slice( idx, dir )
  enddo

  call stop_tprof( 'copy slices' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_slice

subroutine copy_gc_f1( this, bnd_ax )

  implicit none

  class( field ), intent(inout) :: this
  logical, intent(in) :: bnd_ax

  integer :: i
  character(len=20), save :: sname = "copy_gc_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_gc_f1( bnd_ax )
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_gc_f1( bnd_ax )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_f1

subroutine copy_gc_f2( this, bnd_ax )

  implicit none

  class( field ), intent(inout) :: this
  logical, intent(in) :: bnd_ax

  integer :: i
  character(len=20), save :: sname = "copy_gc_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_gc_f2( bnd_ax )
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_gc_f2( bnd_ax )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_f2

subroutine acopy_gc_f1( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "acopy_gc_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%acopy_gc_f1()
    if ( i == 0 ) cycle
    call this%rf_im(i)%acopy_gc_f1()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine acopy_gc_f1

subroutine acopy_gc_f2( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "acopy_gc_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%acopy_gc_f2()
    if ( i == 0 ) cycle
    call this%rf_im(i)%acopy_gc_f2()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine acopy_gc_f2

subroutine pipe_gc_send( this, tag, sid )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: tag
  integer, intent(inout) :: sid

  integer :: i, j, k, m
  integer :: idproc, idproc_des, stageid, lnvp, comm
  integer :: n1p, count, dtype, ierr
  integer, dimension(2) :: gc1, gc2
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:,:,:,:), allocatable, save :: buf
  character(len=20), save :: sname = "pipe_gc_send"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  lnvp       = this%pp%getlnvp()
  stageid    = this%pp%getstageid()
  idproc     = this%pp%getidproc()
  idproc_des = idproc - lnvp
  n1p        = size(this%rf_re(0)%f1, 2)
  comm       = this%pp%getlworld()
  dtype      = this%pp%getmreal()

  gc1 = this%rf_re(0)%get_gc_num(1)
  gc2 = this%rf_re(0)%get_gc_num(2)

  if ( gc2(p_upper) <= 0 ) call write_err( 'No guard cells for pipeline guard cell copy!' )

  if ( stageid == 0 ) then
    sid = MPI_REQUEST_NULL
    call stop_tprof( 'pipeline' )
    call write_dbg( cls_name, sname, cls_level, 'ends' )
    return
  endif

  if ( .not. allocated(buf) ) then
    allocate( buf(this%dim, n1p, gc2(p_upper), 0:2*this%num_modes) )
  endif

  ! pack m=0 mode
  do k = 1, gc2(p_upper)
    do j = 1, n1p
      do i = 1, this%dim
        buf(i,j,k,0) = this%rf_re(0)%f2(i, j-gc1(p_lower), k)
      enddo
    enddo
  enddo

  ! pack m>0 modes
  if (this%num_modes > 0) then
    do m = 1, this%num_modes
      do k = 1, gc2(p_upper)
        do j = 1, n1p
          do i = 1, this%dim
            buf(i,j,k,2*m-1) = this%rf_re(m)%f2(i, j-gc1(p_lower), k)
            buf(i,j,k,2*m  ) = this%rf_im(m)%f2(i, j-gc1(p_lower), k)
          enddo
        enddo
      enddo
    enddo
  endif

  count = size(buf)
  call MPI_ISEND( buf, count, dtype, idproc_des, tag, comm, sid, ierr )
  ! check for error
  if (ierr /= 0) call write_err('MPI_ISEND failed.')

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_gc_send

subroutine pipe_gc_recv( this, tag )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: tag
  ! integer, intent(inout) :: rid

  integer :: i, j, k, m
  integer :: idproc, idproc_src, stageid, nstage, lnvp, comm
  integer :: n1p, nzp, count, dtype, ierr
  integer, dimension(2) :: gc1, gc2
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:,:,:,:), allocatable, save :: buf
  character(len=20), save :: sname = "pipe_gc_recv"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  lnvp       = this%pp%getlnvp()
  stageid    = this%pp%getstageid()
  nstage     = this%pp%getnstage()
  idproc     = this%pp%getidproc()
  idproc_src = idproc + lnvp
  n1p        = size(this%rf_re(0)%f1, 2)
  nzp        = this%gp%get_ndp(2)
  comm       = this%pp%getlworld()
  dtype      = this%pp%getmreal()

  gc1 = this%rf_re(0)%get_gc_num(1)
  gc2 = this%rf_re(0)%get_gc_num(2)

  if ( gc2(p_upper) <= 0 ) call write_err( 'No guard cells for pipeline guard cell copy!' )

  if ( stageid == nstage-1 ) then
    call stop_tprof( 'pipeline' )
    call write_dbg( cls_name, sname, cls_level, 'ends' )
    return
  endif

  if ( .not. allocated(buf) ) then
    allocate( buf(this%dim, n1p, gc2(p_upper), 0:2*this%num_modes) )
  endif

  count = size(buf)
  call MPI_RECV( buf, count, dtype, idproc_src, tag, comm, stat, ierr )
  ! check for error
  if (ierr /= 0) call write_err('MPI_RECV failed.')

  ! unpack m=0 mode
  do k = 1, gc2(p_upper)
    do j = 1, n1p
      do i = 1, this%dim
        this%rf_re(0)%f2(i, j-gc1(p_lower), nzp+k) = buf(i,j,k,0)
      enddo
    enddo
  enddo

  ! unpack m>0 modes
  if (this%num_modes > 0) then
    do m = 1, this%num_modes
      do k = 1, gc2(p_upper)
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_re(m)%f2(i, j-gc1(p_lower), nzp+k) = buf(i,j,k,2*m-1)
            this%rf_im(m)%f2(i, j-gc1(p_lower), nzp+k) = buf(i,j,k,2*m  )
          enddo
        enddo
      enddo
    enddo
  endif

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_gc_recv

! subroutine copy_gc_pipe( this, rtag, stag, rid, sid )

!   implicit none

!   class( field ), intent(inout) :: this
!   integer, intent(in) :: rtag, stag
!   integer, intent(inout) :: rid, sid

!   integer :: i, dir
!   integer :: idproc, idproc_next, idproc_last, nstage, stageid, nvp, comm
!   integer :: nrp, nzp, count, dtype, gc_num, m
!   integer :: ierr
!   integer, dimension(MPI_STATUS_SIZE) :: stat
!   real, dimension(:,:,:), allocatable, save :: rbuf, sbuf
!   character(len=20), save :: sname = "copy_gc_pipe"

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   nvp = this%gp%nvp(1)
!   nstage = this%pp%getnstage()
!   stageid = this%pp%getstageid()
!   idproc = this%pp%getidproc()
!   idproc_last = idproc - nvp
!   idproc_next = idproc + nvp
!   nrp = this%gp%ndp(1)
!   nzp = this%gp%ndp(2)
!   comm = this%pp%getlworld()
!   dtype = this%pp%getmreal()

!   gc_num = this%rf_re(0)%gc_num(p_upper,2)

!   if ( gc_num > 0 ) then

!     if (this%num_modes == 0) then
!       m = 1
!     else
!       m = 2 * this%num_modes + 1
!     endif
!     count = this%dim * nrp * gc_num * m

!     if ( .not. allocated(rbuf) ) then
!       allocate( rbuf(this%dim, nrp, gc_num, m), sbuf(this%dim, nrp, gc_num, m) )
!     endif

!     ! receiver
!     if ( stageid < nstage-1 ) then
!       call MPI_IRECV( this%f2(1,1-this%gc_num(p_lower,1),nzp+1), &
!         count, dtype, idproc_next, rtag, comm, rid, ierr )
!     else
!       rid = MPI_REQUEST_NULL
!     endif
!     ! sender
!     if ( stageid > 0 ) then
!       call MPI_ISEND( this%f2(1,1-this%gc_num(p_lower,1),1), &
!         count, dtype, idproc_last, stag, comm, sid, ierr )
!     else
!       sid = MPI_REQUEST_NULL
!     endif

!   else
!     call write_err( 'No guard cells for backward copy!' )
!   endif

!   do i = 0, this%num_modes
!     call this%rf_re(i)%copy_gc_pipe( dir, rtag, stag, rid, sid )
!     if ( i == 0 ) cycle
!     call this%rf_im(i)%copy_gc_pipe( dir, rtag, stag, rid, sid )
!   enddo

!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine copy_gc_pipe

subroutine pipe_send( this, stag, id, nslice )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: stag
  integer, intent(inout) :: id
  integer, intent(in), optional :: nslice

  integer :: i, j, k, m, ns
  integer :: idproc, idproc_des, lnvp, comm, stageid, nstage
  integer :: nzp, n1p, count, dtype, ierr
  integer, dimension(2) :: gc
  character(len=20), save :: sname = "pipe_send"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  ns = 1
  if ( present(nslice) ) ns = nslice

  lnvp       = this%pp%getlnvp()
  nstage     = this%pp%getnstage()
  stageid    = this%pp%getstageid()
  idproc     = this%pp%getidproc()
  idproc_des = idproc + lnvp
  nzp        = this%gp%get_ndp(2)
  n1p        = size(this%rf_re(0)%f1,2)
  comm       = this%pp%getlworld()
  dtype      = this%pp%getmreal()
  gc         = this%rf_re(0)%get_gc_num(1)

  if ( stageid == nstage-1 ) then
    id = MPI_REQUEST_NULL
    call stop_tprof( 'pipeline' )
    call write_dbg( cls_name, sname, cls_level, 'ends' )
    return
  endif

  if ( .not. allocated( this%psend_buf ) ) then
    allocate( this%psend_buf( this%dim, n1p, ns, 0:2*this%num_modes ) )
  endif

  ! copy m=0 mode
  do k = 1, ns
    do j = 1, n1p
      do i = 1, this%dim
        this%psend_buf(i,j,k,0) = this%rf_re(0)%f2( i, j-gc(1), nzp+1-ns+k )
      enddo
    enddo
  enddo

  ! copy m>0 mode
  if ( this%num_modes > 0 ) then
    do m = 1, this%num_modes
      do k = 1, ns
        do j = 1, n1p
          do i = 1, this%dim
            this%psend_buf(i,j,k,2*m-1) = this%rf_re(m)%f2( i, j-gc(1), nzp+1-ns+k )
            this%psend_buf(i,j,k,2*m  ) = this%rf_im(m)%f2( i, j-gc(1), nzp+1-ns+k )
          enddo
        enddo
      enddo
    enddo
  endif

  count = size( this%psend_buf )
  call MPI_ISEND( this%psend_buf, count, dtype, idproc_des, stag, comm, id, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_ISEND failed.' )
  endif

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_send

subroutine pipe_recv( this, rtag, nslice )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: rtag
  integer, intent(in), optional :: nslice

  integer :: i, j, k, m, ns
  integer :: idproc, idproc_src, lnvp, n1p, comm, stageid
  integer :: count, dtype, ierr, id
  integer, dimension(2) :: gc
  integer, dimension(MPI_STATUS_SIZE) :: stat
  character(len=20), save :: sname = "pipe_precv"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  ns = 1
  if ( present(nslice) ) ns = nslice

  lnvp       = this%pp%getlnvp()
  stageid    = this%pp%getstageid()
  idproc     = this%pp%getidproc()
  idproc_src = idproc - lnvp
  comm       = this%pp%getlworld()
  dtype      = this%pp%getmreal()
  n1p        = size(this%rf_re(0)%f1,2)
  gc         = this%rf_re(0)%get_gc_num(1)

  if ( stageid == 0 ) then
    call stop_tprof( 'pipeline' )
    call write_dbg( cls_name, sname, cls_level, 'ends' )
    return
  endif

  if ( .not. allocated( this%precv_buf ) ) then
    allocate( this%precv_buf( this%dim, n1p, ns, 0:2*this%num_modes ) )
  endif

  count = size( this%precv_buf )
  call MPI_RECV( this%precv_buf, count, dtype, idproc_src, rtag, comm, stat, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_RECV failed.' )
  endif

  ! copy m=0 mode
  do k = 1, ns
    do j = 1, n1p
      do i = 1, this%dim
        this%rf_re(0)%f2( i, j-gc(1), 1-ns+k ) = this%precv_buf(i,j,k,0)
      enddo
    enddo
  enddo

  ! copy m>0 mode
  if ( this%num_modes > 0 ) then
    do m = 1, this%num_modes
      do k = 1, ns
        do j = 1, n1p
          do i = 1, this%dim
            this%rf_re(m)%f2( i, j-gc(1), 1-ns+k ) = this%precv_buf(i,j,k,2*m-1)
            this%rf_im(m)%f2( i, j-gc(1), 1-ns+k ) = this%precv_buf(i,j,k,2*m)
          enddo
        enddo
      enddo
    enddo
  endif

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine pipe_recv

subroutine write_hdf5_single( this, files, dim )

  implicit none

  class( field ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim

  integer :: i
  character(len=32), save :: sname = 'write_hdf5_single'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes

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

  do i = 0, this%num_modes

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

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%smooth%smooth_f1( this%rf_re(i), i )
      cycle
    endif

    call this%smooth%smooth_f1( this%rf_re(i), i )
    call this%smooth%smooth_f1( this%rf_im(i), i )

  enddo

  call this%copy_gc_f1( bnd_ax = .false. )

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

function get_num_modes( this )

  implicit none

  class( field ), intent(in) :: this
  integer :: get_num_modes

  get_num_modes = this%num_modes

end function get_num_modes

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

    do i = 0, this%num_modes
      this%rf_re(i) = that
      if (i == 0) cycle
      this%rf_im(i) = that
    enddo

  class is (field)

    do i = 0, this%num_modes
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

    do i = 0, this%num_modes
      call this%rf_re(i)%as( that )
      if (i == 0) cycle
      call this%rf_im(i)%as( that )
    enddo

  class is (field)

    do i = 0, this%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a2%num_modes
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

  do i = 0, a2%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
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

  do i = 0, a1%num_modes
    call sub_f2( a1%rf_re(i), ua2_re(i) )
    if (i==0) cycle
    call sub_f2( a1%rf_im(i), ua2_im(i) )
  enddo

end subroutine sub_f2_unitary

end module field_class