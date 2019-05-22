module ufield_class

use parallel_pipe_class
use grid_class
use hdf5io_class
use param
use sys
use mpi

implicit none

private

character(len=20), parameter :: cls_name = "ufield"
integer, parameter :: cls_level = 4

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
end interface


public :: ufield
public :: add_f1, add_f2, sub_f1, sub_f2, dot_f1

! real, dimension(:), pointer, save :: buf => null() ! data buffer used for MPI

type :: ufield

  private

  class( parallel_pipe ), pointer, public :: pp => null()
  real, dimension(:,:), pointer, public :: f1 => null()
  real, dimension(:,:,:), pointer, public :: f2 => null()
  integer, dimension(2) :: nd ! number of global grid points
  integer, dimension(2) :: ndp ! number of local grid points
  integer, dimension(2) :: nvp
  integer :: dim ! dimension of array (guard cell excluded)
  integer, dimension(2) :: noff
  integer, dimension(2,2) :: gc_num ! number of guard cells
  integer :: mode
  logical :: has_2d

  contains

  generic :: new        => init_ufield, init_ufield_cp
  generic :: del        => end_ufield
  generic :: get_nd     => get_nd_all, get_nd_dim
  generic :: get_ndp    => get_ndp_all, get_ndp_dim
  generic :: get_gc_num => get_gc_num_all, get_gc_num_dim
  generic :: get_nvp    => get_nvp_all, get_nvp_dim
  generic :: get_noff   => get_noff_all, get_noff_dim
  generic :: write_hdf5 => write_hdf5_single, write_hdf5_pipe
  procedure :: get_dim
  procedure :: copy_slice
  procedure :: get_f1
  procedure :: get_f2
  procedure :: copy_gc_f1, copy_gc_f2, copy_gc_stage
  procedure :: acopy_gc_f1, acopy_gc_f2, acopy_gc_stage
  procedure :: has2d

  generic :: assignment(=)   => assign_f1
  generic :: as              => assign_f2
  ! generic :: operator(+)     => add_f1_v1, add_f1_v2
  ! generic :: operator(-)     => sub_f1_v1, sub_f1_v2
  ! generic :: operator(*)     => dot_f1_v1, dot_f1_v2
  ! generic :: operator(.add.) => add_f2_v1, add_f2_v2
  ! generic :: operator(.sub.) => sub_f2_v1, sub_f2_v2
  ! generic :: operator(.dot.) => dot_f2_v1, dot_f2_v2

  procedure, private :: init_ufield, init_ufield_cp, end_ufield
  procedure, private :: get_nd_all, get_nd_dim
  procedure, private :: get_ndp_all, get_ndp_dim
  procedure, private :: get_gc_num_all, get_gc_num_dim
  procedure, private :: get_nvp_all, get_nvp_dim
  procedure, private :: get_noff_all, get_noff_dim
  procedure, private :: write_hdf5_single, write_hdf5_pipe

  ! procedure, private, pass(a1) :: add_f1_v1, add_f1_v2
  ! procedure, private, pass(a1) :: dot_f1_v1, dot_f1_v2
  ! procedure, private, pass(a1) :: sub_f1_v1, sub_f1_v2
  ! procedure, private, pass(a1) :: add_f2_v1, add_f2_v2
  ! procedure, private, pass(a1) :: dot_f2_v1, dot_f2_v2
  ! procedure, private, pass(a1) :: sub_f2_v1, sub_f2_v2
  procedure, private :: assign_f1, assign_f2

end type ufield

contains

subroutine init_ufield( this, pp, gp, dim, mode, gc_num, has_2d )

  implicit none

  class( ufield ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: dim, mode
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
  this%mode = mode

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

subroutine init_ufield_cp( this, that )

  implicit none

  class( ufield ), intent(inout) :: this
  class( ufield ), intent(in) :: that

  this%pp => that%pp
  this%nd = that%get_nd()
  this%nvp = that%get_nvp()
  this%dim = that%dim
  this%ndp = that%get_ndp()
  this%noff = that%get_noff()
  this%gc_num = that%get_gc_num()
  this%mode = that%mode

  allocate( this%f1( this%dim, 1-this%gc_num(p_lower,1):this%ndp(1)+this%gc_num(p_upper,1) ) )
  this%f1 = 0.0

  if ( that%has2d() ) then
    this%has_2d = .true.
    allocate( this%f2( this%dim, &
              1-this%gc_num(p_lower,1):this%ndp(1)+this%gc_num(p_upper,1), &
              1-this%gc_num(p_lower,2):this%ndp(2)+this%gc_num(p_upper,2) ) )
    this%f2 = 0.0
  else
    this%has_2d = .false.
  endif

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

subroutine copy_gc_f1( this, bnd_ax )

  implicit none

  class( ufield ), intent(inout) :: this
  logical, intent(in) :: bnd_ax

  integer :: idproc, idproc_left, idproc_right, nvp, comm
  integer :: nrp, count, dtype
  integer :: tag = 1, msgid, ierr, i, j
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real :: pha


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
    else

      if ( bnd_ax ) then
        if ( mod(this%mode,2) == 0 ) then
          pha = 1.0
        else
          pha = -1.0
        endif
        do j = 1, this%gc_num(p_lower,1)
          do i = 1, this%dim
            this%f1(i,1-j) = pha * this%f1(i,j)
          enddo
        enddo
      endif

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

subroutine copy_gc_f2( this, bnd_ax )

  implicit none

  class( ufield ), intent(inout) :: this
  logical, intent(in) :: bnd_ax

  integer :: idproc, idproc_left, idproc_right, nvp, comm
  integer :: nrp, nzp, count, dtype, lgc1, ugc1, lgc2, ugc2
  integer :: tag = 1, msgid, ierr, i, j, k
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:,:,:), allocatable :: buf1, buf2
  real :: pha


  idproc = this%pp%getlidproc()
  idproc_left =  idproc - 1
  idproc_right = idproc + 1
  nrp = this%ndp(1)
  nzp = this%ndp(2)
  nvp = this%pp%getlnvp()
  comm = this%pp%getlgrp()
  dtype = this%pp%getmreal()

  allocate( buf1( this%dim, maxval(this%gc_num(:,1)), size(this%f2,3) ) )
  allocate( buf2( this%dim, maxval(this%gc_num(:,1)), size(this%f2,3) ) )

  lgc1 = this%gc_num(p_lower,1)
  ugc1 = this%gc_num(p_upper,1)
  lgc2 = this%gc_num(p_lower,2)
  ugc2 = this%gc_num(p_upper,2)

  ! forward message passing
  if ( lgc1 > 0 ) then
    count = this%dim * lgc1 * size(this%f2,3)
    ! receiver
    if ( idproc > 0 ) then
      buf1 = 0.0
      call MPI_IRECV( buf1, count, dtype, idproc_left, tag, comm, msgid, ierr )
    endif
    ! sender
    if ( idproc < nvp-1 ) then
      do k = 1, size(this%f2,3)
        do j = 1, lgc1
          do i = 1, this%dim
            buf2(i,j,k) = this%f2( i, nrp+j-lgc1, k-lgc2 )
          enddo
        enddo
      enddo
      call MPI_SEND( buf2, count, dtype, idproc_right, tag, comm, ierr )
    endif
    ! wait receiving finish
    if ( idproc > 0 ) then
      call MPI_WAIT( msgid, stat, ierr )
      do k = 1, size(this%f2,3)
        do j = 1, lgc1
          do i = 1, this%dim
            this%f2( i, j-lgc1, k-lgc2 ) = buf1(i,j,k)
          enddo
        enddo
      enddo

    else

      if ( bnd_ax ) then
        if ( mod(this%mode,2) == 0 ) then
          pha = 1.0
        else
          pha = -1.0
        endif
        do k = 1, size(this%f2,3)
          do j = 1, lgc1
            do i = 1, this%dim
              this%f2(i,1-j,k) = pha * this%f2(i,j,k)
            enddo
          enddo
        enddo
      endif

    endif
  endif

  ! backward message passing
  if ( ugc1 > 0 ) then
    count = this%dim * ugc1 * size(this%f2,3)
    ! receiver
    if ( idproc < nvp-1 ) then
      buf1 = 0.0
      call MPI_IRECV( buf1, count, dtype, idproc_right, tag, comm, msgid, ierr )
    endif
    ! sender
    if ( idproc > 0 ) then
      do k = 1, size(this%f2,3)
        do j = 1, ugc1
          do i = 1, this%dim
            buf2(i,j,k) = this%f2( i, j, k-lgc2 )
          enddo
        enddo
      enddo
      call MPI_SEND( buf2, count, dtype, idproc_left, tag, comm, ierr )
    endif
    ! wait receiving finish
    if ( idproc < nvp-1 ) then
      call MPI_WAIT( msgid, stat, ierr )
      do k = 1, size(this%f2,3)
        do j = 1, ugc1
          do i = 1, this%dim
            this%f2( i, nrp+j, k-lgc2 ) = buf1(i,j,k)
          enddo
        enddo
      enddo
    endif
  endif

  deallocate( buf1 )
  deallocate( buf2 )

end subroutine copy_gc_f2

subroutine copy_gc_stage( this, dir )

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

end subroutine copy_gc_stage

subroutine acopy_gc_f1( this )

  implicit none

  class( ufield ), intent(inout) :: this

  integer :: idproc, idproc_left, idproc_right, nvp, comm
  integer :: nrp, count, dtype
  integer :: tag = 1, msgid, ierr, i
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:), allocatable :: buf

  idproc = this%pp%getlidproc()
  idproc_left =  idproc - 1
  idproc_right = idproc + 1
  nrp = this%ndp(1)
  nvp = this%pp%getlnvp()
  comm = this%pp%getlgrp()
  dtype = this%pp%getmreal()

  count = this%dim
  allocate( buf(count) )

  if ( this%gc_num(p_upper,1) == 0 ) then
    call write_err( 'Upper guard cells must be set up for deposition' )
  endif
  ! forward message passing
  ! receiver
  if ( idproc > 0 ) then
    call MPI_IRECV( buf, count, dtype, &
      idproc_left, tag, comm, msgid, ierr )
  endif
  ! sender
  if ( idproc < nvp-1 ) then
    call MPI_SEND( this%f1(1,nrp+1), count, dtype, &
      idproc_right, tag, comm, ierr )
  endif
  ! wait receiving finish and add up guard cells
  if ( idproc > 0 ) then
    call MPI_WAIT( msgid, stat, ierr )
    do i = 1, count
      this%f1(i,1) = this%f1(i,1) + buf(i)
    enddo
  else
    do i = 1, count
      this%f1(i,1) = this%f1(i,1) - this%f1(i,0)
      ! this%f1(i,0) = this%f1(i,1)
      ! this%f1(i,0) = this%f1(i,0)
    enddo
  endif

  ! ! backward message passing
  ! ! receiver
  ! if ( idproc < nvp-1 ) then
  !   call MPI_IRECV( this%f1(1,nrp+1), count, dtype, &
  !     idproc_right, tag, comm, msgid, ierr )
  ! endif
  ! ! sender
  ! if ( idproc > 0 ) then
  !   call MPI_SEND( this%f1(1,1), count, dtype, &
  !     idproc_left, tag, comm, ierr )
  ! endif
  ! ! wait receiving finish
  ! if ( idproc < nvp-1 ) then
  !   call MPI_WAIT( msgid, stat, ierr )
  ! endif

  deallocate( buf )

end subroutine acopy_gc_f1

subroutine acopy_gc_f2( this )

  implicit none

  class( ufield ), intent(inout) :: this

  integer :: idproc, idproc_left, idproc_right, nvp, comm
  integer :: nrp, nzp, count, dtype
  integer :: tag = 1, msgid, ierr, i, j
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:,:), allocatable :: buf1, buf2

  idproc = this%pp%getlidproc()
  idproc_left =  idproc - 1
  idproc_right = idproc + 1
  nrp = this%ndp(1)
  nzp = this%ndp(2)
  nvp = this%pp%getlnvp()
  comm = this%pp%getlgrp()
  dtype = this%pp%getmreal()

  count = this%dim * (nzp+1)
  allocate( buf1( this%dim, nzp+1 ) )
  allocate( buf2( this%dim, nzp+1 ) )

  if ( this%gc_num(p_upper,1) == 0 ) then
    call write_err( 'Upper guard cells must be set up for deposition' )
  endif
  ! forward message passing
  ! receiver
  if ( idproc > 0 ) then
    buf1 = 0.0
    call MPI_IRECV( buf1, count, dtype, &
      idproc_left, tag, comm, msgid, ierr )
  endif
  ! sender
  if ( idproc < nvp-1 ) then
    do j = 1, nzp+1
      do i = 1, this%dim
        buf2(i,j) = this%f2(i,nrp+1,j)
      enddo
    enddo
    call MPI_SEND( buf2, count, dtype, &
      idproc_right, tag, comm, ierr )
  endif
  ! wait receiving finish and add up guard cells
  if ( idproc > 0 ) then
    call MPI_WAIT( msgid, stat, ierr )
    do j = 1, nzp+1
      do i = 1, this%dim
        this%f2(i,1,j) = this%f2(i,1,j) + buf1(i,j)
      enddo
    enddo
  else
    do j = 1, nzp+1
      do i = 1, this%dim
        this%f2(i,1,j) = this%f2(i,1,j) - this%f2(i,0,j)
        ! this%f2(i,0,j) = this%f2(i,1,j)
        ! this%f2(i,0,j) = 2.0 * this%f2(i,0,j)
      enddo
    enddo
  endif

  ! backward message passing
  ! receiver
  ! if ( idproc < nvp-1 ) then
  !   buf1 = 0.0
  !   call MPI_IRECV( buf1, count, dtype, &
  !     idproc_right, tag, comm, msgid, ierr )
  ! endif
  ! ! sender
  ! if ( idproc > 0 ) then
  !   do j = 1, nzp+1
  !     do i = 1, this%dim
  !       buf2(i,j) = this%f2(i,1,j)
  !     enddo
  !   enddo
  !   call MPI_SEND( buf2, count, dtype, &
  !     idproc_left, tag, comm, ierr )
  ! endif
  ! ! wait receiving finish
  ! if ( idproc < nvp-1 ) then
  !   call MPI_WAIT( msgid, stat, ierr )
  !   do j = 1, nzp+1
  !     do i = 1, this%dim
  !       this%f2(i,nrp+1,j) = buf1(i,j)
  !     enddo
  !   enddo
  ! endif

  deallocate( buf1 )
  deallocate( buf2 )

end subroutine acopy_gc_f2

subroutine acopy_gc_stage( this )

  implicit none

  class( ufield ), intent(inout) :: this

  integer :: idproc, idproc_next, idproc_last, nstage, stageid, nvp, comm
  integer :: nzp, nrp, count, dtype
  integer :: tag = 1, msgid, ierr, i, j
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:,:), allocatable :: buf1, buf2

  nvp = this%nvp(1)
  nstage = this%pp%getnstage()
  stageid = this%pp%getstageid()
  idproc = this%pp%getidproc()
  idproc_last = idproc - nvp
  idproc_next = idproc + nvp
  nrp = this%ndp(1)
  nzp = this%ndp(2)
  comm = this%pp%getlworld()
  dtype = this%pp%getmreal()

  count = this%dim * (nrp+1)
  allocate( buf1( this%dim, nrp+1 ) )
  allocate( buf2( this%dim, nrp+1 ) )
        
  ! forward message passing
  ! receiver
  if ( stageid > 0 ) then
    buf1 = 0.0
    call MPI_IRECV( buf1, count, dtype, idproc_last, tag, comm, msgid, ierr )
  endif
  ! sender
  if ( stageid < nstage-1 ) then
    do j = 1, nrp+1
      do i = 1, this%dim
        buf2 = this%f2(i,j,nzp+1)
      enddo
    enddo
    call MPI_SEND( buf2, count, dtype, idproc_next, tag, comm, ierr )
  endif
  ! wait receiving finish
  if ( stageid > 0 ) then
    call MPI_WAIT( msgid, stat, ierr )
    do j = 1, nrp+1
      do i = 1, this%dim
        this%f2(i,j,1) = this%f2(i,j,1) + buf1(i,j)
      enddo
    enddo
  endif

  ! backward message passing
  ! receiver
  if ( stageid < nstage-1 ) then
    buf1 = 0.0
    call MPI_IRECV( buf1, count, dtype, idproc_next, tag, comm, msgid, ierr )
  endif
  ! sender
  if ( stageid > 0 ) then
    do j = 1, nrp+1
      do i = 1, this%dim
        buf2(i,j) = this%f2(i,j,1)
      enddo
    enddo
    call MPI_SEND( buf2, count, dtype, idproc_last, tag, comm, ierr )
  endif
  ! wait receiving finish
  if ( stageid < nstage-1 ) then
    call MPI_WAIT( msgid, stat, ierr )
    do j = 1, nrp+1
      do i = 1, this%dim
        this%f2(i,j,nrp+1) = buf1(i,j)
      enddo
    enddo
  endif

  deallocate( buf1 )
  deallocate( buf2 )

end subroutine acopy_gc_stage

subroutine assign_f1( this, that )

  implicit none

  class( ufield ), intent(inout) :: this
  class(*), intent(in) :: that

  integer :: i, j

  select type (that)
    
    type is (real)

      this%f1 = that

    class is (ufield)

      if ( all( this%gc_num(:,1)==that%gc_num(:,1) ) ) then
        this%f1 = that%f1
      else
        call write_err( "guard cells not matched!" )
      endif

    class default

      call write_err( "invalid assignment type!" )

  end select

end subroutine assign_f1

subroutine assign_f2( this, that )

  implicit none

  class( ufield ), intent(inout) :: this
  class(*), intent(in) :: that

  integer :: i, j, k

  select type (that)
    
    type is (real)

      this%f2 = that

    class is (ufield)

      if ( all( this%gc_num==that%gc_num ) ) then
        this%f2 = that%f2
      else
        call write_err( "guard cells not matched!" )
      endif

    class default

      call write_err( "invalid assignment type!" )

  end select

end subroutine assign_f2

subroutine add_f1_binary_dim( a1, a2, a3, dim1, dim2, dim3 )

  implicit none

  class( ufield ), intent(in) :: a1, a2
  class( ufield ), intent(inout) :: a3
  integer, intent(in), dimension(:) :: dim1, dim2, dim3
  
  integer :: ndim, i

  ndim = size( dim1 )

  if ( all(a1%gc_num(:,1)==a2%gc_num(:,1)) &
    .and. all(a1%gc_num(:,1)==a3%gc_num(:,1))  ) then

    do i = 1, ndim
      a3%f1( dim3(i), : ) = a1%f1( dim1(i), : ) + a2%f1( dim2(i), : )
    enddo

  else
    call write_err( "guard cells not matched!" )
  endif

end subroutine add_f1_binary_dim

subroutine add_f1_unitary_dim( a, b, adim, bdim )

  implicit none

  class( ufield ), intent(in) :: a
  class( ufield ), intent(inout) :: b
  integer, intent(in), dimension(:) :: adim, bdim
  
  integer :: ndim, i

  ndim = size( adim )

  if ( all(a%gc_num(:,1)==b%gc_num(:,1)) ) then

    do i = 1, ndim
      b%f1( bdim(i), : ) = b%f1( bdim(i), : ) + a%f1( adim(i), : )
    enddo

  else
    call write_err( "guard cells not matched!" )
  endif

end subroutine add_f1_unitary_dim

subroutine add_f1_binary( a1, a2, a3 )

  implicit none

  class( * ), intent(in) :: a1, a2
  class( ufield ), intent(inout) :: a3
  
  select type (a1)

    type is (real)

      select type (a2)

        type is (real)

          a3%f1 = a1 + a2

        class is (ufield)

          if ( all(a2%gc_num(:,1)==a3%gc_num(:,1)) ) then
            a3%f1 = a1 + a2%f1
          else
            call write_err( "guard cells not matched!" )
          endif

        class default
          call write_err( "Invalid argument type!" )
      end select

    class is (ufield)

      select type (a2)

        type is (real)

          if ( all(a1%gc_num(:,1)==a3%gc_num(:,1)) ) then
            a3%f1 = a1%f1 + a2
          else
            call write_err( "guard cells not matched!" )
          endif

        class is (ufield)

          if ( all(a1%gc_num(:,1)==a3%gc_num(:,1)) .and. &
               all(a2%gc_num(:,1)==a3%gc_num(:,1)) ) then
            a3%f1 = a1%f1 + a2%f1
          else
            call write_err( "guard cells not matched!" )
          endif

        class default
          call write_err( "Invalid argument type!" )
      end select

    class default
      call write_err( "Invalid argument type!" )
  end select

end subroutine add_f1_binary

subroutine add_f1_unitary( a, b )

  implicit none

  class( * ), intent(in) :: a
  class( ufield ), intent(inout) :: b
  
  select type (a)

    type is (real)

      b%f1 = b%f1 + a

    class is (ufield)

      if ( all(a%gc_num(:,1)==b%gc_num(:,1)) ) then
        b%f1 = b%f1 + a%f1
      else
        call write_err( "guard cells not matched!" )
      endif

    class default
      call write_err( "Invalid argument type!" )
  end select

end subroutine add_f1_unitary

subroutine sub_f1_binary_dim( a1, a2, a3, dim1, dim2, dim3 )

  implicit none

  class( ufield ), intent(in) :: a1, a2
  class( ufield ), intent(inout) :: a3
  integer, intent(in), dimension(:) :: dim1, dim2, dim3
  
  integer :: ndim, i

  ndim = size( dim1 )

  if ( all(a1%gc_num(:,1)==a2%gc_num(:,1)) &
    .and. all(a1%gc_num(:,1)==a3%gc_num(:,1))  ) then

    do i = 1, ndim
      a3%f1( dim3(i), : ) = a1%f1( dim1(i), : ) - a2%f1( dim2(i), : )
    enddo

  else
    call write_err( "guard cells not matched!" )
  endif

end subroutine sub_f1_binary_dim

subroutine sub_f1_unitary_dim( a, b, adim, bdim )

  implicit none

  class( ufield ), intent(in) :: a
  class( ufield ), intent(inout) :: b
  integer, intent(in), dimension(:) :: adim, bdim
  
  integer :: ndim, i

  ndim = size( adim )

  if ( all(a%gc_num(:,1)==b%gc_num(:,1)) ) then

    do i = 1, ndim
      b%f1( bdim(i), : ) = b%f1( bdim(i), : ) - a%f1( adim(i), : )
    enddo

  else
    call write_err( "guard cells not matched!" )
  endif

end subroutine sub_f1_unitary_dim

subroutine sub_f1_binary( a1, a2, a3 )

  implicit none

  class( * ), intent(in) :: a1, a2
  class( ufield ), intent(inout) :: a3
  
  select type (a1)

    type is (real)

      select type (a2)

        type is (real)

          a3%f1 = a1 - a2

        class is (ufield)

          if ( all(a2%gc_num(:,1)==a3%gc_num(:,1)) ) then
            a3%f1 = a1 - a2%f1
          else
            call write_err( "guard cells not matched!" )
          endif

        class default
          call write_err( "Invalid argument type!" )
      end select

    class is (ufield)

      select type (a2)

        type is (real)

          if ( all(a1%gc_num(:,1)==a3%gc_num(:,1)) ) then
            a3%f1 = a1%f1 - a2
          else
            call write_err( "guard cells not matched!" )
          endif

        class is (ufield)

          if ( all(a1%gc_num(:,1)==a3%gc_num(:,1)) .and. &
               all(a2%gc_num(:,1)==a3%gc_num(:,1)) ) then
            a3%f1 = a1%f1 - a2%f1
          else
            call write_err( "guard cells not matched!" )
          endif

        class default
          call write_err( "Invalid argument type!" )
      end select

    class default
      call write_err( "Invalid argument type!" )
  end select

end subroutine sub_f1_binary

subroutine sub_f1_unitary( a, b )

  implicit none

  class( * ), intent(in) :: a
  class( ufield ), intent(inout) :: b
  
  select type (a)

    type is (real)

      b%f1 = b%f1 - a

    class is (ufield)

      if ( all(a%gc_num(:,1)==b%gc_num(:,1)) ) then
        b%f1 = b%f1 - a%f1
      else
        call write_err( "guard cells not matched!" )
      endif

    class default
      call write_err( "Invalid argument type!" )
  end select

end subroutine sub_f1_unitary

subroutine dot_f1_unitary( a, b )

  implicit none

  class( * ), intent(in) :: a
  class( ufield ), intent(inout) :: b
  
  select type (a)

    type is (real)

      b%f1 = b%f1 * a

    class is (ufield)

      if ( all(a%gc_num(:,1)==b%gc_num(:,1)) ) then
        b%f1 = b%f1 * a%f1
      else
        call write_err( "guard cells not matched!" )
      endif

    class default
      call write_err( "Invalid argument type!" )
  end select

end subroutine dot_f1_unitary

subroutine add_f2_binary( a1, a2, a3 )

  implicit none

  class( * ), intent(in) :: a1, a2
  class( ufield ), intent(inout) :: a3
  
  select type (a1)

    type is (real)

      select type (a2)

        type is (real)

          a3%f2 = a1 + a2

        class is (ufield)

          if ( all(a2%gc_num==a3%gc_num) ) then
            a3%f2 = a1 + a2%f2
          else
            call write_err( "guard cells not matched!" )
          endif

        class default
          call write_err( "Invalid argument type!" )
      end select

    class is (ufield)

      select type (a2)

        type is (real)

          if ( all(a1%gc_num==a3%gc_num) ) then
            a3%f2 = a1%f2 + a2
          else
            call write_err( "guard cells not matched!" )
          endif

        class is (ufield)

          if ( all(a1%gc_num==a3%gc_num) .and. &
               all(a2%gc_num==a3%gc_num) ) then
            a3%f2 = a1%f2 + a2%f2
          else
            call write_err( "guard cells not matched!" )
          endif

        class default
          call write_err( "Invalid argument type!" )
      end select

    class default
      call write_err( "Invalid argument type!" )
  end select

end subroutine add_f2_binary

subroutine add_f2_unitary( a, b )

  implicit none

  class( * ), intent(in) :: a
  class( ufield ), intent(inout) :: b
  
  select type (a)

    type is (real)

      b%f2 = b%f2 + a

    class is (ufield)

      if ( all(a%gc_num==b%gc_num) ) then
        b%f2 = b%f2 + a%f2
      else
        call write_err( "guard cells not matched!" )
      endif

    class default
      call write_err( "Invalid argument type!" )
  end select

end subroutine add_f2_unitary

subroutine sub_f2_binary( a1, a2, a3 )

  implicit none

  class( * ), intent(in) :: a1, a2
  class( ufield ), intent(inout) :: a3
  
  select type (a1)

    type is (real)

      select type (a2)

        type is (real)

          a3%f2 = a1 - a2

        class is (ufield)

          if ( all(a2%gc_num==a3%gc_num) ) then
            a3%f2 = a1 - a2%f2
          else
            call write_err( "guard cells not matched!" )
          endif

        class default
          call write_err( "Invalid argument type!" )
      end select

    class is (ufield)

      select type (a2)

        type is (real)

          if ( all(a1%gc_num==a3%gc_num) ) then
            a3%f2 = a1%f2 - a2
          else
            call write_err( "guard cells not matched!" )
          endif

        class is (ufield)

          if ( all(a1%gc_num==a3%gc_num) .and. &
               all(a2%gc_num==a3%gc_num) ) then
            a3%f2 = a1%f2 - a2%f2
          else
            call write_err( "guard cells not matched!" )
          endif

        class default
          call write_err( "Invalid argument type!" )
      end select

    class default
      call write_err( "Invalid argument type!" )
  end select

end subroutine sub_f2_binary

subroutine sub_f2_unitary( a, b )

  implicit none

  class( * ), intent(in) :: a
  class( ufield ), intent(inout) :: b
  
  select type (a)

    type is (real)

      b%f2 = b%f2 - a

    class is (ufield)

      if ( all(a%gc_num==b%gc_num) ) then
        b%f2 = b%f2- a%f2
      else
        call write_err( "guard cells not matched!" )
      endif

    class default
      call write_err( "Invalid argument type!" )
  end select

end subroutine sub_f2_unitary

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

function has2d( this )

  implicit none

  class( ufield ), intent(in) :: this
  logical :: has2d

  has2d = this%has_2d

end function has2d

! function add_f1_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( ufield ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( ufield ), allocatable :: a3

!   integer :: i, j, dim_min

!   allocate( a3 )
!   ! call a3%new(a1)

!   select type (a2)

!   type is (real)

!     call a3%new(a1)
!     a3%f1 = a1%f1 + a2

!   class is (ufield)

!     if ( all(a1%gc_num(:,1)==a2%gc_num(:,1)) ) then

!       dim_min = min( a1%get_dim(), a2%get_dim() )

!       if ( a1%get_dim() > a2%get_dim() ) then
!         call a3%new(a1)
!       else
!         call a3%new(a2)
!       endif

!       do i = 1, dim_min
!         a3%f1(i,:) = a1%f1(i,:) + a2%f1(i,:)
!       enddo

!       if ( a1%get_dim() > a2%get_dim() ) then
!         do i = dim_min+1, a1%dim
!           a3%f1(i,:) = a1%f1(i,:)
!         enddo
!       elseif ( a1%get_dim() < a2%get_dim() ) then
!         do i = dim_min+1, a2%dim
!           a3%f1(i,:) = a2%f1(i,:)
!         enddo
!       endif      

!     else
!       call write_err( "guard cells not matched!" )
!     endif

!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function add_f1_v1

! function add_f2_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( ufield ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( ufield ), allocatable :: a3

!   integer :: i, j, k

!   allocate( a3 )
!   call a3%new(a1)

!   select type (a2)
!   type is (real)
!     a3%f2 = a1%f2 + a2
!   class is (ufield)
!     if ( all( a1%gc_num==a2%gc_num ) ) then
!       a3%f2 = a1%f2 + a2%f2
!     else
!       call write_err( "guard cells not matched!" )
!     endif
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function add_f2_v1

! function add_f1_v2( a2, a1 ) result( a3 )

!   implicit none

!   real, intent(in) :: a2
!   class( ufield ), intent(in) :: a1
!   class( ufield ), allocatable :: a3

!   integer :: i, j

!   allocate( a3 )
!   call a3%new(a1)

!   a3%f1 = a1%f1 + a2
    
! end function add_f1_v2

! function add_f2_v2( a2, a1 ) result( a3 )

!   implicit none

!   real, intent(in) :: a2
!   class( ufield ), intent(in) :: a1
!   class( ufield ), allocatable :: a3

!   integer :: i, j, k

!   allocate( a3 )
!   call a3%new(a1)

!   a3%f2 = a1%f2 + a2

! end function add_f2_v2

! function sub_f1_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( ufield ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( ufield ), allocatable :: a3

!   integer :: i, j

!   allocate( a3 )
!   call a3%new(a1)

!   select type (a2)
!   type is (real)
!     a3%f1 = a1%f1 - a2
!   class is (ufield)
!     if ( all(a1%gc_num(:,1)==a2%gc_num(:,1)) ) then
!       a3%f1 = a1%f1 - a2%f1
!     else
!       call write_err( "guard cells not matched!" )
!     endif
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function sub_f1_v1

! function sub_f2_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( ufield ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( ufield ), allocatable :: a3

!   integer :: i, j, k

!   allocate( a3 )
!   call a3%new(a1)

!   select type (a2)
!   type is (real)
!     a3%f2 = a1%f2 - a2
!   class is (ufield)
!     if ( all(a1%gc_num==a2%gc_num) ) then
!       a3%f2 = a1%f2 - a2%f2
!     else
!       call write_err( "guard cells not matched!" )
!     endif
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function sub_f2_v1

! function sub_f1_v2( a2, a1 ) result( a3 )

!   implicit none

!   real, intent(in) :: a2
!   class( ufield ), intent(in) :: a1
!   class( ufield ), allocatable :: a3

!   integer :: i, j

!   allocate( a3 )
!   call a3%new(a1)

!   a3%f1 = a2 - a1%f1
    
! end function sub_f1_v2

! function sub_f2_v2( a2, a1 ) result( a3 )

!   implicit none

!   real, intent(in) :: a2
!   class( ufield ), intent(in) :: a1
!   class( ufield ), allocatable :: a3

!   integer :: i, j, k

!   allocate( a3 )
!   call a3%new(a1)

!   a3%f2 = a2 - a1%f2
    
! end function sub_f2_v2

! function dot_f1_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( ufield ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( ufield ), allocatable :: a3

!   integer :: i, j

!   allocate( a3 )
!   call a3%new(a1)

!   select type (a2)
!   type is (real)
!     a3%f1 = a1%f1 * a2
!   class is (ufield)
!     if ( all(a1%gc_num(:,1)==a2%gc_num(:,1)) ) then
!       a3%f1 = a1%f1 * a2%f1
!     else
!       call write_err( "guard cells not matched!" )
!     endif
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function dot_f1_v1

! function dot_f2_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( ufield ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( ufield ), allocatable :: a3

!   integer :: i, j, k

!   allocate( a3 )
!   call a3%new(a1)

!   select type (a2)
!   type is (real)
!     a3%f2 = a1%f2 * a2
!   class is (ufield)
!     if ( all(a1%gc_num==a2%gc_num) ) then
!       a3%f2 = a1%f2 * a2%f2
!     else
!       call write_err( "guard cells not matched!" )
!     endif
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function dot_f2_v1

! function dot_f1_v2( a2, a1 ) result( a3 )

!   implicit none

!   real, intent(in) :: a2
!   class( ufield ), intent(in) :: a1
!   class( ufield ), allocatable :: a3

!   integer :: i, j

!   allocate( a3 )
!   call a3%new(a1)

!   a3%f1 = a1%f1 * a2
    
! end function dot_f1_v2

! function dot_f2_v2( a2, a1 ) result( a3 )

!   implicit none

!   real, intent(in) :: a2
!   class( ufield ), intent(in) :: a1
!   class( ufield ), allocatable :: a3

!   integer :: i, j, k

!   allocate( a3 )
!   call a3%new(a1)

!   a3%f2 = a1%f2 * a2
    
! end function dot_f2_v2

end module ufield_class