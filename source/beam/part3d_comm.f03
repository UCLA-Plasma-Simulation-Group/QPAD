module part3d_comm

use sysutil
use param
use mpi
use parallel_module
use options_class
use ufield_class
use part3d_class

implicit none

save

! direction indicators
! outward/inward: positive/negative direction in r
! forward/backward: positive/negative direction in xi
integer, parameter :: p_self = 0, &
                      p_iwd = 1, &
                      p_owd = 2, &
                      p_bwd = 3, &
                      p_fwd = 4

integer, parameter :: iter_max = 1000

! send and receive buffer for MPI
real, dimension(:,:,:), pointer :: send_buf => null()
real, dimension(:,:,:), pointer :: recv_buf => null()

! index of holes
integer(kind=LG), dimension(:), pointer :: ihole => null()

! count of MPI buffers
integer, dimension(4) :: send_cnt = 0
integer, dimension(4) :: recv_cnt = 0

! max dimension of particle coordinates
integer :: dim_max = 6

! current buffer size, same for send and receive buffers
integer :: buf_size = 0

! self and neighbor processor id
integer, dimension(0:4) :: pid

! indicator of if neighbors are physical boundaries
logical, dimension(4) :: phys_bnd

! edges of current MPI node
real, dimension(2) :: redge, zedge

interface init_part3d_comm
  module procedure init_part3d_comm
end interface

interface end_part3d_comm
  module procedure end_part3d_comm
end interface

interface set_part3d_comm
  module procedure set_part3d_comm
end interface

interface move_part3d_comm
  module procedure move_part3d_comm
end interface

public :: init_part3d_comm, end_part3d_comm, set_part3d_comm, move_part3d_comm

contains

subroutine set_part3d_comm( part_dim, npmax )

  implicit none

  integer, intent(in) :: part_dim
  integer(kind=LG), intent(in) :: npmax

  real :: buf_ratio = 0.01 

  dim_max  = max( part_dim, dim_max )
  buf_size = max( int(npmax * buf_ratio), buf_size )

end subroutine set_part3d_comm

subroutine init_part3d_comm( opts )

  implicit none

  type( options ), intent(in) :: opts

  integer :: nvp, nst
  integer, dimension(2) :: noff, ndp
  real :: dr, dz

  allocate( send_buf( dim_max, buf_size, 2 ) )
  allocate( recv_buf( dim_max, buf_size, 2 ) )
  allocate( ihole( buf_size ) )

  nvp = num_procs_loc()
  nst = num_stages()

  ! set neighbor processors' id
  pid( p_self ) = id_proc()
  pid( p_iwd )  = pid( p_self ) - 1
  pid( p_owd )  = pid( p_self ) + 1
  pid( p_bwd )  = pid( p_self ) - nvp
  pid( p_fwd )  = pid( p_self ) + nvp

  ! check if the neighbors are physical boundaries
  phys_bnd = get_phys_bnd( pid(p_self), nvp, nst )

  ! set edges of MPI node
  noff = opts%get_noff()
  ndp  = opts%get_ndp()
  dr   = opts%get_dr()
  dz   = opts%get_dxi()
  
  if ( phys_bnd(p_iwd) ) then
    redge(p_lower) = 0.0
  else
    redge(p_lower) = noff(1) * dr
  endif
  redge(p_upper) = redge(p_lower) + ndp(1) * dr
  zedge(p_lower) = noff(2) * dz
  zedge(p_upper) = zedge(p_lower) + ndp(2) * dz

end subroutine init_part3d_comm

subroutine end_part3d_comm()

  implicit none

  if ( associated( send_buf ) ) deallocate( send_buf )
  if ( associated( recv_buf ) ) deallocate( recv_buf )
  if ( associated( ihole ) ) deallocate( ihole )

end subroutine end_part3d_comm

subroutine move_part3d_comm( part, rtag, stag, id )

  implicit none

  class( part3d ), intent(inout) :: part
  integer, intent(in) :: rtag, stag
  integer, intent(inout) :: id

  integer :: bsize, rsize, ssize, ierr, part_dim, iter
  integer :: loc_grp, world
  integer, dimension(MPI_STATUS_SIZE) :: istat
  integer, dimension(2) :: sid, rid
  integer, dimension(4) :: scnt_max

  world    = comm_world()
  loc_grp  = comm_loc()
  part_dim = part%part_dim
  bsize    = buf_size * dim_max

  call start_tprof( 'move 3D particles' )

  ! receive particles from backward node (pipeline communication)
  if ( .not. phys_bnd(p_bwd) ) then

    call mpi_recv( recv_buf(:,:,p_lower), bsize, p_dtype_real, pid(p_bwd), &
      rtag, world, istat, ierr )

    call mpi_get_count( istat, p_dtype_real, rsize, ierr )
    recv_cnt(p_bwd) = rsize / part_dim

    call unpack_particles( part, zedge, recv_buf(:,:,p_lower), recv_cnt(p_bwd), p_bwd, &
      part%pbuff, send_cnt(p_fwd), p_fwd )

  endif
  recv_cnt(p_bwd) = 0

  ! send particles to forward node (pipeline communication)
  if ( .not. phys_bnd(p_fwd) ) then

    call pack_particles( part, zedge, part%pbuff, ihole, send_cnt(p_fwd), p_fwd )
    ssize = send_cnt(p_fwd) * part_dim

    call mpi_isend( part%pbuff, ssize, p_dtype_real, pid(p_fwd), &
      stag, world, id, ierr )

  else
    id = MPI_REQUEST_NULL
  endif
  send_cnt(p_fwd) = 0

  ! iteration of moving particles in r direction
  iter = 0
  do

    if ( iter == iter_max ) then
      call write_err( 'Iteration overflow' )
    endif

    iter = iter + 1

    ! receive from inward and unpack particles
    if ( .not. phys_bnd(p_iwd) ) then
      call mpi_irecv( recv_buf(:,:,p_lower), bsize, p_dtype_real, pid(p_iwd), &
        iter, world, rid(1), ierr )
    else
      rid(1) = MPI_REQUEST_NULL
    endif

    ! receive from outward and unpack particles
    if ( .not. phys_bnd(p_owd) ) then
      call mpi_irecv( recv_buf(:,:,p_upper), bsize, p_dtype_real, pid(p_owd), &
        iter, world, rid(2), ierr )
    else
      rid(2) = MPI_REQUEST_NULL
    endif

    ! pack particles to inward/outward send buffers
    call pack_particles( part, redge, send_buf(:,:,p_lower), ihole, &
      send_cnt(p_iwd), p_iwd )
    call pack_particles( part, redge, send_buf(:,:,p_upper), ihole, &
      send_cnt(p_owd), p_owd )

    ! send inward
    if ( .not. phys_bnd(p_iwd) ) then

      ssize = send_cnt(p_iwd) * part_dim
      call mpi_isend( send_buf(:,:,p_lower), ssize, p_dtype_real, pid(p_iwd), &
        iter, world, sid(1), ierr )
      
    else
      sid(1) = MPI_REQUEST_NULL
    endif
    send_cnt(p_iwd) = 0
    
    ! send outward
    if ( .not. phys_bnd(p_owd) ) then

      ssize = send_cnt(p_owd) * part_dim
      call mpi_isend( send_buf(:,:,p_upper), ssize, p_dtype_real, pid(p_owd), &
        iter, world, sid(2), ierr )

    else
      sid(2) = MPI_REQUEST_NULL
    endif
    send_cnt(p_owd) = 0

    ! wait receiving finish
    call mpi_wait( rid(1), istat, ierr )
    call mpi_get_count( istat, p_dtype_real, rsize, ierr )
    recv_cnt(p_iwd) = rsize / part_dim

    call mpi_wait( rid(2), istat, ierr )
    call mpi_get_count( istat, p_dtype_real, rsize, ierr )
    recv_cnt(p_owd) = rsize / part_dim

    ! wait sending finsh
    call mpi_wait( sid(1), istat, ierr )
    call mpi_wait( sid(2), istat, ierr )

    ! unpack particles in inward receive buffer
    if ( .not. phys_bnd(p_iwd) ) then
      call unpack_particles( part, redge, recv_buf(:,:,p_lower), recv_cnt(p_iwd), p_iwd, &
        send_buf(:,:,p_upper), send_cnt(p_owd), p_owd )
    endif

    ! unpack particles in outward receive buffer
    if ( .not. phys_bnd(p_owd) ) then
      call unpack_particles( part, redge, recv_buf(:,:,p_upper), recv_cnt(p_owd), p_owd, &
        send_buf(:,:,p_lower), send_cnt(p_iwd), p_iwd )
    endif

    ! check if need move particles further
    call mpi_allreduce( send_cnt, scnt_max, 4, p_dtype_int, MPI_MAX, loc_grp, ierr )
    if ( scnt_max(p_iwd) == 0 .and. scnt_max(p_owd) == 0 ) then
      exit
    endif

  enddo

  call stop_tprof( 'move 3D particles' )

end subroutine move_part3d_comm

! -----------------------------------------------------------------------------
! Copy the particles in the receive buffers belonging to this node to the
! particle array and those not belong to this node to other send buffers.
! This routine only moves particles inside current stage.
! -----------------------------------------------------------------------------
subroutine unpack_particles( part, edge, rbuf, rcnt, src, sbuf, scnt, des )

  implicit none

  class(part3d), intent(inout) :: part
  real, intent(in), dimension(2) :: edge
  real, intent(in), dimension(:,:) :: rbuf
  real, intent(inout), dimension(:,:) :: sbuf
  integer, intent(inout) :: scnt, rcnt
  integer, intent(in) :: src, des

  integer :: i, npp, part_dim, stay_cnt
  real, dimension(3) :: x
  logical :: has_spin

  npp      = part%npp
  part_dim = part%part_dim
  has_spin = part%has_spin

  if ( npp + rcnt > part%npmax ) then
    call write_err( 'Particle overflow' )
    ! TODO: resize part3d memory
  endif

  stay_cnt = 0

  do i = 1, rcnt

    x = rbuf( 1:3, i )

    ! check if particle belongs to here
    if ( goto_here( x, edge, src ) ) then

      stay_cnt = stay_cnt + 1

      part%x( 1:3, npp + stay_cnt ) = rbuf( 1:3, i )
      part%p( 1:3, npp + stay_cnt ) = rbuf( 4:6, i )
      part%q( npp + stay_cnt )      = rbuf(   7, i )
      if ( has_spin ) then
        part%s( 1:3, npp + stay_cnt ) = rbuf( 8:10, i )
      endif

      rcnt = rcnt - 1

    ! check if particle goes to the destination neighbor node
    elseif ( goto_des( x, edge, des ) ) then

      scnt = scnt + 1
      rcnt = rcnt - 1

      sbuf( 1:part_dim, scnt ) = rbuf( 1:part_dim, i )

    ! particle goes back to where it comes from, should not happen
    elseif ( goto_des( x, edge, src ) ) then

      call write_err( 'There is returning particle in local particle &
        &communication! ' )

    endif

  enddo

  ! check if all the particles in the receive buffer have been unpacked
  if ( rcnt /= 0 ) then
    call write_err( 'There are unpacked particles!' )
  endif

  part%npp = npp + stay_cnt

end subroutine unpack_particles

subroutine pack_particles( part, edge, sbuf, ihole, scnt, des )

  implicit none

  class(part3d), intent(inout) :: part
  real, intent(in), dimension(2) :: edge
  real, intent(inout), dimension(:,:) :: sbuf
  integer(kind=LG), intent(inout), dimension(:) :: ihole
  integer, intent(inout) :: scnt
  integer, intent(in) :: des

  integer :: i, npp, go_cnt, part_dim, buf_size
  real, dimension(3) :: x
  logical :: has_spin

  npp = part%npp
  part_dim = part%part_dim
  has_spin = part%has_spin
  buf_size = size( sbuf, 2 )

  go_cnt = 0

  do i = 1, npp

    x = part%x(:,i)

    if ( goto_des( x, edge, des ) ) then

      if ( scnt >= buf_size ) then
        call write_err( '3D particle MPI buffer overflow!' )
      endif

      scnt = scnt + 1
      go_cnt = go_cnt + 1
      ihole(go_cnt) = i

      sbuf( 1:3, scnt ) = part%x(1:3,i)
      sbuf( 4:6, scnt ) = part%p(1:3,i)
      sbuf(   7, scnt ) = part%q(i)
      if ( has_spin ) then
        sbuf( 8:10, scnt ) = part%s(1:3,i)
      endif

    endif
  enddo

  ! fill the holes inversely
  do i = go_cnt, 1, -1

    part%x( 1:3, ihole(i) ) = part%x( 1:3, npp )
    part%p( 1:3, ihole(i) ) = part%p( 1:3, npp )
    part%q( ihole(i) )      = part%q(npp)
    if ( has_spin ) then
      part%s( 1:3, ihole(i) ) = part%s( 1:3, npp )
    endif
    
    npp = npp - 1

  enddo

  part%npp = npp

end subroutine pack_particles

subroutine resize_buf( ratio )
! this subroutine is currently not in use.

  implicit none

  real, intent(in), optional :: ratio

  real :: r = 2.0
  real, dimension(:,:,:), allocatable :: buf_tmp
  integer :: i, j

  if ( present(ratio) ) r = ratio

  allocate( buf_tmp( dim_max, buf_size, 4 ) )
  buf_size = int( buf_size * r )

  ! resize receive buffer
  do j = 1, 4
    do i = 1, recv_cnt(j)
      buf_tmp(:,i,j) = recv_buf(:,i,j)
    enddo
  enddo
  deallocate( recv_buf )
  allocate( recv_buf( dim_max, buf_size, 4 ) )
  do j = 1, 4
    do i = 1, recv_cnt(j)
      recv_buf(:,i,j) = buf_tmp(:,i,j)
    enddo
  enddo

  ! resize send buffer
  do j = 1, 4
    do i = 1, send_cnt(j)
      buf_tmp(:,i,j) = send_buf(:,i,j)
    enddo
  enddo
  deallocate( send_buf )
  allocate( send_buf( dim_max, buf_size, 4 ) )
  do j = 1, 4
    do i = 1, recv_cnt(j)
      send_buf(:,i,j) = buf_tmp(:,i,j)
    enddo
  enddo

  ! resize index of holes
  deallocate( ihole )
  allocate( ihole( buf_size ) )

  deallocate( buf_tmp )

end subroutine resize_buf

function get_phys_bnd( proc_id, nvp, nst ) result( res )

  implicit none

  integer, intent(in) :: proc_id, nvp, nst
  logical, dimension(4) :: res

  integer :: loc_id, st_id

  res = .false.

  loc_id = mod( proc_id, nvp )
  st_id = int( proc_id / nvp )

  if ( loc_id == 0 )     res( p_iwd ) = .true.
  if ( loc_id == nvp-1 ) res( p_owd ) = .true.
  if ( st_id == 0 )      res( p_bwd ) = .true.
  if ( st_id == nst-1 )  res( p_fwd ) = .true.

end function get_phys_bnd

function goto_here( x, edge, src )

  implicit none

  real, intent(in), dimension(3) :: x
  real, intent(in), dimension(2) :: edge
  integer, intent(in) :: src
  logical :: goto_here

  real :: pos

  goto_here = .false.

  select case ( src )
  case ( p_fwd, p_bwd )
    pos = x(3)
    goto_here = pos >= edge(p_lower) .and. pos < edge(p_upper)
  case ( p_iwd, p_owd )
    pos = sqrt( x(1)*x(1) + x(2)*x(2) )
    goto_here = pos >= edge(p_lower) .and. pos < edge(p_upper)
  end select

end function goto_here

function goto_des( x, edge, des )

  implicit none

  real, intent(in), dimension(3) :: x
  real, intent(in), dimension(2) :: edge
  integer, intent(in) :: des
  logical :: goto_des

  real :: pos
  
  goto_des = .false.

  select case ( des )
  case ( p_fwd )
    pos = x(3)
    goto_des = pos >= edge(p_upper)
  case ( p_bwd )
    pos = x(3)
    goto_des = pos < edge(p_lower)
  case ( p_iwd )
    pos = sqrt( x(1)*x(1) + x(2)*x(2) )
    goto_des = pos < edge(p_lower)
  case ( p_owd )
    pos = sqrt( x(1)*x(1) + x(2)*x(2) )
    goto_des = pos >= edge(p_upper)
  end select

end function goto_des

end module part3d_comm