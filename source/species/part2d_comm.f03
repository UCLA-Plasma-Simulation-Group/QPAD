module part2d_comm

use sysutil_module
use param
use mpi
use parallel_module
use options_class
use ufield_class
use part2d_class

implicit none

save

! direction indicators
! outward/inward: positive/negative direction in r
integer, parameter :: p_self = 0, &
                      p_iwd = 1, &
                      p_owd = 2

integer, parameter :: iter_max = 1000

! send and receive buffer for MPI
real, dimension(:), pointer :: send_buf_lower => null()
real, dimension(:), pointer :: send_buf_upper => null()
real, dimension(:), pointer :: recv_buf_lower => null()
real, dimension(:), pointer :: recv_buf_upper => null()

! index of holes
integer(kind=LG), dimension(:), pointer :: ihole => null()

! count of MPI buffers
integer, dimension(2) :: send_cnt = 0
integer, dimension(2) :: recv_cnt = 0

! max dimension of particle coordinates
integer :: dim_max = 8

! current buffer size, same for send and receive buffers
integer :: buf_size = 0

! self and neighbor processor id
integer, dimension(0:2) :: pid

! indicator of if neighbors are physical boundaries
logical, dimension(2) :: phys_bnd

! edges of current MPI node
real, dimension(2) :: redge

interface init_part2d_comm
  module procedure init_part2d_comm
end interface

interface end_part2d_comm
  module procedure end_part2d_comm
end interface

interface set_part2d_comm
  module procedure set_part2d_comm
end interface

interface move_part2d_comm
  module procedure move_part2d_comm
end interface

public :: init_part2d_comm, end_part2d_comm, set_part2d_comm, move_part2d_comm

contains

subroutine set_part2d_comm( part_dim, npmax )

  implicit none

  integer, intent(in) :: part_dim
  integer, intent(in) :: npmax

  real :: buf_ratio = 0.2

  dim_max  = max( part_dim, dim_max )
  buf_size = max( int(npmax * buf_ratio), buf_size )

end subroutine set_part2d_comm

subroutine init_part2d_comm( opts )

  implicit none

  class( options), intent(in) :: opts

  integer :: nvp, noff, nrp
  real :: dr

  allocate( send_buf_lower( dim_max * buf_size ) )
  allocate( send_buf_upper( dim_max * buf_size ) )
  allocate( recv_buf_lower( dim_max * buf_size ) )
  allocate( recv_buf_upper( dim_max * buf_size ) )
  allocate( ihole( buf_size ) )

  nvp = num_procs_loc()

  ! set neighbor processors' id
  pid( p_self ) = id_proc()
  pid( p_iwd )  = pid( p_self ) - 1
  pid( p_owd )  = pid( p_self ) + 1

  ! check if the neighbors are physical boundaries
  phys_bnd = get_phys_bnd( pid(p_self), nvp )

  ! set edges of MPI node
  noff = opts%get_noff(1)
  nrp  = opts%get_ndp(1)
  dr   = opts%get_dr()
  
  if ( phys_bnd(p_iwd) ) then
    redge(p_lower) = 0.0
  else
    redge(p_lower) = noff * dr
  endif
  redge(p_upper) = redge(p_lower) + nrp * dr

end subroutine init_part2d_comm

subroutine end_part2d_comm()

  implicit none

  if ( associated( send_buf_lower ) ) deallocate( send_buf_lower )
  if ( associated( send_buf_upper ) ) deallocate( send_buf_upper )
  if ( associated( recv_buf_lower ) ) deallocate( recv_buf_lower )
  if ( associated( recv_buf_upper ) ) deallocate( recv_buf_upper )
  if ( associated( ihole ) ) deallocate( ihole )

end subroutine end_part2d_comm

subroutine move_part2d_comm( part )

  implicit none

  class( part2d ), intent(inout) :: part

  integer :: bsize, rsize, ssize, ierr, part_dim, iter, loc_grp, world
  integer, dimension(MPI_STATUS_SIZE) :: istat
  integer, dimension(2) :: sid, rid
  integer, dimension(2) :: scnt_max

  character(len=18), save :: sname = 'move_part2d_comm'

  world    = comm_world()
  loc_grp  = comm_loc()
  part_dim = part%part_dim
  bsize    = buf_size * dim_max

  call write_dbg("part2d_comm", sname, 2, 'starts')
  call start_tprof( 'move 2D particles' )

  ! pack particles to inward/outward send buffers
  call pack_particles( part, redge, send_buf_lower, ihole, &
    send_cnt(p_iwd), p_iwd )
  call pack_particles( part, redge, send_buf_upper, ihole, &
    send_cnt(p_owd), p_owd )

  ! iteration of moving particles in r direction
  iter = 0
  do

    if ( iter == iter_max ) then
      call write_err( 'Iteration overflow' )
    endif

    iter = iter + 1

    ! send inward
    if ( .not. phys_bnd(p_iwd) ) then

      ssize = send_cnt(p_iwd) * part_dim
      call mpi_isend( send_buf_lower, ssize, p_dtype_real, pid(p_iwd), &
        iter, world, sid(1), ierr )
      
    else
      sid(1) = MPI_REQUEST_NULL
    endif
    
    ! send outward
    if ( .not. phys_bnd(p_owd) ) then

      ssize = send_cnt(p_owd) * part_dim
      call mpi_isend( send_buf_upper, ssize, p_dtype_real, pid(p_owd), &
        iter, world, sid(2), ierr )

    else
      sid(2) = MPI_REQUEST_NULL
    endif

    ! receive from inward and unpack particles
    if ( .not. phys_bnd(p_iwd) ) then
      call mpi_irecv( recv_buf_lower, bsize, p_dtype_real, pid(p_iwd), &
        iter, world, rid(1), ierr )
    else
      rid(1) = MPI_REQUEST_NULL
    endif

    ! receive from outward and unpack particles
    if ( .not. phys_bnd(p_owd) ) then
      call mpi_irecv( recv_buf_upper, bsize, p_dtype_real, pid(p_owd), &
        iter, world, rid(2), ierr )
    else
      rid(2) = MPI_REQUEST_NULL
    endif

    ! wait receiving finish
    call mpi_wait( rid(1), istat, ierr )
    call mpi_get_count( istat, p_dtype_real, rsize, ierr )
    recv_cnt(p_iwd) = rsize / part_dim

    call mpi_wait( rid(2), istat, ierr )
    call mpi_get_count( istat, p_dtype_real, rsize, ierr )
    recv_cnt(p_owd) = rsize / part_dim

    ! wait sending finsh
    call mpi_wait( sid(1), istat, ierr ); send_cnt(p_iwd) = 0
    call mpi_wait( sid(2), istat, ierr ); send_cnt(p_owd) = 0

    ! unpack particles in inward receive buffer
    if ( .not. phys_bnd(p_iwd) ) then
      call unpack_relay_particles( part, redge, recv_buf_lower, recv_cnt(p_iwd), p_iwd, &
        send_buf_upper, send_cnt(p_owd), p_owd )
    endif

    ! unpack particles in outward receive buffer
    if ( .not. phys_bnd(p_owd) ) then
      call unpack_relay_particles( part, redge, recv_buf_upper, recv_cnt(p_owd), p_owd, &
        send_buf_lower, send_cnt(p_iwd), p_iwd )
    endif

    ! check if need move particles further
    call mpi_allreduce( send_cnt, scnt_max, 2, p_dtype_int, MPI_MAX, loc_grp, ierr )
    if ( scnt_max(p_iwd) == 0 .and. scnt_max(p_owd) == 0 ) then
      exit
    endif

  enddo

  call stop_tprof( 'move 2D particles' )
  call write_dbg("part2d_comm", sname, 2, 'ends')

end subroutine move_part2d_comm

! -----------------------------------------------------------------------------
! Copy the particles in the receive buffers belonging to this node to the
! particle array and those not belong to this node to other send buffers.
! This routine only moves particles inside current stage.
! -----------------------------------------------------------------------------
subroutine unpack_relay_particles( part, edge, rbuf, rcnt, src, sbuf, scnt, des )

  implicit none

  class(part2d), intent(inout) :: part
  real, intent(in), dimension(2) :: edge
  real, intent(in), dimension(:) :: rbuf
  real, intent(inout), dimension(:) :: sbuf
  integer, intent(inout) :: scnt, rcnt
  integer, intent(in) :: src, des

  integer :: i, j, npp, part_dim, stay_cnt, go_cnt, stride
  real, dimension(2) :: x

  npp      = part%npp
  part_dim = part%part_dim

  if ( npp + rcnt > part%npmax ) then
    call write_err( 'Particle overflow' )
    ! TODO: resize part2d memory
  endif

  stay_cnt = 0
  go_cnt = 0

  do i = 1, rcnt

    stride = (i - 1) * part_dim

    x(1) = rbuf( 1 + stride )
    x(2) = rbuf( 2 + stride )

    ! check if particle belongs to here
    if ( goto_here( x, edge ) ) then

      stay_cnt = stay_cnt + 1

      part%x( 1, npp + stay_cnt )  = rbuf( 1 + stride )
      part%x( 2, npp + stay_cnt )  = rbuf( 2 + stride )
      part%p( 1, npp + stay_cnt )  = rbuf( 3 + stride )
      part%p( 2, npp + stay_cnt )  = rbuf( 4 + stride )
      part%p( 3, npp + stay_cnt )  = rbuf( 5 + stride )
      part%gamma( npp + stay_cnt ) = rbuf( 6 + stride )
      part%psi( npp + stay_cnt )   = rbuf( 7 + stride )
      part%q( npp + stay_cnt )     = rbuf( 8 + stride )

    ! check if particle goes to the destination neighbor node
    elseif ( goto_des( x, edge, des ) ) then

      scnt = scnt + 1
      go_cnt = go_cnt + 1

      do j = 1, part_dim
        sbuf( j + (scnt - 1) * part_dim ) = rbuf( j + stride )
      enddo

    ! particle goes back to where it comes from, should not happen
    elseif ( goto_des( x, edge, src ) ) then

      call write_err( 'There is returning particle in local particle &
        &communication! ' )

    endif

  enddo

  ! check if all the particles in the receive buffer have been unpacked
  rcnt = rcnt - stay_cnt - go_cnt
  if ( rcnt /= 0 ) then
    call write_err( 'There are unpacked particles!' )
  endif

  part%npp = npp + stay_cnt

end subroutine unpack_relay_particles

subroutine pack_particles( part, edge, sbuf, ihole, scnt, des )

  implicit none

  class(part2d), intent(inout) :: part
  real, intent(in), dimension(2) :: edge
  real, intent(inout), dimension(:) :: sbuf
  integer(kind=LG), intent(inout), dimension(:) :: ihole
  integer, intent(inout) :: scnt
  integer, intent(in) :: des

  integer :: i, npp, go_cnt, part_dim, buf_size, stride
  real, dimension(2) :: x

  npp = part%npp
  part_dim = part%part_dim
  buf_size = size( sbuf ) / part_dim

  go_cnt = 0

  do i = 1, npp

    x = part%x(:,i)

    if ( goto_des( x, edge, des ) ) then

      if ( scnt >= buf_size ) then
        call write_err( '2D particle MPI buffer overflow!' )
      endif

      scnt = scnt + 1
      go_cnt = go_cnt + 1
      ihole(go_cnt) = i

      stride = (scnt - 1) * part_dim
      sbuf( 1 + stride ) = part%x(1,i)
      sbuf( 2 + stride ) = part%x(2,i)
      sbuf( 3 + stride ) = part%p(1,i)
      sbuf( 4 + stride ) = part%p(2,i)
      sbuf( 5 + stride ) = part%p(3,i)
      sbuf( 6 + stride ) = part%gamma(i)
      sbuf( 7 + stride ) = part%psi(i)
      sbuf( 8 + stride ) = part%q(i)

    endif
  enddo

  ! fill the holes inversely
  do i = go_cnt, 1, -1

    part%x( 1:2, ihole(i) ) = part%x( 1:2, npp )
    part%p( 1:3, ihole(i) ) = part%p( 1:3, npp )
    part%gamma( ihole(i) )  = part%gamma(npp)
    part%psi( ihole(i) )    = part%psi(npp)
    part%q( ihole(i) )      = part%q(npp)
    
    npp = npp - 1

  enddo

  part%npp = npp

end subroutine pack_particles

subroutine resize_buf( ratio )
! this subroutine is currently not in use.

  implicit none

  real, intent(in), optional :: ratio

  real, save :: r = 1.5
  real, dimension(:), allocatable :: buf_tmp
  integer :: i, len_old, len_new

  if ( present(ratio) ) r = ratio

  len_old = dim_max * buf_size
  allocate( buf_tmp( len_old ) )
  len_new = dim_max * int( buf_size * r )

  ! resize lower receive buffer
  do i = 1, len_old
    buf_tmp(i) = recv_buf_lower(i)
  enddo
  deallocate( recv_buf_lower )
  allocate( recv_buf_lower( len_new ) )
  recv_buf_lower = 0.0
  do i = 1, len_old
    recv_buf_lower(i) = buf_tmp(i)
  enddo

  ! resize upper receive buffer
  do i = 1, len_old
    buf_tmp(i) = recv_buf_upper(i)
  enddo
  deallocate( recv_buf_upper )
  allocate( recv_buf_upper( len_new ) )
  recv_buf_upper = 0.0
  do i = 1, len_old
    recv_buf_upper(i) = buf_tmp(i)
  enddo

  ! resize lower send buffer
  do i = 1, len_old
    buf_tmp(i) = send_buf_lower(i)
  enddo
  deallocate( send_buf_lower )
  allocate( send_buf_lower( len_new ) )
  send_buf_lower = 0.0
  do i = 1, len_old
    send_buf_lower(i) = buf_tmp(i)
  enddo

  ! resize upper send buffer
  do i = 1, len_old
    buf_tmp(i) = send_buf_upper(i)
  enddo
  deallocate( send_buf_upper )
  allocate( send_buf_upper( len_new ) )
  send_buf_upper = 0.0
  do i = 1, len_old
    send_buf_upper(i) = buf_tmp(i)
  enddo

  ! resize index of holes
  deallocate( ihole )
  allocate( ihole( len_new / dim_max ) )

  deallocate( buf_tmp )

end subroutine resize_buf

function get_phys_bnd( proc_id, nvp ) result( res )

  implicit none

  integer, intent(in) :: proc_id, nvp
  logical, dimension(2) :: res

  integer :: loc_id

  res = .false.

  loc_id = mod( proc_id, nvp )

  if ( loc_id == 0 )     res( p_iwd ) = .true.
  if ( loc_id == nvp-1 ) res( p_owd ) = .true.

end function get_phys_bnd

function goto_here( x, edge )

  implicit none

  real, intent(in), dimension(2) :: x
  real, intent(in), dimension(2) :: edge
  logical :: goto_here

  real :: pos

  pos = sqrt( x(1)*x(1) + x(2)*x(2) )
  goto_here = pos >= edge(p_lower) .and. pos < edge(p_upper)

end function goto_here

function goto_des( x, edge, des )

  implicit none

  real, intent(in), dimension(2) :: x
  real, intent(in), dimension(2) :: edge
  integer, intent(in) :: des
  logical :: goto_des

  real :: pos
  
  goto_des = .false.

  select case ( des )
  case ( p_iwd )
    pos = sqrt( x(1)*x(1) + x(2)*x(2) )
    goto_des = pos < edge(p_lower)
  case ( p_owd )
    pos = sqrt( x(1)*x(1) + x(2)*x(2) )
    goto_des = pos >= edge(p_upper)
  end select

end function goto_des

end module part2d_comm
