module part2d_comm

use sysutil_module
use param
use mpi
use parallel_module
use options_class
use ufield_class
use part2d_class

implicit none

private

save

! direction indicators
! outward/inward: positive/negative direction in r
integer, parameter :: p_self = 0, &
                      p_iwd = 1, &
                      p_owd = 2

integer, parameter :: iter_max = 1000

! send and receive buffer for MPI
real, dimension(:), allocatable :: send_buf_lower
real, dimension(:), allocatable :: send_buf_upper
real, dimension(:), allocatable :: recv_buf_lower
real, dimension(:), allocatable :: recv_buf_upper

! index of holes
integer(kind=LG), dimension(:), allocatable :: ihole

! temporary arrays used to reallocate buffer
real, dimension(:), allocatable :: tmp_real
integer(kind=LG), dimension(:), allocatable :: tmp_int

! count of MPI buffers
integer, dimension(2) :: send_cnt = 0
integer, dimension(2) :: recv_cnt = 0

! max dimension of particle coordinates
integer :: dim_max = 9

! current buffer size, same for send and receive buffers
integer :: buf_size = 0

! increment for buffer reallocation
real, parameter :: p_buf_incr = 1.5

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

  if ( allocated( send_buf_lower ) ) deallocate( send_buf_lower )
  if ( allocated( send_buf_upper ) ) deallocate( send_buf_upper )
  if ( allocated( recv_buf_lower ) ) deallocate( recv_buf_lower )
  if ( allocated( recv_buf_upper ) ) deallocate( recv_buf_upper )
  if ( allocated( ihole ) ) deallocate( ihole )
  if ( allocated( tmp_real ) ) deallocate( tmp_real )
  if ( allocated( tmp_int ) ) deallocate( tmp_int )

end subroutine end_part2d_comm

subroutine move_part2d_comm( part )

  implicit none

  class( part2d ), intent(inout) :: part

  integer :: rsize, ssize, ierr, part_dim, iter, loc_grp, world
  integer, dimension(MPI_STATUS_SIZE) :: istat
  integer, dimension(2) :: sid, rid
  integer, dimension(2) :: scnt_max

  character(len=18), save :: sname = 'move_part2d_comm'

  world    = comm_world()
  loc_grp  = comm_loc()
  part_dim = part%part_dim

  call write_dbg("part2d_comm", sname, 2, 'starts')
  call start_tprof( 'move 2D particles' )

  ! pack particles to inward/outward send buffers
  call pack_particles( part )

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

    ! receive from outward and unpack particles
    if ( .not. phys_bnd(p_owd) ) then

      ! get the message size and resize receiving buffer if necessary
      call mpi_probe( pid(p_owd), iter, world, istat, ierr )
      call mpi_get_count( istat, p_dtype_real, recv_cnt(p_owd), ierr )
      recv_cnt(p_owd) = recv_cnt(p_owd) / part_dim

      rsize = size( recv_buf_upper )
      if ( recv_cnt(p_owd) > rsize / part_dim ) then
        call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
          &particle upper MPI receiving buffer!', only_root = .false. )
        rsize = int( recv_cnt(p_owd) * p_buf_incr ) * part_dim
        deallocate( recv_buf_upper )
        allocate( recv_buf_upper( rsize ) )
      endif

      call mpi_irecv( recv_buf_upper, rsize, p_dtype_real, pid(p_owd), iter, &
        world, rid(2), ierr )
    else
      rid(2) = MPI_REQUEST_NULL
    endif

    ! receive from inward and unpack particles
    if ( .not. phys_bnd(p_iwd) ) then

      ! get the message size and resize receiving buffer if necessary
      call mpi_probe( pid(p_iwd), iter, world, istat, ierr )
      call mpi_get_count( istat, p_dtype_real, recv_cnt(p_iwd), ierr )
      recv_cnt(p_iwd) = recv_cnt(p_iwd) / part_dim

      rsize = size( recv_buf_lower )
      if ( recv_cnt(p_iwd) > rsize / part_dim ) then
        call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
          &particle lower MPI receiving buffer!', only_root = .false. )
        rsize = int( recv_cnt(p_iwd) * p_buf_incr ) * part_dim
        deallocate( recv_buf_lower )
        allocate( recv_buf_lower( rsize ) )
      endif

      call mpi_irecv( recv_buf_lower, rsize, p_dtype_real, pid(p_iwd), &
        iter, world, rid(1), ierr )
    else
      rid(1) = MPI_REQUEST_NULL
    endif

    ! wait sending finish
    call mpi_wait( sid(1), MPI_STATUS_IGNORE, ierr ); send_cnt(p_iwd) = 0
    call mpi_wait( sid(2), MPI_STATUS_IGNORE, ierr ); send_cnt(p_owd) = 0

    ! wait receiving finish
    call mpi_wait( rid(2), MPI_STATUS_IGNORE, ierr )
    call mpi_wait( rid(1), MPI_STATUS_IGNORE, ierr )

    ! unpack particles in inward/outward receive buffer
    call unpack_relay_particles( part )

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
subroutine unpack_relay_particles( part )

  implicit none

  class(part2d), intent(inout) :: part

  integer :: i, j, npp, part_dim, stay_cnt, go_cnt, stride
  integer, dimension(2) :: sbuf_cnt
  real, dimension(2) :: x
  real :: ratio

  npp      = part%npp
  part_dim = part%part_dim
  sbuf_cnt(p_iwd) = size( send_buf_lower ) / part_dim
  sbuf_cnt(p_owd) = size( send_buf_upper ) / part_dim

  ! if all the received particles stay in this partition, check particle buffer size
  ratio = real( npp + sum(recv_cnt) ) / part%npmax
  if ( ratio > 1.0 ) then
    call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
      &particle buffer!', only_root = .false. )
    call part%realloc( ratio = ratio * p_buf_incr )
  endif

  ! if all the particles in the lower receiving buffer go to the outward partition,
  ! check the upper sending buffer size
  if ( recv_cnt(p_iwd) + send_cnt(p_owd) > sbuf_cnt(p_owd) ) then
    call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
      &particle sending buffer!', only_root = .false. )
    sbuf_cnt(p_owd) = ( recv_cnt(p_iwd) + send_cnt(p_owd) ) * p_buf_incr
    deallocate( send_buf_upper )
    allocate( send_buf_upper( sbuf_cnt(p_owd) * part_dim ) )
  endif

  ! if all the particles in the upper receiving buffer go to the inward partition,
  ! check the lower sending buffer size
  if ( recv_cnt(p_owd) + send_cnt(p_iwd) > sbuf_cnt(p_iwd) ) then
    call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
      &particle sending buffer!', only_root = .false. )
    sbuf_cnt(p_iwd) = ( recv_cnt(p_owd) + send_cnt(p_iwd) ) * p_buf_incr
    deallocate( send_buf_lower )
    allocate( send_buf_lower( sbuf_cnt(p_iwd) * part_dim ) )
  endif

  stay_cnt = 0; go_cnt = 0

  ! unpack and relay the particle in the lower receiving buffer
  do i = 1, recv_cnt(p_iwd)

    stride = (i - 1) * part_dim

    x(1) = recv_buf_lower( 1 + stride )
    x(2) = recv_buf_lower( 2 + stride )

    ! check if particle belongs to here
    if ( goto_here( x, redge ) ) then

      stay_cnt = stay_cnt + 1

      part%x_l( 1, npp + stay_cnt )  = recv_buf_lower( 1 + stride )
      part%x_l( 2, npp + stay_cnt )  = recv_buf_lower( 2 + stride )
      part%x( 1, npp + stay_cnt )  = recv_buf_lower( 3 + stride )
      part%x( 2, npp + stay_cnt )  = recv_buf_lower( 4 + stride )
!       part%x_r( 1, npp + stay_cnt )  = recv_buf_lower( 5 + stride )
!       part%x_r( 2, npp + stay_cnt )  = recv_buf_lower( 6 + stride )
      part%p_l( 1, npp + stay_cnt )  = recv_buf_lower( 5 + stride )
      part%p_l( 2, npp + stay_cnt )  = recv_buf_lower( 6 + stride )
      part%p_l( 3, npp + stay_cnt )  = recv_buf_lower( 7 + stride )
      part%p( 1, npp + stay_cnt )  = recv_buf_lower( 8 + stride )
      part%p( 2, npp + stay_cnt )  = recv_buf_lower( 9 + stride )
      part%p( 3, npp + stay_cnt )  = recv_buf_lower( 10 + stride )
      part%gamma( npp + stay_cnt ) = recv_buf_lower( 11 + stride )
      part%psi( npp + stay_cnt )   = recv_buf_lower( 12 + stride )
      part%q( npp + stay_cnt )     = recv_buf_lower( 13 + stride )
      part%w( npp + stay_cnt )     = recv_buf_lower( 14 + stride )
      part%w0( npp + stay_cnt )     = recv_buf_lower( 15 + stride )

    ! check if particle goes to the outward neighbor node
    elseif ( goto_des( x, redge, p_owd ) ) then

      do j = 1, part_dim
        send_buf_upper( j + send_cnt(p_owd) * part_dim ) = recv_buf_lower( j + stride )
      enddo
      send_cnt(p_owd) = send_cnt(p_owd) + 1
      go_cnt = go_cnt + 1

    ! particle goes back to where it comes from, should not happen
    elseif ( goto_des( x, redge, p_iwd ) ) then

      call write_err( 'There is returning particle in local particle &
        &communication! ' )

    endif

  enddo

  ! check if the particles in the receive buffer have been unpacked
  if ( stay_cnt + go_cnt /= recv_cnt(p_iwd) ) then
    call write_err( 'There are unpacked particles!' )
  else
    recv_cnt(p_iwd) = 0
    npp = npp + stay_cnt
  endif

  stay_cnt = 0; go_cnt = 0

  ! unpack and relay the particle in the upper receiving buffer
  do i = 1, recv_cnt(p_owd)

    stride = (i - 1) * part_dim

    x(1) = recv_buf_upper( 1 + stride )
    x(2) = recv_buf_upper( 2 + stride )

    ! check if particle belongs to here
    if ( goto_here( x, redge ) ) then

      stay_cnt = stay_cnt + 1

      part%x_l( 1, npp + stay_cnt )  = recv_buf_upper( 1 + stride )
      part%x_l( 2, npp + stay_cnt )  = recv_buf_upper( 2 + stride )
      part%x( 1, npp + stay_cnt )  = recv_buf_upper( 3 + stride )
      part%x( 2, npp + stay_cnt )  = recv_buf_upper( 4 + stride )
!       part%x_r( 1, npp + stay_cnt )  = recv_buf_upper( 5 + stride )
!       part%x_r( 2, npp + stay_cnt )  = recv_buf_upper( 6 + stride )
      part%p_l( 1, npp + stay_cnt )  = recv_buf_upper( 5 + stride )
      part%p_l( 2, npp + stay_cnt )  = recv_buf_upper( 6 + stride )
      part%p_l( 3, npp + stay_cnt )  = recv_buf_upper( 7 + stride )
      part%p( 1, npp + stay_cnt )  = recv_buf_upper( 8 + stride )
      part%p( 2, npp + stay_cnt )  = recv_buf_upper( 9 + stride )
      part%p( 3, npp + stay_cnt )  = recv_buf_upper( 10 + stride )
      part%gamma( npp + stay_cnt ) = recv_buf_upper( 11 + stride )
      part%psi( npp + stay_cnt )   = recv_buf_upper( 12 + stride )
      part%q( npp + stay_cnt )     = recv_buf_upper( 13 + stride )
      part%w( npp + stay_cnt )     = recv_buf_upper( 14 + stride )
      part%w0( npp + stay_cnt )    = recv_buf_upper( 15 + stride )

    ! check if particle goes to the inward neighbor node
    elseif ( goto_des( x, redge, p_iwd ) ) then

      do j = 1, part_dim
        send_buf_lower( j + send_cnt(p_iwd) * part_dim ) = recv_buf_upper( j + stride )
      enddo
      send_cnt(p_iwd) = send_cnt(p_iwd) + 1
      go_cnt = go_cnt + 1

    ! particle goes back to where it comes from, should not happen
    elseif ( goto_des( x, redge, p_owd ) ) then

      call write_err( 'There is returning particle in local particle &
        &communication! ' )

    endif

  enddo

  ! check if the particles in the receive buffer have been unpacked
  if ( stay_cnt + go_cnt /= recv_cnt(p_owd) ) then
    call write_err( 'There are unpacked particles!' )
  else
    recv_cnt(p_owd) = 0
    npp = npp + stay_cnt
  endif

  part%npp = npp

end subroutine unpack_relay_particles

subroutine pack_particles( part )

  implicit none

  class(part2d), intent(inout) :: part

  integer :: i, npp, go_cnt, ihole_cnt, stride, part_dim
  integer, dimension(2) :: sbuf_cnt
  real, dimension(2) :: x

  npp             = part%npp
  part_dim        = part%part_dim
  sbuf_cnt(p_iwd) = size( send_buf_lower ) / part%part_dim
  sbuf_cnt(p_owd) = size( send_buf_upper ) / part%part_dim
  ihole_cnt       = size( ihole )

  go_cnt = 0

  do i = 1, npp

    x = part%x(:,i)

    ! particles go to the inward partition
    if ( goto_des( x, redge, p_iwd ) ) then

      stride = send_cnt(p_iwd) * part_dim
      send_buf_lower( 1 + stride ) = part%x_l(1,i)
      send_buf_lower( 2 + stride ) = part%x_l(2,i)
      send_buf_lower( 3 + stride ) = part%x(1,i)
      send_buf_lower( 4 + stride ) = part%x(2,i)
!       send_buf_lower( 5 + stride ) = part%x_r(1,i)
!       send_buf_lower( 6 + stride ) = part%x_r(2,i)
      send_buf_lower( 5 + stride ) = part%p_l(1,i)
      send_buf_lower( 6 + stride ) = part%p_l(2,i)
      send_buf_lower( 7 + stride ) = part%p_l(3,i)
      send_buf_lower( 8 + stride ) = part%p(1,i)
      send_buf_lower( 9 + stride ) = part%p(2,i)
      send_buf_lower( 10 + stride ) = part%p(3,i)
      send_buf_lower( 11 + stride ) = part%gamma(i)
      send_buf_lower( 12 + stride ) = part%psi(i)
      send_buf_lower( 13 + stride ) = part%q(i)
      send_buf_lower( 14 + stride ) = part%w(i)
      send_buf_lower( 15 + stride ) = part%w0(i)

      send_cnt(p_iwd) = send_cnt(p_iwd) + 1
      go_cnt = go_cnt + 1
      ihole(go_cnt) = i

    ! particles go to the outward partition
    elseif ( goto_des( x, redge, p_owd ) ) then

      stride = send_cnt(p_owd) * part_dim
      send_buf_upper( 1 + stride ) = part%x_l(1,i)
      send_buf_upper( 2 + stride ) = part%x_l(2,i)
      send_buf_upper( 3 + stride ) = part%x(1,i)
      send_buf_upper( 4 + stride ) = part%x(2,i)
!       send_buf_upper( 5 + stride ) = part%x_r(1,i)
!       send_buf_upper( 6 + stride ) = part%x_r(2,i)
      send_buf_upper( 5 + stride ) = part%p_l(1,i)
      send_buf_upper( 6 + stride ) = part%p_l(2,i)
      send_buf_upper( 7 + stride ) = part%p_l(3,i)
      send_buf_upper( 8 + stride ) = part%p(1,i)
      send_buf_upper( 9 + stride ) = part%p(2,i)
      send_buf_upper( 10 + stride ) = part%p(3,i)
      send_buf_upper( 11 + stride ) = part%gamma(i)
      send_buf_upper( 12 + stride ) = part%psi(i)
      send_buf_upper( 13 + stride ) = part%q(i)
      send_buf_upper( 14 + stride ) = part%w(i)
      send_buf_upper( 15 + stride ) = part%w0(i)

      send_cnt(p_owd) = send_cnt(p_owd) + 1
      go_cnt = go_cnt + 1
      ihole(go_cnt) = i

    endif

    ! check lower sending buffer size
    if ( send_cnt(p_iwd) >= sbuf_cnt(p_iwd) ) then
      call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
        &particle MPI sending buffer!', only_root = .false. )
      sbuf_cnt(p_iwd) = send_cnt(p_iwd) * p_buf_incr
      allocate( tmp_real( sbuf_cnt(p_iwd) * part_dim ) )
      tmp_real = 0.0
      tmp_real( 1:size(send_buf_lower) ) = send_buf_lower
      call move_alloc( tmp_real, send_buf_lower )
    endif

    ! check upper sending buffer size
    if ( send_cnt(p_owd) >= sbuf_cnt(p_owd) ) then
      call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
        &particle MPI sending buffer!', only_root = .false. )
      sbuf_cnt(p_owd) = send_cnt(p_owd) * p_buf_incr
      allocate( tmp_real( sbuf_cnt(p_owd) * part_dim ) )
      tmp_real = 0.0
      tmp_real( 1:size(send_buf_upper) ) = send_buf_upper
      call move_alloc( tmp_real, send_buf_upper )
    endif

    ! check ihole array size
    if ( go_cnt >= ihole_cnt ) then
      call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
        &particle MPI buffer!', only_root = .false. )
      ihole_cnt = go_cnt * p_buf_incr
      allocate( tmp_int( ihole_cnt ) )
      tmp_int = 0
      tmp_int( 1:size(ihole) ) = ihole
      call move_alloc( tmp_int, ihole )
    endif

  enddo

  ! fill the holes inversely
  do i = go_cnt, 1, -1

    part%x_l( 1:2, ihole(i) ) = part%x_l( 1:2, npp )
    part%x( 1:2, ihole(i) ) = part%x( 1:2, npp )
!     part%x_r( 1:2, ihole(i) ) = part%x_r( 1:2, npp )
    part%p_l( 1:3, ihole(i) ) = part%p_l( 1:3, npp )
    part%p( 1:3, ihole(i) ) = part%p( 1:3, npp )
    part%gamma( ihole(i) )  = part%gamma(npp)
    part%psi( ihole(i) )    = part%psi(npp)
    part%q( ihole(i) )      = part%q(npp)
    part%w( ihole(i) )      = part%w(npp)
    part%w0( ihole(i) )     = part%w0(npp)


    npp = npp - 1

  enddo

  part%npp = npp

end subroutine pack_particles

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
