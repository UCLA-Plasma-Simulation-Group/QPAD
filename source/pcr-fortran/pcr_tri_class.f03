module pcr_tri_class

use mpi

implicit none
private

type :: pcr_tri

  private

  ! MPI information
  integer :: mpi_comm, myid, num_procs, mpi_dtype

  ! MPI send/recv buffer
  real, dimension(:), pointer :: sbuf, rbuf

  ! Global problem size and number of rows in this partition
  integer :: n_global, n_local

  ! a, b, c are the lower, central and upper diagonal elements
  ! rhs is the right-hand-side of the linear system
  real, dimension(:), allocatable :: a, b, c, rhs, x

  ! total steps of cyclic reduction
  integer :: tot_steps

  ! steps of cyclic redution that only needs neighbor processors' data
  integer :: neighbor_steps

  logical :: share_mpi_buf = .false.

  contains

  procedure :: create            => create_pcr_tri
  procedure :: destroy           => destroy_pcr_tri
  procedure :: solve             => solve_pcr_tri
  generic   :: set_values_rhs    => set_values_rhs_pcr_tri, set_values_rhs_byrow_pcr_tri
  generic   :: set_values_matrix => set_values_matrix_pcr_tri, set_values_matrix_byrow_pcr_tri
  generic   :: get_values_matrix => get_values_matrix_pcr_tri, get_values_matrix_byrow_pcr_tri
  generic   :: get_values_x      => get_values_x_pcr_tri, get_values_x_byrow_pcr_tri
  procedure :: print_x           => print_x_pcr_tri
  procedure :: print_rhs         => print_rhs_pcr_tri
  procedure, private :: pack_data      => pack_data_pcr_tri
  procedure, private :: unpack_data    => unpack_data_pcr_tri
  procedure, private :: set_idle_ghost => set_idle_ghost_pcr_tri
  procedure, private :: set_values_rhs_pcr_tri
  procedure, private :: set_values_rhs_byrow_pcr_tri
  procedure, private :: set_values_matrix_pcr_tri
  procedure, private :: get_values_matrix_pcr_tri
  procedure, private :: set_values_matrix_byrow_pcr_tri
  procedure, private :: get_values_matrix_byrow_pcr_tri
  procedure, private :: get_values_x_pcr_tri
  procedure, private :: get_values_x_byrow_pcr_tri

end type pcr_tri

public :: pcr_tri

integer, parameter :: stdout = 6, stderr = 0
integer, parameter :: p_lower = 0, p_upper = 1

! shared MPI buffers
real, dimension(:), pointer :: mpi_buf => null()
integer, save :: n_shared_solver = 0

contains

subroutine create_pcr_tri( this, n_global, mpi_comm, buf_mode, err )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  integer, intent(in) :: n_global, mpi_comm
  character(len=*), intent(in) :: buf_mode
  integer, intent(out) :: err

  integer :: ierr

  this%n_global = n_global
  this%mpi_comm = mpi_comm

  call MPI_COMM_RANK( this%mpi_comm, this%myid, ierr )
  call MPI_COMM_SIZE( this%mpi_comm, this%num_procs, ierr )

  if ( .not. is_power_of_2(n_global) ) then
    write(stderr, *) "The problem size is not the power of 2."
    err = 1
    return
  endif
  this%tot_steps = nint( log2( real(this%n_global) ) )

  if ( .not. is_power_of_2(this%num_procs) ) then
    write(stderr, *) "The number of processors is not the power of 2."
    err = 1
    return
  endif

  this%n_local = this%n_global / this%num_procs
  this%neighbor_steps = nint( log2( real(this%n_local) ) ) + 1
  allocate( this%a(1-this%n_local:2*this%n_local) ); this%a = 0.0
  allocate( this%b(1-this%n_local:2*this%n_local) ); this%b = 0.0
  allocate( this%c(1-this%n_local:2*this%n_local) ); this%c = 0.0
  allocate( this%x(this%n_local) ); this%x = 0.0
  allocate( this%rhs(1-this%n_local:2*this%n_local) ); this%rhs = 0.0

  allocate( this%sbuf(4*this%n_local), this%rbuf(4*this%n_local) )

  if ( trim(buf_mode) == 'shared' ) then

    if ( .not. associated(mpi_buf) ) then
      allocate( mpi_buf( 8 * this%n_local ) )
    else if ( associated(mpi_buf) .and. size(mpi_buf) < 8 * this%n_local ) then
      deallocate(mpi_buf)
      allocate( mpi_buf( 8 * this%n_local ) )
    endif
    this%sbuf => mpi_buf(                1 : 4*this%n_local )
    this%rbuf => mpi_buf( 4*this%n_local+1 : 8*this%n_local )
    this%share_mpi_buf = .true.
    n_shared_solver = n_shared_solver + 1

  else if ( trim(buf_mode) == 'exclusive' ) then
    this%share_mpi_buf = .false.
    allocate( this%sbuf(4*this%n_local), this%rbuf(4*this%n_local) )
  else
    write(stderr, *) "Invalid MPI buffer type!"
    err = 1
    return
  endif

  if ( digits(1.0) > 24 ) then
    this%mpi_dtype = MPI_DOUBLE_PRECISION
  else
    this%mpi_dtype = MPI_REAL
  endif

  err = 0

end subroutine create_pcr_tri

subroutine destroy_pcr_tri( this )

  implicit none
  class( pcr_tri ), intent(inout) :: this

  deallocate( this%a, this%b, this%c, this%x, this%rhs )
  if ( this%share_mpi_buf ) then
    n_shared_solver = n_shared_solver - 1
    if ( n_shared_solver == 0 ) deallocate( mpi_buf )
  else
    deallocate( this%sbuf, this%rbuf )
  endif

end subroutine destroy_pcr_tri

subroutine solve_pcr_tri( this )

  implicit none
  class( pcr_tri ), intent(inout) :: this

  integer :: step, stride, pstride, i, pid_lower, pid_upper, comm_size
  integer :: sid, rid, ierr
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:), pointer :: a_tmp => null(), b_tmp => null(), c_tmp => null(), rhs_tmp => null()
  real :: coefm, coefp
  integer :: im, ip

  stride = 1; pstride = 1
  do step = 1, this%tot_steps

    if ( step <= this%neighbor_steps ) then
      pstride =  1
      comm_size = stride
    else
      pstride = pstride * 2
      comm_size = this%n_local
    endif

    pid_lower = this%myid - pstride
    pid_upper = this%myid + pstride

    ! send data to lower neighbor processor
    if ( pid_lower >= 0 ) then
      call this%pack_data( 1, comm_size )
      call MPI_ISEND( this%sbuf, 4 * comm_size, this%mpi_dtype, pid_lower, 1, this%mpi_comm, sid, ierr )
    endif
    if ( pid_upper <= this%num_procs - 1 ) then
      call MPI_RECV( this%rbuf, 4 * comm_size, this%mpi_dtype, pid_upper, 1, this%mpi_comm, stat, ierr )
      call this%unpack_data( this%n_local + 1, this%n_local + comm_size )
    else
      call this%set_idle_ghost( p_upper )
    endif
    if ( pid_lower >= 0 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! send data to upper neighbor processor
    if ( pid_upper <= this%num_procs - 1 ) then
      call this%pack_data( this%n_local - comm_size + 1, this%n_local )
      call MPI_ISEND( this%sbuf, 4 * comm_size, this%mpi_dtype, pid_upper, 1, this%mpi_comm, sid, ierr )
    endif
    if ( pid_lower >= 0 ) then
      call MPI_RECV( this%rbuf, 4 * comm_size, this%mpi_dtype, pid_lower, 1, this%mpi_comm, stat, ierr )
      call this%unpack_data( 1 - comm_size, 0 )
    else
      call this%set_idle_ghost( p_lower )
    endif
    if ( pid_upper <= this%num_procs - 1 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! use sbuf as the temporary arrays
    this%sbuf = 0.0
    a_tmp   => this%sbuf(                1 :   this%n_local )
    b_tmp   => this%sbuf(   this%n_local+1 : 2*this%n_local )
    c_tmp   => this%sbuf( 2*this%n_local+1 : 3*this%n_local )
    rhs_tmp => this%sbuf( 3*this%n_local+1 : 4*this%n_local )

    do i = 1, this%n_local
      im = i - comm_size
      ip = i + comm_size
      coefm = -this%a(i) / this%b(im)
      coefp = -this%c(i) / this%b(ip)
      a_tmp(i) = coefm * this%a(im)
      b_tmp(i) = coefm * this%c(im) + this%b(i) + coefp * this%a(ip)
      c_tmp(i) = coefp * this%c(ip)
      rhs_tmp(i) = this%rhs(i) + coefm * this%rhs(im) + coefp * this%rhs(ip)
    enddo

    this%a(1:this%n_local) = a_tmp(1:this%n_local)
    this%b(1:this%n_local) = b_tmp(1:this%n_local)
    this%c(1:this%n_local) = c_tmp(1:this%n_local)
    this%rhs(1:this%n_local) = rhs_tmp(1:this%n_local)

    stride = stride * 2

  enddo

  do i = 1, this%n_local
    this%x(i) = this%rhs(i) / this%b(i)
  enddo

end subroutine solve_pcr_tri

subroutine pack_data_pcr_tri( this, i0, i1 )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  integer, intent(in) :: i0, i1

  integer :: i, num
  num = i1 - i0 + 1
  this%sbuf(      1:  num) = this%a(i0:i1)
  this%sbuf(  num+1:2*num) = this%b(i0:i1)
  this%sbuf(2*num+1:3*num) = this%c(i0:i1)
  this%sbuf(3*num+1:4*num) = this%rhs(i0:i1)

end subroutine pack_data_pcr_tri

subroutine unpack_data_pcr_tri( this, i0, i1 )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  integer, intent(in) :: i0, i1

  integer :: i, num
  num = i1 - i0 + 1
  this%a(i0:i1)   = this%rbuf(      1:  num)
  this%b(i0:i1)   = this%rbuf(  num+1:2*num)
  this%c(i0:i1)   = this%rbuf(2*num+1:3*num)
  this%rhs(i0:i1) = this%rbuf(3*num+1:4*num)

end subroutine unpack_data_pcr_tri

subroutine set_idle_ghost_pcr_tri( this, bnd )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  integer, intent(in) :: bnd

  integer :: i0, i1

  if ( bnd == p_lower ) then
    i0 = 1 - this%n_local
    i1 = 0
    this%a(i0:i1) = 0.0
    this%b(i0:i1) = 1.0
    this%c(i0:i1) = 0.0
    this%rhs(i0:i1) = 0.0
  endif

  if ( bnd == p_upper ) then
    i0 = this%n_local + 1
    i1 = 2 * this%n_local
    this%a(i0:i1) = 0.0
    this%b(i0:i1) = 1.0
    this%c(i0:i1) = 0.0
    this%rhs(i0:i1) = 0.0
  endif

end subroutine set_idle_ghost_pcr_tri

subroutine set_values_rhs_pcr_tri( this, rhs )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  real, intent(in), dimension(:) :: rhs
  integer :: i

  do i = 1, this%n_local
    this%rhs(i) = rhs(i)
  enddo

end subroutine set_values_rhs_pcr_tri

subroutine set_values_rhs_byrow_pcr_tri( this, rhs, i_loc )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  real, intent(in) :: rhs
  integer, intent(in) :: i_loc
  this%rhs(i_loc) = rhs

end subroutine set_values_rhs_byrow_pcr_tri

subroutine set_values_matrix_pcr_tri( this, a, b, c )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  real, intent(in), dimension(:) :: a, b, c
  integer :: i

  do i = 1, this%n_local
    this%a(i) = a(i)
    this%b(i) = b(i)
    this%c(i) = c(i)
  enddo
  if ( this%myid == 0 ) this%a(1) = 0.0
  if ( this%myid == this%num_procs - 1 ) this%c(this%n_local) = 0.0

end subroutine set_values_matrix_pcr_tri

subroutine set_values_matrix_byrow_pcr_tri( this, a, b, c, i_loc )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  real, intent(in) :: a, b, c
  integer, intent(in) :: i_loc
  integer :: i

  this%a(i_loc) = a
  this%b(i_loc) = b
  this%c(i_loc) = c
  if ( this%myid == 0 .and. i_loc == 1 ) this%a(i_loc) = 0.0
  if ( this%myid == this%num_procs - 1 .and. i_loc == this%n_local ) this%c(i_loc) = 0.0

end subroutine set_values_matrix_byrow_pcr_tri

subroutine get_values_matrix_pcr_tri( this, a, b, c )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  real, intent(out), dimension(:) :: a, b, c
  integer :: i

  do i = 1, this%n_local
    a(i) = this%a(i)
    b(i) = this%b(i)
    c(i) = this%c(i)
  enddo

end subroutine get_values_matrix_pcr_tri

subroutine get_values_matrix_byrow_pcr_tri( this, a, b, c, i_loc )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  real, intent(out) :: a, b, c
  integer, intent(in) :: i_loc
  integer :: i

  a = this%a(i_loc)
  b = this%b(i_loc)
  c = this%c(i_loc)

end subroutine get_values_matrix_byrow_pcr_tri

subroutine get_values_x_pcr_tri( this, x )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  real, intent(out), dimension(:) :: x
  integer :: i

  do i = 1, this%n_local
    x(i) = this%x(i)
  enddo

end subroutine get_values_x_pcr_tri

subroutine get_values_x_byrow_pcr_tri( this, x, i_loc )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  real, intent(out) :: x
  integer, intent(in) :: i_loc
  x = this%x(i_loc)

end subroutine get_values_x_byrow_pcr_tri

subroutine print_x_pcr_tri( this, filename )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  character(len=*), intent(in) :: filename
  integer :: unit, i
  character(len=128) :: fstr

  write( fstr, '(I0.4)' ) this%myid
  fstr = trim(filename) // '.' // trim(fstr)
  unit = 2
  open( unit, file=trim(fstr) )
  do i = 1, this%n_local
    write( unit, '(ES30.15)' ) this%x(i)
  enddo

end subroutine print_x_pcr_tri

subroutine print_rhs_pcr_tri( this, filename )

  implicit none
  class( pcr_tri ), intent(inout) :: this
  character(len=*), intent(in) :: filename
  integer :: unit, i
  character(len=128) :: fstr

  write( fstr, '(I0.4)' ) this%myid
  fstr = trim(filename) // '.' // trim(fstr)
  unit = 2
  open( unit, file=trim(fstr) )
  do i = 1, this%n_local
    write( unit, '(ES30.15)' ) this%rhs(i)
  enddo

end subroutine print_rhs_pcr_tri

function log2(x)
  implicit none
  real, intent(in) :: x
  real :: log2
  log2 = log(x) / log(2.0)
end function log2

function is_power_of_2(x)
  implicit none
  integer, intent(in) :: x
  logical :: is_power_of_2
  is_power_of_2 = iand(x, x-1) == 0
end function is_power_of_2

end module pcr_tri_class