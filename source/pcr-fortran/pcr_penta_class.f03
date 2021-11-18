module pcr_penta_class

use mpi

implicit none
private

type :: pcr_penta

  private

  ! MPI information
  integer :: mpi_comm, myid, num_procs, mpi_dtype

  ! MPI send/recv buffer
  real, dimension(:), pointer :: sbuf, rbuf

  ! Global problem size and number of rows in this partition
  integer :: n_global, n_local

  ! a, b, c are the lower, central and upper diagonal elements
  ! rhs is the right-hand-side of the linear system
  real, dimension(:), allocatable :: a, b, c, d, e, rhs, x

  ! total steps of cyclic reduction
  integer :: tot_steps

  ! steps of cyclic redution that only needs neighbor processors' data
  integer :: neighbor_steps

  logical :: share_mpi_buf = .false.

  contains

  procedure :: create            => create_pcr_penta
  procedure :: destroy           => destroy_pcr_penta
  procedure :: solve             => solve_pcr_penta
  procedure :: set_values_rhs    => set_values_rhs_pcr_penta
  procedure :: set_values_matrix => set_values_matrix_pcr_penta
  procedure :: get_values_x      => get_values_x_pcr_penta
  procedure :: print_x           => print_x_pcr_penta
  procedure :: print_rhs         => print_rhs_pcr_penta
  procedure, private :: pack_data      => pack_data_pcr_penta
  procedure, private :: unpack_data    => unpack_data_pcr_penta
  procedure, private :: set_idle_ghost => set_idle_ghost_pcr_penta

end type pcr_penta

public :: pcr_penta

integer, parameter :: stdout = 6, stderr = 0
integer, parameter :: p_lower = 0, p_upper = 1
integer, parameter :: p_ghost_inner = 0, p_ghost_outter = 1

! shared MPI buffers
real, dimension(:), pointer :: mpi_buf => null()
integer, save :: n_shared_solver = 0

contains

subroutine create_pcr_penta( this, n_global, mpi_comm, buf_mode, err )

  implicit none
  class( pcr_penta ), intent(inout) :: this
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
  allocate( this%a(1-2*this%n_local:3*this%n_local) ); this%a = 0.0
  allocate( this%b(1-2*this%n_local:3*this%n_local) ); this%b = 0.0
  allocate( this%c(1-2*this%n_local:3*this%n_local) ); this%c = 0.0
  allocate( this%d(1-2*this%n_local:3*this%n_local) ); this%d = 0.0
  allocate( this%e(1-2*this%n_local:3*this%n_local) ); this%e = 0.0
  allocate( this%x(this%n_local) ); this%x = 0.0
  allocate( this%rhs(1-2*this%n_local:3*this%n_local) ); this%rhs = 0.0

  if ( trim(buf_mode) == 'shared' ) then

    if ( .not. associated(mpi_buf) ) then
      allocate( mpi_buf( 12 * this%n_local ) )
    else if ( associated(mpi_buf) .and. size(mpi_buf) < 12 * this%n_local ) then
      deallocate(mpi_buf)
      allocate( mpi_buf( 12 * this%n_local ) )
    endif
    this%sbuf => mpi_buf(                1 :  6*this%n_local )
    this%rbuf => mpi_buf( 6*this%n_local+1 : 12*this%n_local )
    this%share_mpi_buf = .true.
    n_shared_solver = n_shared_solver + 1

  else if ( trim(buf_mode) == 'exclusive' ) then
    this%share_mpi_buf = .false.
    allocate( this%sbuf(6*this%n_local), this%rbuf(6*this%n_local) )
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

end subroutine create_pcr_penta

subroutine destroy_pcr_penta( this )

  implicit none
  class( pcr_penta ), intent(inout) :: this

  deallocate( this%a, this%b, this%c, this%d, this%e, this%x, this%rhs )
  if ( this%share_mpi_buf ) then
    n_shared_solver = n_shared_solver - 1
    if ( n_shared_solver == 0 ) deallocate( mpi_buf )
  else
    deallocate( this%sbuf, this%rbuf )
  endif
  
end subroutine destroy_pcr_penta

subroutine solve_pcr_penta( this )

  implicit none
  class( pcr_penta ), intent(inout) :: this

  integer :: step, i, comm_size
  integer :: stride1, stride2, pstride1, pstride2
  integer :: pid_lower1, pid_lower2, pid_upper1, pid_upper2
  integer :: sid, rid, ierr
  integer, dimension(2) :: i1_l, i2_l, i1_u, i2_u
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:), pointer :: a_tmp => null(), &
                                 b_tmp => null(), &
                                 c_tmp => null(), &
                                 d_tmp => null(), &
                                 e_tmp => null(), &
                                 rhs_tmp => null()
  real :: coefm1, coefp1, coefm2, coefp2, tmp
  real :: bm2, am1, dm2, cm1, ap1, em1, cp1, bp2, ep1, dp2, b0, d0
  integer :: im1, ip1, im2, ip2
  
  a_tmp   => this%sbuf(                    1 :     this%n_local )
  b_tmp   => this%sbuf(     this%n_local + 1 : 2 * this%n_local )
  c_tmp   => this%sbuf( 2 * this%n_local + 1 : 3 * this%n_local )
  d_tmp   => this%sbuf( 3 * this%n_local + 1 : 4 * this%n_local )
  e_tmp   => this%sbuf( 4 * this%n_local + 1 : 5 * this%n_local )
  rhs_tmp => this%sbuf( 5 * this%n_local + 1 : 6 * this%n_local )

  stride1 = 1; stride2 = 2
  do step = 1, this%tot_steps

    if ( step <= this%neighbor_steps - 1 ) then
      pstride1 = 1
      pstride2 = 1
      comm_size = stride1
      i1_l = (/ 1, comm_size /)
      i2_l = (/ comm_size + 1, 2 * comm_size /)
      i1_u = (/ this%n_local - comm_size + 1, this%n_local /)
      i2_u = (/ this%n_local - 2 * comm_size + 1, this%n_local - comm_size /)
    else if ( step == this%neighbor_steps ) then
      pstride1 = 1
      pstride2 = 2
      comm_size = stride1
      i1_l = (/ 1, this%n_local /)
      i2_l = (/ 1, this%n_local /)
      i1_u = (/ 1, this%n_local /)
      i2_u = (/ 1, this%n_local /)
    else
      pstride1 = pstride1 * 2
      pstride2 = pstride2 * 2
      comm_size = this%n_local
      i1_l = (/ 1, this%n_local /)
      i2_l = (/ 1, this%n_local /)
      i1_u = (/ 1, this%n_local /)
      i2_u = (/ 1, this%n_local /)
    endif

    pid_lower1 = this%myid - pstride1
    pid_upper1 = this%myid + pstride1
    pid_lower2 = this%myid - pstride2
    pid_upper2 = this%myid + pstride2

    ! send data to lower neighbor processor
    if ( pid_lower1 >= 0 ) then
      call this%pack_data( i1_l(1), i1_l(2) )
      call MPI_ISEND( this%sbuf, 6 * comm_size, this%mpi_dtype, pid_lower1, 1, this%mpi_comm, sid, ierr )
    endif
    if ( pid_upper1 <= this%num_procs - 1 ) then
      call MPI_RECV( this%rbuf, 6 * comm_size, this%mpi_dtype, pid_upper1, 1, this%mpi_comm, stat, ierr )
      call this%unpack_data( this%n_local + 1, this%n_local + comm_size )
    else
      call this%set_idle_ghost( this%n_local + 1, this%n_local + comm_size )
    endif
    if ( pid_lower1 >= 0 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    if ( pid_lower2 >= 0 ) then
      call this%pack_data( i2_l(1), i2_l(2) )
      call MPI_ISEND( this%sbuf, 6 * comm_size, this%mpi_dtype, pid_lower2, 1, this%mpi_comm, sid, ierr )
    endif
    if ( pid_upper2 <= this%num_procs - 1 ) then
      call MPI_RECV( this%rbuf, 6 * comm_size, this%mpi_dtype, pid_upper2, 1, this%mpi_comm, stat, ierr )
      call this%unpack_data( this%n_local + comm_size + 1, this%n_local + 2 * comm_size )
    else
      call this%set_idle_ghost( this%n_local + comm_size + 1, this%n_local + 2 * comm_size )
    endif
    if ( pid_lower2 >= 0 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! send data to upper neighbor processor
    if ( pid_upper1 <= this%num_procs - 1 ) then
      call this%pack_data( i1_u(1), i1_u(2) )
      call MPI_ISEND( this%sbuf, 6 * comm_size, this%mpi_dtype, pid_upper1, 1, this%mpi_comm, sid, ierr )
    endif
    if ( pid_lower1 >= 0 ) then
      call MPI_RECV( this%rbuf, 6 * comm_size, this%mpi_dtype, pid_lower1, 1, this%mpi_comm, stat, ierr )
      call this%unpack_data( 1 - comm_size, 0 )
    else
      call this%set_idle_ghost( 1 - comm_size, 0 )
    endif
    if ( pid_upper1 <= this%num_procs - 1 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    if ( pid_upper2 <= this%num_procs - 1 ) then
      call this%pack_data( i2_u(1), i2_u(2) )
      call MPI_ISEND( this%sbuf, 6 * comm_size, this%mpi_dtype, pid_upper2, 1, this%mpi_comm, sid, ierr )
    endif
    if ( pid_lower2 >= 0 ) then
      call MPI_RECV( this%rbuf, 6 * comm_size, this%mpi_dtype, pid_lower2, 1, this%mpi_comm, stat, ierr )
      call this%unpack_data( 1 - 2 * comm_size, -comm_size )
    else
      call this%set_idle_ghost( 1 - 2 * comm_size, -comm_size )
    endif
    if ( pid_upper2 <= this%num_procs - 1 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! use sbuf as the temporary arrays
    this%sbuf = 0.0
    do i = 1, this%n_local

      im1 = i - comm_size
      ip1 = i + comm_size
      im2 = i - 2 * comm_size
      ip2 = i + 2 * comm_size

      bm2 = this%b(im2); bp2 = this%b(ip2)
      dm2 = this%d(im2); dp2 = this%d(ip2)
      am1 = this%a(im1); ap1 = this%a(ip1)
      cm1 = this%c(im1); cp1 = this%c(ip1)
      em1 = this%e(im1); ep1 = this%e(ip1)
      d0  = this%d(i);   b0  = this%b(i)

      tmp = am1*cp1*dm2*dp2 - am1*bp2*dm2*ep1 + ap1*bm2*dp2*em1 + bm2*bp2*cm1*ep1 - bm2*cm1*cp1*dp2
      ! print *, "    tmp = ", tmp
      tmp = 1.0 / tmp 
      coefm2 =  am1 * ( ap1*d0*dp2 + b0*bp2*ep1 - b0*cp1*dp2 ) * tmp
      coefm1 = -bm2 * ( ap1*d0*dp2 + b0*bp2*ep1 - b0*cp1*dp2 ) * tmp
      coefp1 = -dp2 * ( am1*d0*dm2 + b0*bm2*em1 - d0*bm2*cm1 ) * tmp
      coefp2 =  ep1 * ( am1*d0*dm2 + b0*bm2*em1 - d0*bm2*cm1 ) * tmp

      a_tmp(i) = coefm2 * this%a(im2)
      b_tmp(i) = coefm2 * this%c(im2) + coefm1 * this%b(im1) + this%a(i)
      c_tmp(i) = coefm2 * this%e(im2) + coefm1 * this%d(im1) + this%c(i) + coefp1 * this%b(ip1) + coefp2 * this%a(ip2)
      d_tmp(i) =                                               this%e(i) + coefp1 * this%d(ip1) + coefp2 * this%c(ip2)
      e_tmp(i) =                                                                                + coefp2 * this%e(ip2)
      rhs_tmp(i) = coefm2 * this%rhs(im2) + coefm1 * this%rhs(im1) + this%rhs(i) + coefp1 * this%rhs(ip1) + coefp2 * this%rhs(ip2)

    enddo

    this%a(1:this%n_local) = a_tmp(1:this%n_local)
    this%b(1:this%n_local) = b_tmp(1:this%n_local)
    this%c(1:this%n_local) = c_tmp(1:this%n_local)
    this%d(1:this%n_local) = d_tmp(1:this%n_local)
    this%e(1:this%n_local) = e_tmp(1:this%n_local)
    this%rhs(1:this%n_local) = rhs_tmp(1:this%n_local)

    stride1 = stride1 * 2
    stride2 = stride2 * 2

  enddo

  do i = 1, this%n_local
    this%x(i) = this%rhs(i) / this%c(i)
  enddo

end subroutine solve_pcr_penta

subroutine pack_data_pcr_penta( this, i0, i1 )

  implicit none
  class( pcr_penta ), intent(inout) :: this
  integer, intent(in) :: i0, i1

  integer :: i, num
  num = i1 - i0 + 1
  this%sbuf(      1:  num) = this%a(i0:i1)
  this%sbuf(  num+1:2*num) = this%b(i0:i1)
  this%sbuf(2*num+1:3*num) = this%c(i0:i1)
  this%sbuf(3*num+1:4*num) = this%d(i0:i1)
  this%sbuf(4*num+1:5*num) = this%e(i0:i1)
  this%sbuf(5*num+1:6*num) = this%rhs(i0:i1)

end subroutine pack_data_pcr_penta

subroutine unpack_data_pcr_penta( this, i0, i1 )

  implicit none
  class( pcr_penta ), intent(inout) :: this
  integer, intent(in) :: i0, i1

  integer :: i, num
  num = i1 - i0 + 1
  this%a(i0:i1)   = this%rbuf(      1:  num)
  this%b(i0:i1)   = this%rbuf(  num+1:2*num)
  this%c(i0:i1)   = this%rbuf(2*num+1:3*num)
  this%d(i0:i1)   = this%rbuf(3*num+1:4*num)
  this%e(i0:i1)   = this%rbuf(4*num+1:5*num)
  this%rhs(i0:i1) = this%rbuf(5*num+1:6*num)

end subroutine unpack_data_pcr_penta

subroutine set_idle_ghost_pcr_penta( this, i0, i1 )

  implicit none
  class( pcr_penta ), intent(inout) :: this
  integer, intent(in) :: i0, i1

  this%a(i0:i1) = 0.0
  this%b(i0:i1) = 1.0
  this%c(i0:i1) = 2.0
  this%d(i0:i1) = 1.0
  this%e(i0:i1) = 0.0
  this%rhs(i0:i1) = 0.0

end subroutine set_idle_ghost_pcr_penta

subroutine set_values_rhs_pcr_penta( this, rhs )

  implicit none
  class( pcr_penta ), intent(inout) :: this
  real, intent(in), dimension(:) :: rhs
  integer :: i

  do i = 1, this%n_local
    this%rhs(i) = rhs(i)
  enddo

end subroutine set_values_rhs_pcr_penta

subroutine set_values_matrix_pcr_penta( this, a, b, c, d, e )

  implicit none
  class( pcr_penta ), intent(inout) :: this
  real, intent(in), dimension(:) :: a, b, c, d, e
  integer :: i

  do i = 1, this%n_local
    this%a(i) = a(i)
    this%b(i) = b(i)
    this%c(i) = c(i)
    this%d(i) = d(i)
    this%e(i) = e(i)
  enddo

  if (this%myid == 0) then
    this%a(1) = epsilon(1.0)
    this%b(1) = epsilon(1.0)
    this%a(2) = epsilon(1.0)
  endif
  if (this%myid == this%num_procs-1) then
    this%e(this%n_local) = epsilon(1.0)
    this%d(this%n_local) = epsilon(1.0)
    this%e(this%n_local-1) = epsilon(1.0)
  endif

end subroutine set_values_matrix_pcr_penta

subroutine get_values_x_pcr_penta( this, x )

  implicit none
  class( pcr_penta ), intent(inout) :: this
  real, intent(out), dimension(:) :: x
  integer :: i

  do i = 1, this%n_local
    x(i) = this%x(i)
  enddo

end subroutine get_values_x_pcr_penta

subroutine print_x_pcr_penta( this, filename )

  implicit none
  class( pcr_penta ), intent(inout) :: this
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

end subroutine print_x_pcr_penta

subroutine print_rhs_pcr_penta( this, filename )

  implicit none
  class( pcr_penta ), intent(inout) :: this
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

end subroutine print_rhs_pcr_penta

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

end module pcr_penta_class