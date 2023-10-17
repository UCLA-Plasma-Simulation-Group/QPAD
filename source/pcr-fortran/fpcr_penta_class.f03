module fpcr_penta_class

use mpi

implicit none
private

type :: fpcr_penta

  private

  ! MPI information
  integer :: mpi_comm, myid, num_procs, mpi_dtype

  ! Global problem size and number of rows in this partition
  integer :: n_global, n_local

  ! a, b, c are the lower, central and upper diagonal elements
  ! rhs is the right-hand-side of the linear system
  real, dimension(:), allocatable :: a, b, c, d, e, rhs, x

  ! total steps of cyclic reduction
  integer :: tot_steps

  ! steps of cyclic redution that only needs neighbor processors' data
  integer :: neighbor_steps

  ! CR coefficients
  real, dimension(:,:), allocatable :: coefm1, coefp1, coefm2, coefp2

  contains

  procedure :: create            => create_fpcr_penta
  procedure :: destroy           => destroy_fpcr_penta
  procedure :: solve             => solve_fpcr_penta
  generic   :: set_values_rhs    => set_values_rhs_fpcr_penta, set_values_rhs_byrow_fpcr_penta
  generic   :: set_values_matrix => set_values_matrix_fpcr_penta, set_values_matrix_byrow_fpcr_penta
  generic   :: get_values_matrix => get_values_matrix_fpcr_penta, get_values_matrix_byrow_fpcr_penta
  generic   :: get_values_x      => get_values_x_fpcr_penta, get_values_x_byrow_fpcr_penta
  procedure :: print_x           => print_x_fpcr_penta
  procedure :: print_rhs         => print_rhs_fpcr_penta
  procedure :: print_matrix      => print_matrix_fpcr_penta
  procedure :: generate_cr_coef  => generate_cr_coef_fpcr_penta
  procedure, private :: set_values_rhs_fpcr_penta
  procedure, private :: set_values_rhs_byrow_fpcr_penta
  procedure, private :: set_values_matrix_fpcr_penta
  procedure, private :: set_values_matrix_byrow_fpcr_penta
  procedure, private :: get_values_matrix_fpcr_penta
  procedure, private :: get_values_matrix_byrow_fpcr_penta
  procedure, private :: get_values_x_fpcr_penta
  procedure, private :: get_values_x_byrow_fpcr_penta

end type fpcr_penta

public :: fpcr_penta

integer, parameter :: stdout = 6, stderr = 0

contains

subroutine create_fpcr_penta( this, n_global, mpi_comm, err )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
  integer, intent(in) :: n_global, mpi_comm
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
  allocate( this%a( this%n_local ) ); this%a = 0.0
  allocate( this%b( this%n_local ) ); this%b = 0.0
  allocate( this%c( this%n_local ) ); this%c = 0.0
  allocate( this%d( this%n_local ) ); this%d = 0.0
  allocate( this%e( this%n_local ) ); this%e = 0.0
  allocate( this%x( this%n_local ) ); this%x = 0.0
  allocate( this%rhs(1-2*this%n_local:3*this%n_local) ); this%rhs = 0.0

  allocate( this%coefm1( this%n_local, this%tot_steps ) )
  allocate( this%coefm2( this%n_local, this%tot_steps ) )
  allocate( this%coefp1( this%n_local, this%tot_steps ) )
  allocate( this%coefp2( this%n_local, this%tot_steps ) )

  if ( digits(1.0) > 24 ) then
    this%mpi_dtype = MPI_DOUBLE_PRECISION
  else
    this%mpi_dtype = MPI_REAL
  endif

  err = 0

end subroutine create_fpcr_penta

subroutine destroy_fpcr_penta( this )

  implicit none
  class( fpcr_penta ), intent(inout) :: this

  if ( allocated(this%a) ) deallocate(this%a)
  if ( allocated(this%b) ) deallocate(this%b)
  if ( allocated(this%c) ) deallocate(this%c)
  if ( allocated(this%d) ) deallocate(this%d)
  if ( allocated(this%e) ) deallocate(this%e)
  if ( allocated(this%x) ) deallocate(this%x)
  if ( allocated(this%rhs) ) deallocate(this%rhs)

end subroutine destroy_fpcr_penta

subroutine solve_fpcr_penta( this )

  implicit none
  class( fpcr_penta ), intent(inout) :: this

  integer :: step, i, comm_size
  integer :: stride1, stride2, pstride1, pstride2
  integer :: pid_lower1, pid_lower2, pid_upper1, pid_upper2
  integer :: sid, ierr
  integer, dimension(2) :: i1_l, i2_l, i1_u, i2_u
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer :: im1, ip1, im2, ip2, ri0, ri1

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
      call MPI_ISEND( this%rhs( i1_l(1) : i1_l(2) ), comm_size, this%mpi_dtype, pid_lower1, 1, this%mpi_comm, sid, ierr )
    endif
    ri0 = this%n_local + 1
    ri1 = this%n_local + comm_size
    if ( pid_upper1 <= this%num_procs - 1 ) then
      call MPI_RECV( this%rhs(ri0:ri1), comm_size, this%mpi_dtype, pid_upper1, 1, this%mpi_comm, stat, ierr )
    else
      this%rhs(ri0:ri1) = 0.0
    endif
    if ( pid_lower1 >= 0 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    if ( pid_lower2 >= 0 ) then
      call MPI_ISEND( this%rhs( i2_l(1) : i2_l(2) ), comm_size, this%mpi_dtype, pid_lower2, 1, this%mpi_comm, sid, ierr )
    endif
    ri0 = this%n_local + comm_size + 1
    ri1 = this%n_local + 2 * comm_size
    if ( pid_upper2 <= this%num_procs - 1 ) then
      call MPI_RECV( this%rhs(ri0:ri1), comm_size, this%mpi_dtype, pid_upper2, 1, this%mpi_comm, stat, ierr )
    else
      this%rhs(ri0:ri1) = 0.0
    endif
    if ( pid_lower2 >= 0 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! send data to upper neighbor processor
    if ( pid_upper1 <= this%num_procs - 1 ) then
      call MPI_ISEND( this%rhs( i1_u(1) : i1_u(2) ), comm_size, this%mpi_dtype, pid_upper1, 1, this%mpi_comm, sid, ierr )
    endif
    ri0 = 1 - comm_size
    ri1 = 0
    if ( pid_lower1 >= 0 ) then
      call MPI_RECV( this%rhs(ri0:ri1), comm_size, this%mpi_dtype, pid_lower1, 1, this%mpi_comm, stat, ierr )
    else
      this%rhs(ri0:ri1) = 0.0
    endif
    if ( pid_upper1 <= this%num_procs - 1 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    if ( pid_upper2 <= this%num_procs - 1 ) then
      call MPI_ISEND( this%rhs( i2_u(1) : i2_u(2) ), comm_size, this%mpi_dtype, pid_upper2, 1, this%mpi_comm, sid, ierr )
    endif
    ri0 = 1 - 2 * comm_size
    ri1 = -comm_size
    if ( pid_lower2 >= 0 ) then
      call MPI_RECV( this%rhs(ri0:ri1), comm_size, this%mpi_dtype, pid_lower2, 1, this%mpi_comm, stat, ierr )
    else
      this%rhs(ri0:ri1) = 0.0
    endif
    if ( pid_upper2 <= this%num_procs - 1 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! calculate
    ! use x as the temporary arrays
    this%x = 0.0
    do i = 1, this%n_local

      im1 = i - comm_size
      ip1 = i + comm_size
      im2 = i - 2 * comm_size
      ip2 = i + 2 * comm_size

      ! use x as the temporary arrays
      this%x(i) = this%rhs(i)
      this%x(i) = this%x(i) + this%coefm1(i,step) * this%rhs(im1)
      this%x(i) = this%x(i) + this%coefm2(i,step) * this%rhs(im2)
      this%x(i) = this%x(i) + this%coefp1(i,step) * this%rhs(ip1)
      this%x(i) = this%x(i) + this%coefp2(i,step) * this%rhs(ip2)

    enddo

    this%rhs(1:this%n_local) = this%x

    stride1 = stride1 * 2
    stride2 = stride2 * 2

  enddo

  do i = 1, this%n_local
    this%x(i) = this%rhs(i) / this%c(i)
  enddo

end subroutine solve_fpcr_penta

subroutine generate_cr_coef_fpcr_penta( this )

  implicit none
  class( fpcr_penta ), intent(inout) :: this

  integer :: step, i, comm_size
  integer :: stride1, stride2, pstride1, pstride2
  integer :: pid_lower1, pid_lower2, pid_upper1, pid_upper2
  integer :: sid1, sid2, sid3, sid4, sid5, ierr
  integer, dimension(2) :: i1_l, i2_l, i1_u, i2_u
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:), allocatable :: a, b, c, d, e
  real :: tmp
  real :: am1, bm1, cm1, dm1, em1
  real :: am2, bm2, cm2, dm2, em2
  real :: ap1, bp1, cp1, dp1, ep1
  real :: ap2, bp2, cp2, dp2, ep2
  real :: a0,  b0,  c0,  d0,  e0
  integer :: im1, ip1, im2, ip2, ri0, ri1
  
  allocate( a( 1-2*this%n_local:3*this%n_local ) ); a = 0.0
  allocate( b( 1-2*this%n_local:3*this%n_local ) ); b = 0.0
  allocate( c( 1-2*this%n_local:3*this%n_local ) ); c = 0.0
  allocate( d( 1-2*this%n_local:3*this%n_local ) ); d = 0.0
  allocate( e( 1-2*this%n_local:3*this%n_local ) ); e = 0.0
  a(1:this%n_local) = this%a
  b(1:this%n_local) = this%b
  c(1:this%n_local) = this%c
  d(1:this%n_local) = this%d
  e(1:this%n_local) = this%e

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
      call MPI_ISEND( a( i1_l(1):i1_l(2) ), comm_size, this%mpi_dtype, pid_lower1, 1, this%mpi_comm, sid1, ierr )
      call MPI_ISEND( b( i1_l(1):i1_l(2) ), comm_size, this%mpi_dtype, pid_lower1, 2, this%mpi_comm, sid2, ierr )
      call MPI_ISEND( c( i1_l(1):i1_l(2) ), comm_size, this%mpi_dtype, pid_lower1, 3, this%mpi_comm, sid3, ierr )
      call MPI_ISEND( d( i1_l(1):i1_l(2) ), comm_size, this%mpi_dtype, pid_lower1, 4, this%mpi_comm, sid4, ierr )
      call MPI_ISEND( e( i1_l(1):i1_l(2) ), comm_size, this%mpi_dtype, pid_lower1, 5, this%mpi_comm, sid5, ierr )
    endif
    ri0 = this%n_local + 1
    ri1 = this%n_local + comm_size
    if ( pid_upper1 <= this%num_procs - 1 ) then
      call MPI_RECV( a(ri0:ri1), comm_size, this%mpi_dtype, pid_upper1, 1, this%mpi_comm, stat, ierr )
      call MPI_RECV( b(ri0:ri1), comm_size, this%mpi_dtype, pid_upper1, 2, this%mpi_comm, stat, ierr )
      call MPI_RECV( c(ri0:ri1), comm_size, this%mpi_dtype, pid_upper1, 3, this%mpi_comm, stat, ierr )
      call MPI_RECV( d(ri0:ri1), comm_size, this%mpi_dtype, pid_upper1, 4, this%mpi_comm, stat, ierr )
      call MPI_RECV( e(ri0:ri1), comm_size, this%mpi_dtype, pid_upper1, 5, this%mpi_comm, stat, ierr )
    else
      a(ri0:ri1) = 0.0
      b(ri0:ri1) = 1.0
      c(ri0:ri1) = 2.0
      d(ri0:ri1) = 1.0
      e(ri0:ri1) = 0.0
    endif
    if ( pid_lower1 >= 0 ) then
      call MPI_WAIT( sid1, stat, ierr )
      call MPI_WAIT( sid2, stat, ierr )
      call MPI_WAIT( sid3, stat, ierr )
      call MPI_WAIT( sid4, stat, ierr )
      call MPI_WAIT( sid5, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    if ( pid_lower2 >= 0 ) then
      call MPI_ISEND( a( i2_l(1):i2_l(2) ), comm_size, this%mpi_dtype, pid_lower2, 1, this%mpi_comm, sid1, ierr )
      call MPI_ISEND( b( i2_l(1):i2_l(2) ), comm_size, this%mpi_dtype, pid_lower2, 2, this%mpi_comm, sid2, ierr )
      call MPI_ISEND( c( i2_l(1):i2_l(2) ), comm_size, this%mpi_dtype, pid_lower2, 3, this%mpi_comm, sid3, ierr )
      call MPI_ISEND( d( i2_l(1):i2_l(2) ), comm_size, this%mpi_dtype, pid_lower2, 4, this%mpi_comm, sid4, ierr )
      call MPI_ISEND( e( i2_l(1):i2_l(2) ), comm_size, this%mpi_dtype, pid_lower2, 5, this%mpi_comm, sid5, ierr )
    endif
    ri0 = this%n_local + comm_size + 1
    ri1 = this%n_local + 2 * comm_size
    if ( pid_upper2 <= this%num_procs - 1 ) then
      call MPI_RECV( a(ri0:ri1), comm_size, this%mpi_dtype, pid_upper2, 1, this%mpi_comm, stat, ierr )
      call MPI_RECV( b(ri0:ri1), comm_size, this%mpi_dtype, pid_upper2, 2, this%mpi_comm, stat, ierr )
      call MPI_RECV( c(ri0:ri1), comm_size, this%mpi_dtype, pid_upper2, 3, this%mpi_comm, stat, ierr )
      call MPI_RECV( d(ri0:ri1), comm_size, this%mpi_dtype, pid_upper2, 4, this%mpi_comm, stat, ierr )
      call MPI_RECV( e(ri0:ri1), comm_size, this%mpi_dtype, pid_upper2, 5, this%mpi_comm, stat, ierr )
    else
      a(ri0:ri1) = 0.0
      b(ri0:ri1) = 1.0
      c(ri0:ri1) = 2.0
      d(ri0:ri1) = 1.0
      e(ri0:ri1) = 0.0
    endif
    if ( pid_lower2 >= 0 ) then
      call MPI_WAIT( sid1, stat, ierr )
      call MPI_WAIT( sid2, stat, ierr )
      call MPI_WAIT( sid3, stat, ierr )
      call MPI_WAIT( sid4, stat, ierr )
      call MPI_WAIT( sid5, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! send data to upper neighbor processor
    if ( pid_upper1 <= this%num_procs - 1 ) then
      call MPI_ISEND( a( i1_u(1):i1_u(2) ), comm_size, this%mpi_dtype, pid_upper1, 1, this%mpi_comm, sid1, ierr )
      call MPI_ISEND( b( i1_u(1):i1_u(2) ), comm_size, this%mpi_dtype, pid_upper1, 2, this%mpi_comm, sid2, ierr )
      call MPI_ISEND( c( i1_u(1):i1_u(2) ), comm_size, this%mpi_dtype, pid_upper1, 3, this%mpi_comm, sid3, ierr )
      call MPI_ISEND( d( i1_u(1):i1_u(2) ), comm_size, this%mpi_dtype, pid_upper1, 4, this%mpi_comm, sid4, ierr )
      call MPI_ISEND( e( i1_u(1):i1_u(2) ), comm_size, this%mpi_dtype, pid_upper1, 5, this%mpi_comm, sid5, ierr )
    endif
    ri0 = 1 - comm_size
    ri1 = 0
    if ( pid_lower1 >= 0 ) then
      call MPI_RECV( a(ri0:ri1), comm_size, this%mpi_dtype, pid_lower1, 1, this%mpi_comm, stat, ierr )
      call MPI_RECV( b(ri0:ri1), comm_size, this%mpi_dtype, pid_lower1, 2, this%mpi_comm, stat, ierr )
      call MPI_RECV( c(ri0:ri1), comm_size, this%mpi_dtype, pid_lower1, 3, this%mpi_comm, stat, ierr )
      call MPI_RECV( d(ri0:ri1), comm_size, this%mpi_dtype, pid_lower1, 4, this%mpi_comm, stat, ierr )
      call MPI_RECV( e(ri0:ri1), comm_size, this%mpi_dtype, pid_lower1, 5, this%mpi_comm, stat, ierr )
    else
      a(ri0:ri1) = 0.0
      b(ri0:ri1) = 1.0
      c(ri0:ri1) = 2.0
      d(ri0:ri1) = 1.0
      e(ri0:ri1) = 0.0
    endif
    if ( pid_upper1 <= this%num_procs - 1 ) then
      call MPI_WAIT( sid1, stat, ierr )
      call MPI_WAIT( sid2, stat, ierr )
      call MPI_WAIT( sid3, stat, ierr )
      call MPI_WAIT( sid4, stat, ierr )
      call MPI_WAIT( sid5, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    if ( pid_upper2 <= this%num_procs - 1 ) then
      call MPI_ISEND( a( i2_u(1):i2_u(2) ), comm_size, this%mpi_dtype, pid_upper2, 1, this%mpi_comm, sid1, ierr )
      call MPI_ISEND( b( i2_u(1):i2_u(2) ), comm_size, this%mpi_dtype, pid_upper2, 2, this%mpi_comm, sid2, ierr )
      call MPI_ISEND( c( i2_u(1):i2_u(2) ), comm_size, this%mpi_dtype, pid_upper2, 3, this%mpi_comm, sid3, ierr )
      call MPI_ISEND( d( i2_u(1):i2_u(2) ), comm_size, this%mpi_dtype, pid_upper2, 4, this%mpi_comm, sid4, ierr )
      call MPI_ISEND( e( i2_u(1):i2_u(2) ), comm_size, this%mpi_dtype, pid_upper2, 5, this%mpi_comm, sid5, ierr )
    endif
    ri0 = 1 - 2 * comm_size
    ri1 = -comm_size
    if ( pid_lower2 >= 0 ) then
      call MPI_RECV( a(ri0:ri1), comm_size, this%mpi_dtype, pid_lower2, 1, this%mpi_comm, stat, ierr )
      call MPI_RECV( b(ri0:ri1), comm_size, this%mpi_dtype, pid_lower2, 2, this%mpi_comm, stat, ierr )
      call MPI_RECV( c(ri0:ri1), comm_size, this%mpi_dtype, pid_lower2, 3, this%mpi_comm, stat, ierr )
      call MPI_RECV( d(ri0:ri1), comm_size, this%mpi_dtype, pid_lower2, 4, this%mpi_comm, stat, ierr )
      call MPI_RECV( e(ri0:ri1), comm_size, this%mpi_dtype, pid_lower2, 5, this%mpi_comm, stat, ierr )
    else
      a(ri0:ri1) = 0.0
      b(ri0:ri1) = 1.0
      c(ri0:ri1) = 2.0
      d(ri0:ri1) = 1.0
      e(ri0:ri1) = 0.0
    endif
    if ( pid_upper2 <= this%num_procs - 1 ) then
      call MPI_WAIT( sid1, stat, ierr )
      call MPI_WAIT( sid2, stat, ierr )
      call MPI_WAIT( sid3, stat, ierr )
      call MPI_WAIT( sid4, stat, ierr )
      call MPI_WAIT( sid5, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! calculate the CR coefficient
    do i = 1, this%n_local

      im1 = i - comm_size
      ip1 = i + comm_size
      im2 = i - 2 * comm_size
      ip2 = i + 2 * comm_size

      am1 = a(im1); bm1 = b(im1); cm1 = c(im1); dm1 = d(im1); em1 = e(im1)
      am2 = a(im2); bm2 = b(im2); cm2 = c(im2); dm2 = d(im2); em2 = e(im2)
      ap1 = a(ip1); bp1 = b(ip1); cp1 = c(ip1); dp1 = d(ip1); ep1 = e(ip1)
      ap2 = a(ip2); bp2 = b(ip2); cp2 = c(ip2); dp2 = d(ip2); ep2 = e(ip2)
      a0 = a(i);    b0  = b(i);   c0  = c(i);   d0  = d(i);   e0  = e(i)

      tmp = am1*cp1*dm2*dp2 - am1*bp2*dm2*ep1 + ap1*bm2*dp2*em1 + bm2*bp2*cm1*ep1 - bm2*cm1*cp1*dp2
      ! print *, "    tmp = ", tmp

      if( abs(tmp) > 1e-40 ) then 
         tmp = 1.0 / tmp
         this%coefm2(i, step) =  am1 * ( ap1*d0*dp2 + b0*bp2*ep1 - b0*cp1*dp2 ) * tmp
         this%coefm1(i, step) = -bm2 * ( ap1*d0*dp2 + b0*bp2*ep1 - b0*cp1*dp2 ) * tmp
         this%coefp1(i, step) = -dp2 * ( am1*d0*dm2 + b0*bm2*em1 - d0*bm2*cm1 ) * tmp
         this%coefp2(i, step) =  ep1 * ( am1*d0*dm2 + b0*bm2*em1 - d0*bm2*cm1 ) * tmp
      else
         this%coefm2(i, step) = 0
         this%coefm1(i, step) = 0
         this%coefp1(i, step) = 0
         this%coefp2(i, step) = 0
      endif
         
      ! use class members as temporary arrays
      this%a(i) = this%coefm2(i,step) * am2
      this%b(i) = this%coefm2(i,step) * cm2 + this%coefm1(i,step) * bm1 + a0
      this%c(i) = this%coefm2(i,step) * em2 + this%coefm1(i,step) * dm1 + c0 + this%coefp1(i,step) * bp1 + this%coefp2(i,step) * ap2
      this%d(i) =                                                         e0 + this%coefp1(i,step) * dp1 + this%coefp2(i,step) * cp2
      this%e(i) =                                                                                          this%coefp2(i,step) * ep2

    enddo

    a(1:this%n_local) = this%a
    b(1:this%n_local) = this%b
    c(1:this%n_local) = this%c
    d(1:this%n_local) = this%d
    e(1:this%n_local) = this%e

    stride1 = stride1 * 2
    stride2 = stride2 * 2

  enddo

  this%c = c(1:this%n_local)
  deallocate( a, b, c, d, e )

  ! off-diagonal elements are all zeros, can be deleted
  ! deallocate( this%a, this%b, this%d, this%e )

end subroutine generate_cr_coef_fpcr_penta

subroutine set_values_rhs_fpcr_penta( this, rhs )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
  real, intent(in), dimension(:) :: rhs
  integer :: i

  do i = 1, this%n_local
    this%rhs(i) = rhs(i)
  enddo

end subroutine set_values_rhs_fpcr_penta

subroutine set_values_rhs_byrow_fpcr_penta( this, rhs, i_loc )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
  real, intent(in) :: rhs
  integer, intent(in) :: i_loc
  this%rhs(i_loc) = rhs
  
end subroutine set_values_rhs_byrow_fpcr_penta

subroutine set_values_matrix_fpcr_penta( this, a, b, c, d, e )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
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

end subroutine set_values_matrix_fpcr_penta

subroutine set_values_matrix_byrow_fpcr_penta( this, a, b, c, d, e, i_loc )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
  real, intent(in) :: a, b, c, d, e
  integer, intent(in) :: i_loc
  integer :: i

  this%a(i_loc) = a
  this%b(i_loc) = b
  this%c(i_loc) = c
  this%d(i_loc) = d
  this%e(i_loc) = e

  if ( this%myid == 0 .and. i_loc == 1 ) then
    this%a(i_loc) = epsilon(1.0)
    this%b(i_loc) = epsilon(1.0)
  else if ( this%myid == 0 .and. i_loc == 2 ) then
    this%a(i_loc) = epsilon(1.0)
  endif
  if ( this%myid == this%num_procs-1 .and. i_loc == this%n_local ) then
    this%e(i_loc) = epsilon(1.0)
    this%d(i_loc) = epsilon(1.0)
  else if ( this%myid == this%num_procs-1 .and. i_loc == this%n_local - 1 ) then
    this%e(i_loc) = epsilon(1.0)
  endif

end subroutine set_values_matrix_byrow_fpcr_penta

subroutine get_values_matrix_fpcr_penta( this, a, b, c, d, e )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
  real, intent(out), dimension(:) :: a, b, c, d, e
  integer :: i

  do i = 1, this%n_local
    a(i) = this%a(i)
    b(i) = this%b(i)
    c(i) = this%c(i)
    d(i) = this%d(i)
    e(i) = this%e(i)
  enddo

end subroutine get_values_matrix_fpcr_penta

subroutine get_values_matrix_byrow_fpcr_penta( this, a, b, c, d, e, i_loc )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
  real, intent(out) :: a, b, c, d, e
  integer, intent(in) :: i_loc
  integer :: i

  a = this%a(i_loc)
  b = this%b(i_loc)
  c = this%c(i_loc)
  d = this%d(i_loc)
  e = this%e(i_loc)

end subroutine get_values_matrix_byrow_fpcr_penta

subroutine get_values_x_fpcr_penta( this, x )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
  real, intent(out), dimension(:) :: x
  integer :: i

  do i = 1, this%n_local
    x(i) = this%x(i)
  enddo

end subroutine get_values_x_fpcr_penta

subroutine get_values_x_byrow_fpcr_penta( this, x, i_loc )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
  real, intent(out) :: x
  integer, intent(in) :: i_loc
  x = this%x(i_loc)

end subroutine get_values_x_byrow_fpcr_penta

subroutine print_x_fpcr_penta( this, filename )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
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

end subroutine print_x_fpcr_penta

subroutine print_rhs_fpcr_penta( this, filename )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
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

end subroutine print_rhs_fpcr_penta

subroutine print_matrix_fpcr_penta( this, filename )

  implicit none
  class( fpcr_penta ), intent(inout) :: this
  character(len=*), intent(in) :: filename
  integer :: unit, i
  character(len=128) :: fstr

  write( fstr, '(I0.4)' ) this%myid
  fstr = trim(filename) // '.' // trim(fstr)
  unit = 2
  open( unit, file=trim(fstr) )
  do i = 1, this%n_local
    write( unit, '(5ES30.15)' ) this%a(i), this%b(i), this%c(i), this%d(i), this%e(i)
  enddo

end subroutine print_matrix_fpcr_penta

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

end module fpcr_penta_class
