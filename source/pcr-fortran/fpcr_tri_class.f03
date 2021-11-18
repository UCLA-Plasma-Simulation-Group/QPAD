module fpcr_tri_class

use mpi

implicit none
private

type :: fpcr_tri

  private

  ! MPI information
  integer :: mpi_comm, myid, num_procs, mpi_dtype

  ! Global problem size and number of rows in this partition
  integer :: n_global, n_local

  ! a, b, c are the lower, central and upper diagonal elements
  ! rhs is the right-hand-side of the linear system
  real, dimension(:), allocatable :: a, b, c, rhs, x

  ! total steps of cyclic reduction
  integer :: tot_steps

  ! steps of cyclic redution that only needs neighbor processors' data
  integer :: neighbor_steps

  ! CR coefficients
  real, dimension(:,:), allocatable :: coefm, coefp

  contains

  procedure :: create            => create_fpcr_tri
  procedure :: destroy           => destroy_fpcr_tri
  procedure :: solve             => solve_fpcr_tri
  procedure :: set_values_rhs    => set_values_rhs_fpcr_tri
  generic   :: set_values_matrix => set_values_matrix_fpcr_tri, set_values_matrix_byrow_fpcr_tri
  procedure :: get_values_x      => get_values_x_fpcr_tri
  procedure :: print_x           => print_x_fpcr_tri
  procedure :: print_rhs         => print_rhs_fpcr_tri
  procedure :: generate_cr_coef  => generate_cr_coef_fpcr_tri
  procedure, private :: set_values_matrix_fpcr_tri
  procedure, private :: set_values_matrix_byrow_fpcr_tri

end type fpcr_tri

public :: fpcr_tri

integer, parameter :: stdout = 6, stderr = 0

contains

subroutine create_fpcr_tri( this, n_global, mpi_comm, err )

  implicit none
  class( fpcr_tri ), intent(inout) :: this
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
  allocate( this%x( this%n_local ) ); this%x = 0.0
  allocate( this%rhs(1-this%n_local:2*this%n_local) ); this%rhs = 0.0

  allocate( this%coefm( this%n_local, this%tot_steps ) )
  allocate( this%coefp( this%n_local, this%tot_steps ) )

  if ( digits(1.0) > 24 ) then
    this%mpi_dtype = MPI_DOUBLE_PRECISION
  else
    this%mpi_dtype = MPI_REAL
  endif

  err = 0

end subroutine create_fpcr_tri

subroutine destroy_fpcr_tri( this )

  implicit none
  class( fpcr_tri ), intent(inout) :: this

  if ( allocated(this%a) ) deallocate(this%a)
  if ( allocated(this%b) ) deallocate(this%b)
  if ( allocated(this%c) ) deallocate(this%c)
  if ( allocated(this%x) ) deallocate(this%x)
  if ( allocated(this%rhs) ) deallocate(this%rhs)

end subroutine destroy_fpcr_tri

subroutine solve_fpcr_tri( this )

  implicit none
  class( fpcr_tri ), intent(inout) :: this

  integer :: step, stride, pstride, i, pid_lower, pid_upper, comm_size
  integer :: sid, ierr
  integer, dimension(2) :: il, iu
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer :: im, ip, ri0, ri1

  stride = 1; pstride = 1
  do step = 1, this%tot_steps

    if ( step <= this%neighbor_steps ) then
      pstride =  1
      comm_size = stride
      il = (/ 1, comm_size /)
      iu = (/ this%n_local - comm_size + 1, this%n_local /)
    else
      pstride = pstride * 2
      comm_size = this%n_local
      il = (/ 1, this%n_local /)
      iu = (/ 1, this%n_local /)
    endif

    pid_lower = this%myid - pstride
    pid_upper = this%myid + pstride

    ! send data to lower neighbor processor
    if ( pid_lower >= 0 ) then
      call MPI_ISEND( this%rhs( il(1):il(2) ), comm_size, this%mpi_dtype, pid_lower, 1, this%mpi_comm, sid, ierr )
    endif
    ri0 = this%n_local + 1
    ri1 = this%n_local + comm_size
    if ( pid_upper <= this%num_procs - 1 ) then
      call MPI_RECV( this%rhs(ri0:ri1), comm_size, this%mpi_dtype, pid_upper, 1, this%mpi_comm, stat, ierr )
    else
      this%rhs(ri0:ri1) = 0.0
    endif
    if ( pid_lower >= 0 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! send data to upper neighbor processor
    if ( pid_upper <= this%num_procs - 1 ) then
      call MPI_ISEND( this%rhs( iu(1):iu(2) ), comm_size, this%mpi_dtype, pid_upper, 1, this%mpi_comm, sid, ierr )
    endif
    ri0 = 1 - comm_size
    ri1 = 0
    if ( pid_lower >= 0 ) then
      call MPI_RECV( this%rhs(ri0:ri1), comm_size, this%mpi_dtype, pid_lower, 1, this%mpi_comm, stat, ierr )
    else
      this%rhs(ri0:ri1) = 0.0
    endif
    if ( pid_upper <= this%num_procs - 1 ) then
      call MPI_WAIT( sid, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! calculate
    ! use x as the temporary arrays
    this%x = 0.0
    do i = 1, this%n_local

      im = i - comm_size
      ip = i + comm_size
      
      ! use x as the temporary arrays
      this%x(i) = this%rhs(i)
      this%x(i) = this%x(i) + this%coefm(i,step) * this%rhs(im)
      this%x(i) = this%x(i) + this%coefp(i,step) * this%rhs(ip)

    enddo

    this%rhs(1:this%n_local) = this%x

    stride = stride * 2

  enddo

  do i = 1, this%n_local
    this%x(i) = this%rhs(i) / this%b(i)
  enddo

end subroutine solve_fpcr_tri

subroutine generate_cr_coef_fpcr_tri( this )

  implicit none
  class( fpcr_tri ), intent(inout) :: this

  integer :: step, i, comm_size
  integer :: stride, pstride
  integer :: pid_lower, pid_upper
  integer :: sid1, sid2, sid3, ierr
  integer, dimension(2) :: il, iu
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, dimension(:), allocatable :: a, b, c
  integer :: im, ip, ri0, ri1
  
  allocate( a( 1-this%n_local:2*this%n_local ) ); a = 0.0
  allocate( b( 1-this%n_local:2*this%n_local ) ); b = 0.0
  allocate( c( 1-this%n_local:2*this%n_local ) ); c = 0.0
  a(1:this%n_local) = this%a
  b(1:this%n_local) = this%b
  c(1:this%n_local) = this%c

  stride = 1
  do step = 1, this%tot_steps

    if ( step <= this%neighbor_steps ) then
      pstride =  1
      comm_size = stride
      il = (/ 1, comm_size /)
      iu = (/ this%n_local - comm_size + 1, this%n_local /)
    else
      pstride = pstride * 2
      comm_size = this%n_local
      il = (/ 1, this%n_local /)
      iu = (/ 1, this%n_local /)
    endif

    pid_lower = this%myid - pstride
    pid_upper = this%myid + pstride

    ! send data to lower neighbor processor
    if ( pid_lower >= 0 ) then
      call MPI_ISEND( a( il(1):il(2) ), comm_size, this%mpi_dtype, pid_lower, 1, this%mpi_comm, sid1, ierr )
      call MPI_ISEND( b( il(1):il(2) ), comm_size, this%mpi_dtype, pid_lower, 2, this%mpi_comm, sid2, ierr )
      call MPI_ISEND( c( il(1):il(2) ), comm_size, this%mpi_dtype, pid_lower, 3, this%mpi_comm, sid3, ierr )
    endif
    ri0 = this%n_local + 1
    ri1 = this%n_local + comm_size
    if ( pid_upper <= this%num_procs - 1 ) then
      call MPI_RECV( a(ri0:ri1), comm_size, this%mpi_dtype, pid_upper, 1, this%mpi_comm, stat, ierr )
      call MPI_RECV( b(ri0:ri1), comm_size, this%mpi_dtype, pid_upper, 2, this%mpi_comm, stat, ierr )
      call MPI_RECV( c(ri0:ri1), comm_size, this%mpi_dtype, pid_upper, 3, this%mpi_comm, stat, ierr )
    else
      a(ri0:ri1) = 0.0
      b(ri0:ri1) = 1.0
      c(ri0:ri1) = 0.0
    endif
    if ( pid_lower >= 0 ) then
      call MPI_WAIT( sid1, stat, ierr )
      call MPI_WAIT( sid2, stat, ierr )
      call MPI_WAIT( sid3, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! send data to upper neighbor processor
    if ( pid_upper <= this%num_procs - 1 ) then
      call MPI_ISEND( a( iu(1):iu(2) ), comm_size, this%mpi_dtype, pid_upper, 1, this%mpi_comm, sid1, ierr )
      call MPI_ISEND( b( iu(1):iu(2) ), comm_size, this%mpi_dtype, pid_upper, 2, this%mpi_comm, sid2, ierr )
      call MPI_ISEND( c( iu(1):iu(2) ), comm_size, this%mpi_dtype, pid_upper, 3, this%mpi_comm, sid3, ierr )
    endif
    ri0 = 1 - comm_size
    ri1 = 0
    if ( pid_lower >= 0 ) then
      call MPI_RECV( a(ri0:ri1), comm_size, this%mpi_dtype, pid_lower, 1, this%mpi_comm, stat, ierr )
      call MPI_RECV( b(ri0:ri1), comm_size, this%mpi_dtype, pid_lower, 2, this%mpi_comm, stat, ierr )
      call MPI_RECV( c(ri0:ri1), comm_size, this%mpi_dtype, pid_lower, 3, this%mpi_comm, stat, ierr )
    else
      a(ri0:ri1) = 0.0
      b(ri0:ri1) = 1.0
      c(ri0:ri1) = 0.0
    endif
    if ( pid_upper <= this%num_procs - 1 ) then
      call MPI_WAIT( sid1, stat, ierr )
      call MPI_WAIT( sid2, stat, ierr )
      call MPI_WAIT( sid3, stat, ierr )
    endif
    call MPI_BARRIER( this%mpi_comm, ierr )

    ! calculate the CR coefficient
    do i = 1, this%n_local

      im = i - comm_size
      ip = i + comm_size
 
      this%coefm(i, step) = -a(i) / b(im)
      this%coefp(i, step) = -c(i) / b(ip)

      ! use class members as temporary arrays
      this%a(i) = this%coefm(i,step) * a(im)
      this%b(i) = this%coefm(i,step) * c(im) + this%coefp(i,step) * a(ip) + b(i)
      this%c(i) = this%coefp(i,step) * c(ip)

    enddo

    a(1:this%n_local) = this%a
    b(1:this%n_local) = this%b
    c(1:this%n_local) = this%c

    stride = stride * 2

  enddo

  this%b = b(1:this%n_local)
  deallocate( a, b, c )

  ! off-diagonal elements are all zeros, can be deleted
  ! deallocate( this%a, this%c )

end subroutine generate_cr_coef_fpcr_tri

subroutine set_values_rhs_fpcr_tri( this, rhs )

  implicit none
  class( fpcr_tri ), intent(inout) :: this
  real, intent(in), dimension(:) :: rhs
  integer :: i

  do i = 1, this%n_local
    this%rhs(i) = rhs(i)
  enddo

end subroutine set_values_rhs_fpcr_tri

subroutine set_values_matrix_fpcr_tri( this, a, b, c )

  implicit none
  class( fpcr_tri ), intent(inout) :: this
  real, intent(in), dimension(:) :: a, b, c
  integer :: i

  do i = 1, this%n_local
    this%a(i) = a(i)
    this%b(i) = b(i)
    this%c(i) = c(i)
  enddo
  if ( this%myid == 0 ) this%a(1) = 0.0
  if ( this%myid == this%num_procs - 1 ) this%c(this%n_local) = 0.0

end subroutine set_values_matrix_fpcr_tri

subroutine set_values_matrix_byrow_fpcr_tri( this, a, b, c, i_loc )

  implicit none
  class( fpcr_tri ), intent(inout) :: this
  real, intent(in) :: a, b, c
  integer, intent(in) :: i_loc
  integer :: i

  this%a(i_loc) = a
  this%b(i_loc) = b
  this%c(i_loc) = c
  if ( this%myid == 0 .and. i_loc == 1 ) this%a(i_loc) = 0.0
  if ( this%myid == this%num_procs - 1 .and. i_loc == this%n_local ) this%c(i_loc) = 0.0

end subroutine set_values_matrix_byrow_fpcr_tri

subroutine get_values_x_fpcr_tri( this, x )

  implicit none
  class( fpcr_tri ), intent(inout) :: this
  real, intent(out), dimension(:) :: x
  integer :: i

  do i = 1, this%n_local
    x(i) = this%x(i)
  enddo

end subroutine get_values_x_fpcr_tri

subroutine print_x_fpcr_tri( this, filename )

  implicit none
  class( fpcr_tri ), intent(inout) :: this
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

end subroutine print_x_fpcr_tri

subroutine print_rhs_fpcr_tri( this, filename )

  implicit none
  class( fpcr_tri ), intent(inout) :: this
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

end subroutine print_rhs_fpcr_tri

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

end module fpcr_tri_class