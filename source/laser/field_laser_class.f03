module field_laser_class

use parallel_module
use options_class
use input_class
use field_complex_class
use ufield_class
use param
use kwargs_class
use sysutil_module
use mpi
use fpcr_penta_class
use profile_laser_class

implicit none

private

character(len=32), parameter :: cls_name = "field_laser"
integer, parameter :: cls_level = 3

public :: field_laser

type, extends( field_complex ) :: field_laser

  private

  type( fpcr_penta ), dimension(:), pointer :: pgc_solver => null()
  real, dimension(:), pointer :: buf_re => null(), buf_im => null()

  class( ufield ), dimension(:), pointer, public :: sr_re => null()
  class( ufield ), dimension(:), pointer, public :: sr_im => null()
  class( ufield ), dimension(:), pointer, public :: si_re => null()
  class( ufield ), dimension(:), pointer, public :: si_im => null()

  class( ufield ), dimension(:), pointer :: axir_re => null()
  class( ufield ), dimension(:), pointer :: axir_im => null()
  class( ufield ), dimension(:), pointer :: axii_re => null()
  class( ufield ), dimension(:), pointer :: axii_im => null()

  ! central laser frequency
  real, public :: k0 = 10

  ! iteration times
  integer, public :: iter = 1

  ! intensity profile
  class( profile_laser ), allocatable :: profile

  contains

  procedure :: new   => init_field_laser
  procedure :: alloc => alloc_field_laser
  procedure :: del   => end_field_laser
  procedure :: solve => solve_field_laser
  procedure :: set_rhs => set_rhs_field_laser
  procedure, private :: init_solver

end type field_laser

contains

subroutine alloc_field_laser( this, input, opts, id )

  implicit none

  class( field_laser ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  integer, intent(in) :: id

  integer :: max_mode

  call input%get( 'simulation.max_mode', max_mode )
  if ( .not. associated( this%pgc_solver ) ) then
    allocate( this%pgc_solver(0:max_mode) )
  endif

  ! initialize profile
  allocate( this%profile )
  call this%profile%new( input, opts, id )

end subroutine alloc_field_laser

subroutine init_field_laser( this, opts, dim, max_mode, gc_num, only_f1, kwargs )

  implicit none
  class( field_laser ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: dim, max_mode
  integer, intent(in), dimension(2,2) :: gc_num
  logical, intent(in), optional :: only_f1
  type( kw_list ), intent(in), optional :: kwargs

  integer :: i, nrp, nr, noff
  real :: dr, dt, dz
  integer :: id, ierr
  integer, dimension(MPI_STATUS_SIZE) :: istat
  character(len=32), save :: sname = "init_field_laser"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = opts%get_ndp(1)
  nr  = opts%get_nd(1)
  noff= opts%get_noff(1)
  dr  = opts%get_dr()
  dz  = opts%get_dxi()
  dt  = opts%get_dt()

  call kwargs%get( 'k0', this%k0 )
  call kwargs%get( 'iter', this%iter )

  ! initialize the rhs
  allocate( this%sr_re(0:max_mode) )
  allocate( this%si_re(0:max_mode) )
  allocate( this%sr_im(max_mode) )
  allocate( this%si_im(max_mode) )
  allocate( this%axir_re(0:max_mode) )
  allocate( this%axii_re(0:max_mode) )
  allocate( this%axir_im(max_mode) )
  allocate( this%axii_im(max_mode) )
  do i = 0, max_mode
    call this%sr_re(i)%new( opts, dim, i, gc_num, has_2d=.true. )
    call this%si_re(i)%new( opts, dim, i, gc_num, has_2d=.true. )
    call this%axir_re(i)%new( opts, dim, i, gc_num, has_2d=.false. )
    call this%axii_re(i)%new( opts, dim, i, gc_num, has_2d=.false. )
    if (i==0) cycle
    call this%sr_im(i)%new( opts, dim, i, gc_num, has_2d=.true. )
    call this%si_im(i)%new( opts, dim, i, gc_num, has_2d=.true. )
    call this%axir_im(i)%new( opts, dim, i, gc_num, has_2d=.false. )
    call this%axii_im(i)%new( opts, dim, i, gc_num, has_2d=.false. )
  enddo

  ! call initialization routine of the parent class
  call this%field_complex%new( opts, dim, max_mode, gc_num )

  ! initialize solver
  call this%init_solver( nr, nrp, noff, this%k0, dt, dr, dz )

  ! launch laser  
  call this%profile%launch( this%cfr_re, this%cfr_im, this%cfi_re, this%cfi_im )
  call this%copy_gc_f2()
  call this%pipe_send_f2( 1, id, 'forward', 'guard' )
  call mpi_wait( id, istat, ierr )
  call this%pipe_send_f2( 1, id, 'backward', 'guard' )
  call mpi_wait( id, istat, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_laser

subroutine end_field_laser( this )

  implicit none
  class( field_laser ), intent(inout) :: this

  integer :: i
  character(len=32), save :: sname = 'end_field_laser'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%pgc_solver(i)%destroy()
    call this%sr_re(i)%del()
    call this%si_re(i)%del()
    call this%axir_re(i)%del()
    call this%axii_re(i)%del()
    if ( i == 0 ) cycle
    call this%sr_im(i)%del()
    call this%si_im(i)%del()
    call this%axir_im(i)%del()
    call this%axii_im(i)%del()
  enddo
  deallocate( this%pgc_solver )
  deallocate( this%sr_re, this%sr_im, this%si_re, this%si_im )
  deallocate( this%axir_re, this%axir_im, this%axii_re, this%axii_im )
  if ( associated( this%buf_re ) ) deallocate( this%buf_re )
  if ( associated( this%buf_im ) ) deallocate( this%buf_im )
  call this%field_complex%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_laser

subroutine init_solver( this, nr, nrp, noff, k0, ds, dr, dz )

  implicit none
  class( field_laser ), intent(inout) :: this
  integer, intent(in) :: nr, nrp, noff
  real, intent(in) :: k0, ds, dr, dz

  integer :: ierr, m, i, j, local_size, idproc, nvp
  real :: a, b, c, d, e, m2, j2

  local_size = 2 * nrp
  nvp = num_procs_loc()
  idproc = id_proc_loc()

  do m = 0, this%max_mode
    call this%pgc_solver(m)%create( 2*nr, comm_loc(), ierr )

    m2 = m * m
    ! set matrix element
    do i = 1, local_size, 2

      j = (noff + i + 1) / 2 - 1
      j2 = j * j

      a = -0.25 * ds * ( 1.0 - 0.5 / j )
      b = 0.0
      c = 0.25 * ds * ( 2.0 + m2 / j2 )
      d = -k0 * dr**2
      e = -0.25 * ds * ( 1.0 + 0.5 / j )
      call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, i )

      a = -0.25 * ds * ( 1.0 - 0.5 / j )
      b = k0 * dr**2
      c = 0.25 * ds * ( 2.0 + m2 / j2 )
      d = 0.0
      e = -0.25 * ds * ( 1.0 + 0.5 / j )
      call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, i+1 )

    enddo

    ! set the axial boundary condition
    if (idproc == 0) then

      if ( m == 0 ) then

        a = 0.0
        b = epsilon(1.0)
        c = ds
        d = -k0 * dr**2
        e = -ds
        call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, 1 )

        a = epsilon(1.0)
        b = k0 * dr**2
        c = ds
        d = 0.0
        e = -ds
        call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, 2 )

      else

        ! matrix elements of row 1 and 2 are given arbitrarily to make sure the matrix
        ! is not singular.
        a = 0.0
        b = epsilon(1.0)
        c = 1.0
        d = 0.0
        e = 0.0
        call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, 1 )

        a = epsilon(1.0)
        b = 0.0
        c = 1.0
        d = 0.0
        e = 0.0
        call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, 2 )

        ! first matrix elements of row 3 and 4 are given zeros indicating the
        ! on-axis values are zeros
        a = 0.0
        b = 0.0
        c = 0.25 * ds * ( 2.0 + m2 )
        d = -k0 * dr**2
        e = -0.375 * ds
        call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, 3 )

        a = 0.0
        b = k0 * dr**2
        c = 0.25 * ds * ( 2.0 + m2 )
        d = 0.0
        e = -0.375 * ds
        call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, 4 )

      endif

    endif

    ! set the outter boundary condition
    if (idproc == nvp - 1) then

      j = nr - 1
      j2 = j * j

      a = -0.25 * ds * ( 1.0 - 0.5 / j )
      b = 0.0
      c = 0.25 * ds * ( 2.0 + m2 / j2 )
      d = -k0 * dr**2
      e = 0.0
      call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, local_size-1 )

      a = -0.25 * ds * ( 1.0 - 0.5 / j )
      b = k0 * dr**2
      c = 0.25 * ds * ( 2.0 + m2 / j2 )
      d = 0.0
      e = 0.0
      call this%pgc_solver(m)%set_values_matrix( a, b, c, d, e, local_size )

    endif

    ! generate cyclic reduction coefficients
    call this%pgc_solver(m)%generate_cr_coef()

  enddo

end subroutine init_solver

! subroutine set_rhs_field_laser( this, chi )
subroutine set_rhs_field_laser( this )

  implicit none
  class( field_laser ), intent(inout) :: this

  integer :: m, i, j, nrp, nzp
  real :: a1, b1, c1, d1, e1
  real :: a2, b2, c2, d2, e2
  real :: dr2_idzh
  real, dimension(:,:,:), pointer :: ar_re => null(), ar_im => null(), ai_re => null(), ai_im => null()

  nrp    = this%cfr_re(0)%get_ndp(1)
  nzp    = this%cfr_re(0)%get_ndp(2)
  dr2_idzh = 0.5 * this%dr**2 / this%dz
  
  ! mode = 0
  ar_re => this%cfr_re(0)%get_f2()
  ai_re => this%cfi_re(0)%get_f2()
  do j = 2 - this%gc_num(1,2), nzp + this%gc_num(2,2) - 1
    do i = 1, nrp
      call this%pgc_solver(0)%get_values_matrix( a1, b1, c1, d1, e1, 2*i-1 )
      call this%pgc_solver(0)%get_values_matrix( a2, b2, c2, d2, e2, 2*i )

      this%sr_re(0)%f2(1,i,j) = dr2_idzh * ( ar_re(1,i,j+1) - ar_re(1,i,j-1) )
      this%si_re(0)%f2(1,i,j) = dr2_idzh * ( ai_re(1,i,j+1) - ai_re(1,i,j-1) )
      
      this%sr_re(0)%f2(1,i,j) = this%sr_re(0)%f2(1,i,j) + d1 * ai_re(1,i,j) &
                              - a1 * ar_re(1,i-1,j) - c1 * ar_re(1,i,j) - e1 * ar_re(1,i+1,j)
      this%si_re(0)%f2(1,i,j) = this%si_re(0)%f2(1,i,j) + b2 * ar_re(1,i,j) &
                              - a2 * ai_re(1,i-1,j) - c2 * ai_re(1,i,j) - e2 * ai_re(1,i+1,j)                              
    enddo
  enddo

  ! mode > 0
  do m = 1, this%max_mode

    ar_re => this%cfr_re(m)%get_f2()
    ai_re => this%cfi_re(m)%get_f2()
    ar_im => this%cfr_im(m)%get_f2()
    ai_im => this%cfi_im(m)%get_f2()

    do j = 2 - this%gc_num(1,2), nzp + this%gc_num(2,2) - 1
      do i = 1, nrp
        call this%pgc_solver(m)%get_values_matrix( a1, b1, c1, d1, e1, 2*i-1 )
        call this%pgc_solver(m)%get_values_matrix( a2, b2, c2, d2, e2, 2*i )

        this%sr_re(m)%f2(1,i,j) = dr2_idzh * ( ar_re(1,i,j+1) - ar_re(1,i,j-1) )
        this%sr_im(m)%f2(1,i,j) = dr2_idzh * ( ar_im(1,i,j+1) - ar_im(1,i,j-1) )
        this%si_re(m)%f2(1,i,j) = dr2_idzh * ( ai_re(1,i,j+1) - ai_re(1,i,j-1) )
        this%si_im(m)%f2(1,i,j) = dr2_idzh * ( ai_im(1,i,j+1) - ai_im(1,i,j-1) )

        this%sr_re(m)%f2(1,i,j) = this%sr_re(m)%f2(1,i,j) + d1 * ai_re(1,i,j) &
                                - a1 * ar_re(1,i-1,j) - c1 * ar_re(1,i,j) - e1 * ar_re(1,i+1,j)
        this%sr_im(m)%f2(1,i,j) = this%sr_im(m)%f2(1,i,j) + d1 * ai_im(1,i,j) &
                                - a1 * ar_im(1,i-1,j) - c1 * ar_im(1,i,j) - e1 * ar_im(1,i+1,j)
        this%si_re(m)%f2(1,i,j) = this%si_re(m)%f2(1,i,j) + b2 * ar_re(1,i,j) &
                                - a2 * ai_re(1,i-1,j) - c2 * ai_re(1,i,j) - e2 * ai_re(1,i+1,j)
        this%si_im(m)%f2(1,i,j) = this%si_im(m)%f2(1,i,j) + b2 * ar_im(1,i,j) &
                                - a2 * ai_im(1,i-1,j) - c2 * ai_im(1,i,j) - e2 * ai_im(1,i+1,j)
      enddo
    enddo
  enddo

end subroutine set_rhs_field_laser

! subroutine solve_field_laser( this, chi, i_slice )
subroutine solve_field_laser( this, i_slice )

  implicit none
  class( field_laser ), intent(inout) :: this
  integer, intent(in) :: i_slice

  integer :: m, i, j, nrp, nzp
  real :: dr2_idzh, rhs

  nrp    = this%cfr_re(0)%get_ndp(1)
  nzp    = this%cfr_re(0)%get_ndp(2)
  dr2_idzh = 0.5 * this%dr**2 / this%dz

  ! handle the first slice, calculate the da/dxi for the first slice
  if ( i_slice == 2 - this%gc_num(1,2) ) then
    do m = 0, this%max_mode
      do i = 1, nrp
        this%axir_re(m)%f1(1,i) = dr2_idzh * ( this%cfr_re(m)%f2(1,i,i_slice+1) - this%cfr_re(m)%f2(1,i,i_slice-1) )
        this%axii_re(m)%f1(1,i) = dr2_idzh * ( this%cfi_re(m)%f2(1,i,i_slice+1) - this%cfi_re(m)%f2(1,i,i_slice-1) )
      enddo
      if ( m == 0 ) cycle
      do i = 1, nrp
        this%axir_im(m)%f1(1,i) = dr2_idzh * ( this%cfr_im(m)%f2(1,i,i_slice+1) - this%cfr_im(m)%f2(1,i,i_slice-1) )
        this%axii_im(m)%f1(1,i) = dr2_idzh * ( this%cfi_im(m)%f2(1,i,i_slice+1) - this%cfi_im(m)%f2(1,i,i_slice-1) )
      enddo
    enddo
  endif

  ! copy slices
  ! there is no need to copy slice.
  ! call this%copy_slice( i_slice, p_copy_2to1 )
  ! do i = 0, this%max_mode
  !   call this%sr_re(i)%copy_slice( idx, p_copy_2to1 )
  !   call this%si_re(i)%copy_slice( idx, p_copy_2to1 )
  !   if ( i == 0 ) cycle
  !   call this%sr_im(i)%copy_slice( idx, p_copy_2to1 )
  !   call this%si_im(i)%copy_slice( idx, p_copy_2to1 )
  ! enddo

  ! set rhs of the PCR solver and solve
  do m = 0, this%max_mode

    ! set rhs of PCR
    do i = 1, nrp
      rhs = this%sr_re(m)%f2(1,i,i_slice) - this%axir_re(m)%f1(1,i)
      call this%pgc_solver(m)%set_values_rhs( rhs, 2*i-1 )
      rhs = this%si_re(m)%f2(1,i,i_slice) - this%axii_re(m)%f1(1,i)
      call this%pgc_solver(m)%set_values_rhs( rhs, 2*i )
    enddo

    call this%pgc_solver(m)%solve()

    ! calculate the da/dxi at next slice before a is replaced
    if ( i_slice < nzp + this%gc_num(2,2) - 1 ) then
      do i = 1, nrp
        this%axir_re(m)%f1(1,i) = dr2_idzh * ( this%cfr_re(m)%f2(1,i,i_slice+2) - this%cfr_re(m)%f2(1,i,i_slice) )
        this%axii_re(m)%f1(1,i) = dr2_idzh * ( this%cfi_re(m)%f2(1,i,i_slice+2) - this%cfi_re(m)%f2(1,i,i_slice) )
      enddo
    else
      this%axir_re(m)%f1(1,i) = 0.0
      this%axii_re(m)%f1(1,i) = 0.0
    endif

    ! get solution
    do i = 1, nrp
      call this%pgc_solver(m)%get_values_x( this%cfr_re(m)%f2(1,i,i_slice), 2*i-1 )
      call this%pgc_solver(m)%get_values_x( this%cfi_re(m)%f2(1,i,i_slice), 2*i )
    enddo

    if ( m == 0 ) cycle

    ! set rhs of PCR
    do i = 1, nrp
      rhs = this%sr_im(m)%f2(1,i,i_slice) - this%axir_im(m)%f1(1,i)
      call this%pgc_solver(m)%set_values_rhs( rhs, 2*i-1 )
      rhs = this%si_im(m)%f2(1,i,i_slice) - this%axii_im(m)%f1(1,i)
      call this%pgc_solver(m)%set_values_rhs( rhs, 2*i )
    enddo

    call this%pgc_solver(m)%solve()

    ! calculate the da/dxi at next slice before a is replaced
    if ( i_slice < nzp + this%gc_num(2,2) - 1 ) then
      do i = 1, nrp
        this%axir_im(m)%f1(1,i) = dr2_idzh * ( this%cfr_im(m)%f2(1,i,i_slice+2) - this%cfr_im(m)%f2(1,i,i_slice) )
        this%axii_im(m)%f1(1,i) = dr2_idzh * ( this%cfi_im(m)%f2(1,i,i_slice+2) - this%cfi_im(m)%f2(1,i,i_slice) )
      enddo
    else
      this%axir_im(m)%f1(1,i) = 0.0
      this%axii_im(m)%f1(1,i) = 0.0
    endif

    ! get solution
    do i = 1, nrp
      call this%pgc_solver(m)%get_values_x( this%cfr_im(m)%f2(1,i,i_slice), 2*i-1 )
      call this%pgc_solver(m)%get_values_x( this%cfi_im(m)%f2(1,i,i_slice), 2*i )
    enddo

  enddo
    
  call this%copy_gc_f2()

end subroutine solve_field_laser

end module field_laser_class