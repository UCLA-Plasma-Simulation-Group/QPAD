module field_laser_class

use parallel_module
use options_class
use input_class
use field_class
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

  class( ufield ), dimension(:), pointer, public :: sr_re => null()
  class( ufield ), dimension(:), pointer, public :: sr_im => null()
  class( ufield ), dimension(:), pointer, public :: si_re => null()
  class( ufield ), dimension(:), pointer, public :: si_im => null()

  class( ufield ), dimension(:), pointer :: axir_re => null()
  class( ufield ), dimension(:), pointer :: axir_im => null()
  class( ufield ), dimension(:), pointer :: axii_re => null()
  class( ufield ), dimension(:), pointer :: axii_im => null()

  class( ufield ), dimension(:), pointer, public :: ar_grad_re => null()
  class( ufield ), dimension(:), pointer, public :: ar_grad_im => null()
  class( ufield ), dimension(:), pointer, public :: ai_grad_re => null()
  class( ufield ), dimension(:), pointer, public :: ai_grad_im => null()

  ! central laser frequency
  real, public :: k0 = 10

  ! 3D time step
  real :: ds

  ! iteration times
  integer, public :: iter = 1

  ! intensity profile
  class( profile_laser ), allocatable :: profile

  contains

  procedure :: new      => init_field_laser
  procedure :: new_aux  => init_auxiliary_field_laser
  procedure :: alloc    => alloc_field_laser
  procedure :: del      => end_field_laser
  procedure :: del_aux  => end_auxiliary_field_laser
  procedure :: solve    => solve_field_laser
  procedure :: set_rhs  => set_rhs_field_laser
  procedure :: set_grad => set_grad_field_laser
  procedure :: gather   => gather_field_laser
  procedure :: zero     => zero_field_laser
  procedure, private :: init_solver

end type field_laser

real, dimension(:), allocatable :: tmp_r_re, tmp_r_im, tmp_i_re, tmp_i_im

contains

subroutine alloc_field_laser( this, input, opts, id )

  implicit none

  class( field_laser ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  integer, intent(in) :: id

  integer :: max_mode
  character(len=32), save :: sname = 'alloc_field_laser'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.max_mode', max_mode )
  if ( .not. associated( this%pgc_solver ) ) then
    allocate( this%pgc_solver(0:max_mode) )
  endif

  ! initialize profile
  allocate( this%profile )
  call this%profile%new( input, opts, id )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

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
  real :: dr, dz
  integer :: id, ierr
  integer, dimension(MPI_STATUS_SIZE) :: istat
  character(len=32), save :: sname = "init_field_laser"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = opts%get_ndp(1)
  nr  = opts%get_nd(1)
  noff= opts%get_noff(1)
  dr  = opts%get_dr()
  dz  = opts%get_dxi()
  this%ds = opts%get_dt()

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
  allocate( this%ar_grad_re(0:max_mode) )
  allocate( this%ai_grad_re(0:max_mode) )
  allocate( this%ar_grad_im(max_mode) )
  allocate( this%ai_grad_im(max_mode) )
  do i = 0, max_mode
    call this%sr_re(i)%new( opts, 1, i, gc_num, has_2d=.true. )
    call this%si_re(i)%new( opts, 1, i, gc_num, has_2d=.true. )
    call this%axir_re(i)%new( opts, 1, i, gc_num, has_2d=.false. )
    call this%axii_re(i)%new( opts, 1, i, gc_num, has_2d=.false. )
    call this%ar_grad_re(i)%new( opts, 2, i, gc_num, has_2d=.false. )
    call this%ai_grad_re(i)%new( opts, 2, i, gc_num, has_2d=.false. )
    if (i==0) cycle
    call this%sr_im(i)%new( opts, 1, i, gc_num, has_2d=.true. )
    call this%si_im(i)%new( opts, 1, i, gc_num, has_2d=.true. )
    call this%axir_im(i)%new( opts, 1, i, gc_num, has_2d=.false. )
    call this%axii_im(i)%new( opts, 1, i, gc_num, has_2d=.false. )
    call this%ar_grad_im(i)%new( opts, 2, i, gc_num, has_2d=.false. )
    call this%ai_grad_im(i)%new( opts, 2, i, gc_num, has_2d=.false. )
  enddo

  ! call initialization routine of the parent class
  call this%field_complex%new( opts, dim, max_mode, gc_num, only_f1 )

  ! initialize solver
  call this%init_solver( nr, nrp, noff, this%k0, this%ds, dr, dz )

  ! launch laser  
  call this%profile%launch( this%cfr_re, this%cfr_im, this%cfi_re, this%cfi_im )
  call this%copy_gc_f2()
  call this%pipe_send( 1, id, 'forward', 'inner' )
  call this%pipe_recv( 1, 'forward', 'guard', 'replace' )
  call mpi_wait( id, istat, ierr )
  call this%pipe_send( 1, id, 'backward', 'inner' )
  call this%pipe_recv( 1, 'backward', 'guard', 'replace' )
  call mpi_wait( id, istat, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_laser

subroutine init_auxiliary_field_laser( this, opts, dim, max_mode, gc_num, only_f1, kwargs )

  implicit none
  class( field_laser ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: dim, max_mode
  integer, intent(in), dimension(2,2) :: gc_num
  logical, intent(in), optional :: only_f1
  type( kw_list ), intent(in), optional :: kwargs

  integer :: i, nrp, nr, noff
  real :: dr, dz
  integer :: id, ierr
  integer, dimension(MPI_STATUS_SIZE) :: istat
  character(len=32), save :: sname = "init_auxiliary_field_laser"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%ds = opts%get_dt()

  ! initialize the rhs
  allocate( this%ar_grad_re(0:max_mode) )
  allocate( this%ai_grad_re(0:max_mode) )
  allocate( this%ar_grad_im(max_mode) )
  allocate( this%ai_grad_im(max_mode) )
  do i = 0, max_mode
    call this%ar_grad_re(i)%new( opts, 2, i, gc_num, has_2d=.false. )
    call this%ai_grad_re(i)%new( opts, 2, i, gc_num, has_2d=.false. )
    if (i==0) cycle
    call this%ar_grad_im(i)%new( opts, 2, i, gc_num, has_2d=.false. )
    call this%ai_grad_im(i)%new( opts, 2, i, gc_num, has_2d=.false. )
  enddo

  ! call initialization routine of the parent class
  call this%field_complex%new( opts, dim, max_mode, gc_num, only_f1 )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_auxiliary_field_laser

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
    call this%ar_grad_re(i)%del()
    call this%ai_grad_re(i)%del()
    if ( i == 0 ) cycle
    call this%sr_im(i)%del()
    call this%si_im(i)%del()
    call this%axir_im(i)%del()
    call this%axii_im(i)%del()
    call this%ar_grad_im(i)%del()
    call this%ai_grad_im(i)%del()
  enddo
  deallocate( this%pgc_solver )
  deallocate( this%sr_re, this%sr_im, this%si_re, this%si_im )
  deallocate( this%axir_re, this%axir_im, this%axii_re, this%axii_im )
  deallocate( this%ar_grad_re, this%ar_grad_im, this%ai_grad_re, this%ai_grad_im )
  call this%field_complex%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_laser

subroutine end_auxiliary_field_laser( this )

  implicit none
  class( field_laser ), intent(inout) :: this

  integer :: i
  character(len=32), save :: sname = 'end_auxiliary_field_laser'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%ar_grad_re(i)%del()
    call this%ai_grad_re(i)%del()
    if ( i == 0 ) cycle
    call this%ar_grad_im(i)%del()
    call this%ai_grad_im(i)%del()
  enddo
  deallocate( this%ar_grad_re, this%ar_grad_im, this%ai_grad_re, this%ai_grad_im )
  call this%field_complex%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_auxiliary_field_laser

subroutine init_solver( this, nr, nrp, noff, k0, ds, dr, dz )

  implicit none
  class( field_laser ), intent(inout) :: this
  integer, intent(in) :: nr, nrp, noff
  real, intent(in) :: k0, ds, dr, dz

  integer :: ierr, m, i, j, local_size, idproc, nvp
  real :: a, b, c, d, e, m2, j2
  character(len=32), save :: sname = 'init_solver'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  local_size = 2 * nrp
  nvp = num_procs_loc()
  idproc = id_proc_loc()

  do m = 0, this%max_mode
    call this%pgc_solver(m)%create( 2*nr, comm_loc(), ierr )

    m2 = m * m
    ! set matrix element
    do i = 1, local_size, 2

      j = (i + 1) / 2 - 1 + noff
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

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_solver

subroutine set_rhs_field_laser( this, chi )

  implicit none
  class( field_laser ), intent(inout) :: this
  class( field ), intent(in) :: chi

  integer :: m, i, j, k, nrp, nzp, noff
  integer :: nvp, idproc
  real :: beta_m, beta_p, alpha, kappa
  real :: dr2_idzh, ds_qtr, m2, ik
  real, dimension(:,:,:), pointer :: ar_re => null(), ar_im => null(), ai_re => null(), ai_im => null()
  real, dimension(:,:,:), pointer :: chi_re => null(), chi_im => null()
  character(len=32), save :: sname = 'set_rhs_field_laser'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp    = this%cfr_re(0)%get_ndp(1)
  nzp    = this%cfr_re(0)%get_ndp(2)
  noff   = this%cfr_re(0)%get_noff(1)
  dr2_idzh = 0.5 * this%dr**2 / this%dz
  kappa  = this%k0 * this%dr**2
  ds_qtr = 0.25 * this%ds
  idproc = id_proc_loc()
  nvp    = num_procs_loc()
  
  ! mode = 0
  ar_re => this%cfr_re(0)%get_f2()
  ai_re => this%cfi_re(0)%get_f2()
  do j = 2 - this%gc_num(1,2), nzp + this%gc_num(2,2) - 1

    ! calculate the inner cells
    do i = 2, nrp - 1
      ik = 1.0 / (i + noff - 1)
      beta_m = ds_qtr * ( 1.0 - 0.5 * ik )
      beta_p = ds_qtr * ( 1.0 + 0.5 * ik )
      alpha = -ds_qtr * 2.0
      this%sr_re(0)%f2(1,i,j) = dr2_idzh * ( ar_re(1,i,j+1) - ar_re(1,i,j-1) ) - kappa * ai_re(1,i,j) &
                              + beta_m * ar_re(1,i-1,j) + alpha * ar_re(1,i,j) + beta_p * ar_re(1,i+1,j)
      this%si_re(0)%f2(1,i,j) = dr2_idzh * ( ai_re(1,i,j+1) - ai_re(1,i,j-1) ) + kappa * ar_re(1,i,j) &
                              + beta_m * ai_re(1,i-1,j) + alpha * ai_re(1,i,j) + beta_p * ai_re(1,i+1,j)                              
    enddo

    ! calculate the first cell
    if ( idproc == 0 ) then
      i = 1
      beta_m = 0.0
      beta_p = this%ds
      alpha = -this%ds
    else
      i = 1
      ik = 1.0 / noff
      beta_m = ds_qtr * ( 1.0 - 0.5 * ik )
      beta_p = ds_qtr * ( 1.0 + 0.5 * ik )
      alpha = -ds_qtr * 2.0
    endif
    this%sr_re(0)%f2(1,i,j) = dr2_idzh * ( ar_re(1,i,j+1) - ar_re(1,i,j-1) ) - kappa * ai_re(1,i,j) &
                              + beta_m * ar_re(1,i-1,j) + alpha * ar_re(1,i,j) + beta_p * ar_re(1,i+1,j)
    this%si_re(0)%f2(1,i,j) = dr2_idzh * ( ai_re(1,i,j+1) - ai_re(1,i,j-1) ) + kappa * ar_re(1,i,j) &
                              + beta_m * ai_re(1,i-1,j) + alpha * ai_re(1,i,j) + beta_p * ai_re(1,i+1,j)

    ! calculate the last cell
    if ( idproc == nvp - 1 ) then
      i = nrp
      ik = 1.0 / (i + noff - 1)
      beta_m = ds_qtr * ( 1.0 - 0.5 * ik )
      beta_p = 0.0
      alpha = -ds_qtr * 2.0
    else
      i = nrp
      ik = 1.0 / (i + noff - 1)
      beta_m = ds_qtr * ( 1.0 - 0.5 * ik )
      beta_p = ds_qtr * ( 1.0 + 0.5 * ik )
      alpha = -ds_qtr * 2.0
    endif
    this%sr_re(0)%f2(1,i,j) = dr2_idzh * ( ar_re(1,i,j+1) - ar_re(1,i,j-1) ) - kappa * ai_re(1,i,j) &
                              + beta_m * ar_re(1,i-1,j) + alpha * ar_re(1,i,j) + beta_p * ar_re(1,i+1,j)
    this%si_re(0)%f2(1,i,j) = dr2_idzh * ( ai_re(1,i,j+1) - ai_re(1,i,j-1) ) + kappa * ar_re(1,i,j) &
                              + beta_m * ai_re(1,i-1,j) + alpha * ai_re(1,i,j) + beta_p * ai_re(1,i+1,j)
  enddo

  ! mode > 0
  do m = 1, this%max_mode

    ar_re => this%cfr_re(m)%get_f2()
    ai_re => this%cfi_re(m)%get_f2()
    ar_im => this%cfr_im(m)%get_f2()
    ai_im => this%cfi_im(m)%get_f2()
    m2 = m * m

    ! calculate the inner cells
    do j = 2 - this%gc_num(1,2), nzp + this%gc_num(2,2) - 1
      do i = 2, nrp - 1
        ik = 1.0 / (i + noff - 1)
        beta_m = ds_qtr * ( 1.0 - 0.5 * ik )
        beta_p = ds_qtr * ( 1.0 + 0.5 * ik )
        alpha = -ds_qtr * ( 2.0 + m2 * ik * ik )

        this%sr_re(m)%f2(1,i,j) = dr2_idzh * ( ar_re(1,i,j+1) - ar_re(1,i,j-1) ) - kappa * ai_re(1,i,j) &
                                + beta_m * ar_re(1,i-1,j) + alpha * ar_re(1,i,j) + beta_p * ar_re(1,i+1,j)
        this%sr_im(m)%f2(1,i,j) = dr2_idzh * ( ar_im(1,i,j+1) - ar_im(1,i,j-1) ) - kappa * ai_im(1,i,j) &
                                + beta_m * ar_im(1,i-1,j) + alpha * ar_im(1,i,j) + beta_p * ar_im(1,i+1,j)
        this%si_re(m)%f2(1,i,j) = dr2_idzh * ( ai_re(1,i,j+1) - ai_re(1,i,j-1) ) + kappa * ar_re(1,i,j) &
                                + beta_m * ai_re(1,i-1,j) + alpha * ai_re(1,i,j) + beta_p * ai_re(1,i+1,j)
        this%si_im(m)%f2(1,i,j) = dr2_idzh * ( ai_im(1,i,j+1) - ai_im(1,i,j-1) ) + kappa * ar_im(1,i,j) &
                                + beta_m * ai_im(1,i-1,j) + alpha * ai_im(1,i,j) + beta_p * ai_im(1,i+1,j)
      enddo

      ! calculate the first cell
      if ( idproc == 0 ) then
        this%sr_re(m)%f2(1,i,j) = 0.0
        this%sr_im(m)%f2(1,i,j) = 0.0
        this%si_re(m)%f2(1,i,j) = 0.0
        this%si_im(m)%f2(1,i,j) = 0.0
      else
        i = 1
        ik = 1.0 / noff
        beta_m = ds_qtr * ( 1.0 - 0.5 * ik )
        beta_p = ds_qtr * ( 1.0 + 0.5 * ik )
        alpha = -ds_qtr * ( 2.0 + m2 * ik * ik )

        this%sr_re(m)%f2(1,i,j) = dr2_idzh * ( ar_re(1,i,j+1) - ar_re(1,i,j-1) ) - kappa * ai_re(1,i,j) &
                                + beta_m * ar_re(1,i-1,j) + alpha * ar_re(1,i,j) + beta_p * ar_re(1,i+1,j)
        this%sr_im(m)%f2(1,i,j) = dr2_idzh * ( ar_im(1,i,j+1) - ar_im(1,i,j-1) ) - kappa * ai_im(1,i,j) &
                                + beta_m * ar_im(1,i-1,j) + alpha * ar_im(1,i,j) + beta_p * ar_im(1,i+1,j)
        this%si_re(m)%f2(1,i,j) = dr2_idzh * ( ai_re(1,i,j+1) - ai_re(1,i,j-1) ) + kappa * ar_re(1,i,j) &
                                + beta_m * ai_re(1,i-1,j) + alpha * ai_re(1,i,j) + beta_p * ai_re(1,i+1,j)
        this%si_im(m)%f2(1,i,j) = dr2_idzh * ( ai_im(1,i,j+1) - ai_im(1,i,j-1) ) + kappa * ar_im(1,i,j) &
                                + beta_m * ai_im(1,i-1,j) + alpha * ai_im(1,i,j) + beta_p * ai_im(1,i+1,j)
      endif

      ! calculate the last cell
      if ( idproc == nvp - 1 ) then
        i = nrp
        ik = 1.0 / (noff + i - 1)
        beta_m = ds_qtr * ( 1.0 - 0.5 * ik )
        beta_p = 0.0
        alpha = -ds_qtr * ( 2.0 + m2 * ik * ik )
      else
        i = nrp
        ik = 1.0 / (noff + i - 1)
        beta_m = ds_qtr * ( 1.0 - 0.5 * ik )
        beta_p = ds_qtr * ( 1.0 + 0.5 * ik )
        alpha = -ds_qtr * ( 2.0 + m2 * ik * ik )
      endif
      this%sr_re(m)%f2(1,i,j) = dr2_idzh * ( ar_re(1,i,j+1) - ar_re(1,i,j-1) ) - kappa * ai_re(1,i,j) &
                              + beta_m * ar_re(1,i-1,j) + alpha * ar_re(1,i,j) + beta_p * ar_re(1,i+1,j)
      this%sr_im(m)%f2(1,i,j) = dr2_idzh * ( ar_im(1,i,j+1) - ar_im(1,i,j-1) ) - kappa * ai_im(1,i,j) &
                              + beta_m * ar_im(1,i-1,j) + alpha * ar_im(1,i,j) + beta_p * ar_im(1,i+1,j)
      this%si_re(m)%f2(1,i,j) = dr2_idzh * ( ai_re(1,i,j+1) - ai_re(1,i,j-1) ) + kappa * ar_re(1,i,j) &
                              + beta_m * ai_re(1,i-1,j) + alpha * ai_re(1,i,j) + beta_p * ai_re(1,i+1,j)
      this%si_im(m)%f2(1,i,j) = dr2_idzh * ( ai_im(1,i,j+1) - ai_im(1,i,j-1) ) + kappa * ar_im(1,i,j) &
                              + beta_m * ai_im(1,i-1,j) + alpha * ai_im(1,i,j) + beta_p * ai_im(1,i+1,j)
    enddo
  enddo

  ! calculate the contribution from the plasma susceptibility
  do m = 0, this%max_mode
    
    do k = 0, m
      chi_re => chi%rf_re(k)%get_f2()
      ar_re  => this%cfr_re(m-k)%get_f2()
      ai_re  => this%cfi_re(m-k)%get_f2()

      do j = 2 - this%gc_num(1,2), nzp + this%gc_num(2,2) - 1
        do i = 1, nrp
          this%sr_re(m)%f2(1,i,j) = this%sr_re(m)%f2(1,i,j) + ds_qtr * chi_re(1,i,j) * ar_re(1,i,j)
          this%si_re(m)%f2(1,i,j) = this%si_re(m)%f2(1,i,j) + ds_qtr * chi_re(1,i,j) * ai_re(1,i,j)
        enddo
      enddo

      if ( k == 0 .or. k == m ) cycle

      chi_im => chi%rf_im(k)%get_f2()
      ar_im  => this%cfr_im(m-k)%get_f2()
      ai_im  => this%cfi_im(m-k)%get_f2()

      do j = 2 - this%gc_num(1,2), nzp + this%gc_num(2,2) - 1
        do i = 1, nrp
          this%sr_re(m)%f2(1,i,j) = this%sr_re(m)%f2(1,i,j) - ds_qtr * chi_im(1,i,j) * ar_im(1,i,j)
          this%si_re(m)%f2(1,i,j) = this%si_re(m)%f2(1,i,j) - ds_qtr * chi_im(1,i,j) * ai_im(1,i,j)
        enddo
      enddo
    enddo

    if ( m == 0 ) cycle

    do k = 0, m

      if ( k /= 0 ) then
        chi_im => chi%rf_im(k)%get_f2()
        ar_re  => this%cfr_re(m-k)%get_f2()
        ai_re  => this%cfi_re(m-k)%get_f2()

        do j = 2 - this%gc_num(1,2), nzp + this%gc_num(2,2) - 1
          do i = 1, nrp
            this%sr_im(m)%f2(1,i,j) = this%sr_im(m)%f2(1,i,j) + ds_qtr * chi_im(1,i,j) * ar_re(1,i,j)
            this%si_im(m)%f2(1,i,j) = this%si_im(m)%f2(1,i,j) + ds_qtr * chi_im(1,i,j) * ai_re(1,i,j)
          enddo
        enddo
      endif

      if ( k /= m ) then
        chi_re => chi%rf_re(k)%get_f2()
        ar_im  => this%cfr_im(m-k)%get_f2()
        ai_im  => this%cfi_im(m-k)%get_f2()

        do j = 2 - this%gc_num(1,2), nzp + this%gc_num(2,2) - 1
          do i = 1, nrp
            this%sr_im(m)%f2(1,i,j) = this%sr_im(m)%f2(1,i,j) + ds_qtr * chi_re(1,i,j) * ar_im(1,i,j)
            this%si_im(m)%f2(1,i,j) = this%si_im(m)%f2(1,i,j) + ds_qtr * chi_re(1,i,j) * ai_im(1,i,j)
          enddo
        enddo
      endif

    enddo

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_rhs_field_laser

subroutine set_grad_field_laser( this )

  implicit none
  class( field_laser ), intent(inout) :: this

  integer :: m, i, j, nrp, idproc, noff
  real :: idrh, ir
  real, dimension(:,:), pointer :: ar_re => null(), ar_im => null(), ai_re => null(), ai_im => null()
  character(len=32), save :: sname = 'set_grad_field_laser'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp    = this%cfr_re(0)%get_ndp(1)
  noff   = this%cfr_re(0)%get_noff(1)
  idrh   = 0.5 / this%dr
  idproc = id_proc_loc()

  ! m = 0 mode
  ar_re => this%cfr_re(0)%get_f1()
  ai_re => this%cfi_re(0)%get_f1()
  do i = 2, nrp
    this%ar_grad_re(0)%f1(1,i) = idrh * ( ar_re(1,i+1) - ar_re(1,i-1) )
    this%ai_grad_re(0)%f1(1,i) = idrh * ( ai_re(1,i+1) - ai_re(1,i-1) )
    this%ar_grad_re(0)%f1(2,i) = 0.0
    this%ai_grad_re(0)%f1(2,i) = 0.0
  enddo
  if ( idproc == 0 ) then
    this%ar_grad_re(0)%f1(1,1) = 0.0
    this%ai_grad_re(0)%f1(1,1) = 0.0
    this%ar_grad_re(0)%f1(2,1) = 0.0
    this%ai_grad_re(0)%f1(2,1) = 0.0
  else
    this%ar_grad_re(0)%f1(1,1) = idrh * ( ar_re(1,2) - ar_re(1,0) )
    this%ai_grad_re(0)%f1(1,1) = idrh * ( ai_re(1,2) - ai_re(1,0) )
    this%ar_grad_re(0)%f1(2,1) = 0.0
    this%ai_grad_re(0)%f1(2,1) = 0.0
  endif

  ! DEBUG
  ! print *, "sum(ar_re) = ", sum(ar_re(1,1:nrp))
  ! print *, "sum(ai_re) = ", sum(ai_re(1,1:nrp))
  ! print *, "sum(ar_grad_re) = ", sum(this%ar_grad_re(0)%f1(1,1:nrp))
  ! print *, "sum(ai_grad_re) = ", sum(this%ai_grad_re(0)%f1(1,1:nrp))

  ! m > 0 modes
  do m = 1, this%max_mode

    ar_re => this%cfr_re(m)%get_f1()
    ar_im => this%cfr_im(m)%get_f1()
    ai_re => this%cfi_re(m)%get_f1()
    ai_im => this%cfi_im(m)%get_f1()

    do i = 2, nrp
      ir = 1.0 / ( ( noff + i - 1 ) * this%dr )
      this%ar_grad_re(m)%f1(1,i) = idrh * ( ar_re(1,i+1) - ar_re(1,i-1) )
      this%ar_grad_im(m)%f1(1,i) = idrh * ( ar_im(1,i+1) - ar_im(1,i-1) )
      this%ar_grad_re(m)%f1(2,i) = -ir * m * ar_im(1,i)
      this%ar_grad_im(m)%f1(2,i) =  ir * m * ar_re(1,i)
      this%ai_grad_re(m)%f1(1,i) = idrh * ( ai_re(1,i+1) - ai_re(1,i-1) )
      this%ai_grad_im(m)%f1(1,i) = idrh * ( ai_im(1,i+1) - ai_im(1,i-1) )
      this%ai_grad_re(m)%f1(2,i) = -ir * m * ai_im(1,i)
      this%ai_grad_im(m)%f1(2,i) =  ir * m * ai_re(1,i)
    enddo
    if ( idproc == 0 ) then
      if ( mod(m,2) == 1 ) then
        this%ar_grad_re(m)%f1(1,i) = 2.0 * idrh * ar_re(1,2)
        this%ar_grad_im(m)%f1(1,i) = 2.0 * idrh * ar_im(1,2)
        this%ar_grad_re(m)%f1(2,i) = -m * this%ar_grad_im(m)%f1(1,i)
        this%ar_grad_im(m)%f1(2,i) =  m * this%ar_grad_re(m)%f1(1,i)
        this%ai_grad_re(m)%f1(1,i) = 2.0 * idrh * ai_re(1,2)
        this%ai_grad_im(m)%f1(1,i) = 2.0 * idrh * ai_im(1,2)
        this%ai_grad_re(m)%f1(2,i) = -m * this%ai_grad_im(m)%f1(1,i)
        this%ai_grad_im(m)%f1(2,i) =  m * this%ai_grad_re(m)%f1(1,i)
      else
        this%ar_grad_re(m)%f1(1,i) = 0.0
        this%ar_grad_im(m)%f1(1,i) = 0.0
        this%ar_grad_re(m)%f1(2,i) = 0.0
        this%ar_grad_im(m)%f1(2,i) = 0.0
        this%ai_grad_re(m)%f1(1,i) = 0.0
        this%ai_grad_im(m)%f1(1,i) = 0.0
        this%ai_grad_re(m)%f1(2,i) = 0.0
        this%ai_grad_im(m)%f1(2,i) = 0.0
      endif
    else
      ir = 1.0 / ( noff * this%dr )
      this%ar_grad_re(m)%f1(1,i) = idrh * ( ar_re(1,2) - ar_re(1,0) )
      this%ar_grad_im(m)%f1(1,i) = idrh * ( ar_im(1,2) - ar_im(1,0) )
      this%ar_grad_re(m)%f1(2,i) = -ir * m * ar_im(1,1)
      this%ar_grad_im(m)%f1(2,i) =  ir * m * ar_re(1,1)
      this%ai_grad_re(m)%f1(1,i) = idrh * ( ai_re(1,2) - ai_re(1,0) )
      this%ai_grad_im(m)%f1(1,i) = idrh * ( ai_im(1,2) - ai_im(1,0) )
      this%ai_grad_re(m)%f1(2,i) = -ir * m * ai_im(1,1)
      this%ai_grad_im(m)%f1(2,i) =  ir * m * ai_re(1,1)
    endif

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_grad_field_laser

subroutine solve_field_laser( this, chi, i_slice )

  implicit none
  class( field_laser ), intent(inout) :: this
  class( field ), intent(inout) :: chi
  integer, intent(in) :: i_slice

  integer :: m, k, i, j, nrp, nzp
  real :: dr2_idzh, ds_qtr, rhs
  real, dimension(:,:), pointer :: chi_re => null(), chi_im => null()
  real, dimension(:,:), pointer :: ar_re => null(), ar_im => null()
  real, dimension(:,:), pointer :: ai_re => null(), ai_im => null()
  real, dimension(:,:), pointer :: axir_re => null(), axir_im => null()
  real, dimension(:,:), pointer :: axii_re => null(), axii_im => null()
  real, dimension(:,:,:), pointer :: ar_re_f2 => null(), ar_im_f2 => null()
  real, dimension(:,:,:), pointer :: ai_re_f2 => null(), ai_im_f2 => null()
  real, dimension(:,:,:), pointer :: sr_re_f2 => null(), sr_im_f2 => null()
  real, dimension(:,:,:), pointer :: si_re_f2 => null(), si_im_f2 => null()
  character(len=32), save :: sname = 'solve_field_laser'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp    = this%cfr_re(0)%get_ndp(1)
  nzp    = this%cfr_re(0)%get_ndp(2)
  dr2_idzh = 0.5 * this%dr**2 / this%dz
  ds_qtr = 0.25 * this%ds

  if ( .not. allocated(tmp_r_re) ) then
    allocate( tmp_r_re(nrp), tmp_r_im(nrp), tmp_i_re(nrp), tmp_i_im(nrp) )
  endif
  tmp_r_re = 0.0
  tmp_r_im = 0.0
  tmp_i_re = 0.0
  tmp_i_im = 0.0

  ! handle the first slice, calculate the da/dxi for the first slice
  if ( i_slice == 2 - this%gc_num(1,2) ) then
    do m = 0, this%max_mode
      ar_re_f2 => this%cfr_re(m)%get_f2()
      ai_re_f2 => this%cfi_re(m)%get_f2()
      axir_re => this%axir_re(m)%get_f1()
      axii_re => this%axii_re(m)%get_f1()
      do i = 1, nrp
        axir_re(1,i) = dr2_idzh * ( ar_re_f2(1,i,i_slice+1) - ar_re_f2(1,i,i_slice-1) )
        axii_re(1,i) = dr2_idzh * ( ai_re_f2(1,i,i_slice+1) - ai_re_f2(1,i,i_slice-1) )
      enddo
      if ( m == 0 ) cycle
      ar_im_f2 => this%cfr_im(m)%get_f2()
      ai_im_f2 => this%cfi_im(m)%get_f2()
      axir_im => this%axir_im(m)%get_f1()
      axii_im => this%axii_im(m)%get_f1()
      do i = 1, nrp
        axir_im(1,i) = dr2_idzh * ( ar_im_f2(1,i,i_slice+1) - ar_im_f2(1,i,i_slice-1) )
        axii_im(1,i) = dr2_idzh * ( ai_im_f2(1,i,i_slice+1) - ai_im_f2(1,i,i_slice-1) )
      enddo
    enddo
  endif

  ! calculate the contribution from the plasma susceptibility
  do m = 0, this%max_mode
    do k = 0, m
      chi_re => chi%rf_re(k)%get_f1()
      ar_re  => this%cfr_re(m-k)%get_f1()
      ai_re  => this%cfi_re(m-k)%get_f1()

      do i = 1, nrp
        tmp_r_re(i) = ds_qtr * chi_re(1,i) * ar_re(1,i)
        tmp_i_re(i) = ds_qtr * chi_re(1,i) * ai_re(1,i)
      enddo

      if ( k == 0 .or. k == m ) cycle

      chi_im => chi%rf_im(k)%get_f1()
      ar_im  => this%cfr_im(m-k)%get_f1()
      ai_im  => this%cfi_im(m-k)%get_f1()

      do i = 1, nrp
        tmp_r_re(i) = tmp_r_re(i) - ds_qtr * chi_im(1,i) * ar_im(1,i)
        tmp_i_re(i) = tmp_i_re(i) - ds_qtr * chi_im(1,i) * ai_im(1,i)
      enddo
    enddo

    if ( m == 0 ) cycle

    do k = 0, m

      if ( k /= 0 ) then
        chi_im => chi%rf_im(k)%get_f1()
        ar_re  => this%cfr_re(m-k)%get_f1()
        ai_re  => this%cfi_re(m-k)%get_f1()

        do i = 1, nrp
          tmp_r_im(i) = tmp_r_im(i) + ds_qtr * chi_im(1,i) * ar_re(1,i)
          tmp_i_im(i) = tmp_i_im(i) + ds_qtr * chi_im(1,i) * ai_re(1,i)
        enddo
      endif

      if ( k /= m ) then
        chi_re => chi%rf_re(k)%get_f1()
        ar_im  => this%cfr_im(m-k)%get_f1()
        ai_im  => this%cfi_im(m-k)%get_f1()

        do i = 1, nrp
          tmp_r_im(i) = tmp_r_im(i) + ds_qtr * chi_re(1,i) * ar_im(1,i)
          tmp_i_im(i) = tmp_i_im(i) + ds_qtr * chi_re(1,i) * ai_im(1,i)
        enddo
      endif

    enddo
  enddo

  ! set rhs of the PCR solver and solve
  do m = 0, this%max_mode

    sr_re_f2 => this%sr_re(m)%get_f2()
    si_re_f2 => this%si_re(m)%get_f2()
    ar_re_f2 => this%cfr_re(m)%get_f2()
    ai_re_f2 => this%cfi_re(m)%get_f2()
    axir_re  => this%axir_re(m)%get_f1()
    axii_re  => this%axii_re(m)%get_f1()
    ar_re    => this%cfr_re(m)%get_f1()
    ai_re    => this%cfi_re(m)%get_f1()

    ! set rhs of PCR
    do i = 1, nrp
      rhs = sr_re_f2(1,i,i_slice) - axir_re(1,i) + tmp_r_re(i)
      call this%pgc_solver(m)%set_values_rhs( rhs, 2*i-1 )
      rhs = si_re_f2(1,i,i_slice) - axii_re(1,i) + tmp_i_re(i)
      call this%pgc_solver(m)%set_values_rhs( rhs, 2*i )
    enddo

    call this%pgc_solver(m)%solve()

    ! calculate the da/dxi at next slice before a is replaced
    if ( i_slice < nzp + this%gc_num(2,2) - 1 ) then
      do i = 1, nrp
        axir_re(1,i) = dr2_idzh * ( ar_re_f2(1,i,i_slice+2) - ar_re_f2(1,i,i_slice) )
        axii_re(1,i) = dr2_idzh * ( ai_re_f2(1,i,i_slice+2) - ai_re_f2(1,i,i_slice) )
      enddo
    else
      axir_re = 0.0
      axii_re = 0.0
    endif

    ! get solution
    do i = 1, nrp
      call this%pgc_solver(m)%get_values_x( ar_re(1,i), 2*i-1 )
      call this%pgc_solver(m)%get_values_x( ai_re(1,i), 2*i )
    enddo

    if ( m == 0 ) cycle

    sr_im_f2 => this%sr_im(m)%get_f2()
    si_im_f2 => this%si_im(m)%get_f2()
    ar_im_f2 => this%cfr_im(m)%get_f2()
    ai_im_f2 => this%cfi_im(m)%get_f2()
    axir_im  => this%axir_im(m)%get_f1()
    axii_im  => this%axii_im(m)%get_f1()
    ar_im    => this%cfr_im(m)%get_f1()
    ai_im    => this%cfi_im(m)%get_f1()

    ! set rhs of PCR
    do i = 1, nrp
      rhs = sr_im_f2(1,i,i_slice) - axir_im(1,i) + tmp_r_im(i)
      call this%pgc_solver(m)%set_values_rhs( rhs, 2*i-1 )
      rhs = si_im_f2(1,i,i_slice) - axii_im(1,i) + tmp_i_im(i)
      call this%pgc_solver(m)%set_values_rhs( rhs, 2*i )
    enddo

    call this%pgc_solver(m)%solve()

    ! calculate the da/dxi at next slice before a is replaced
    if ( i_slice < nzp + this%gc_num(2,2) - 1 ) then
      do i = 1, nrp
        axir_im(1,i) = dr2_idzh * ( ar_im_f2(1,i,i_slice+2) - ar_im_f2(1,i,i_slice) )
        axii_im(1,i) = dr2_idzh * ( ai_im_f2(1,i,i_slice+2) - ai_im_f2(1,i,i_slice) )
      enddo
    else
      axir_im = 0.0
      axii_im = 0.0
    endif

    ! get solution
    do i = 1, nrp
      call this%pgc_solver(m)%get_values_x( ar_im(1,i), 2*i-1 )
      call this%pgc_solver(m)%get_values_x( ai_im(1,i), 2*i )
    enddo

  enddo
    
  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_laser

subroutine gather_field_laser( this, that )

  implicit none
  class( field_laser ), intent(inout) :: this
  class( field_laser ), intent(in) :: that

  integer :: m

  call add_f1( that, this )
  do m = 0, this%max_mode
    call add_f1( that%ar_grad_re(m), this%ar_grad_re(m) )
    call add_f1( that%ai_grad_re(m), this%ai_grad_re(m) )
    if ( m == 0 ) cycle
    call add_f1( that%ar_grad_im(m), this%ar_grad_im(m) )
    call add_f1( that%ai_grad_im(m), this%ai_grad_im(m) )
  enddo

  ! DEBUG
  ! print *, "sum(cfr_re(0)) = ", sum( this%cfr_re(0)%f1 )
  ! print *, "sum(cfi_re(0)) = ", sum( this%cfi_re(0)%f1 )

end subroutine gather_field_laser

subroutine zero_field_laser( this, only_f1 )

  implicit none
  class( field_laser ), intent(inout) :: this
  logical, intent(in) :: only_f1

  integer :: m

  ! clean up f1
  do m = 0, this%max_mode

    this%cfr_re(m)%f1 = 0.0
    this%cfi_re(m)%f1 = 0.0
    if ( associated( this%sr_re ) ) then
      this%sr_re(m)%f1 = 0.0
      this%si_re(m)%f1 = 0.0
    endif
    if ( associated( this%axir_re ) ) then
      this%axir_re(m)%f1 = 0.0
      this%axii_re(m)%f1 = 0.0
    endif
    if ( associated( this%ar_grad_re ) ) then
      this%ar_grad_re(m)%f1 = 0.0
      this%ai_grad_re(m)%f1 = 0.0
    endif

    if ( m == 0 ) cycle

    this%cfr_im(m)%f1 = 0.0
    this%cfi_im(m)%f1 = 0.0
    if ( associated( this%sr_im ) ) then
      this%sr_im(m)%f1 = 0.0
      this%si_im(m)%f1 = 0.0
    endif
    if ( associated( this%axir_im ) ) then
      this%axir_im(m)%f1 = 0.0
      this%axii_im(m)%f1 = 0.0
    endif
    if ( associated( this%ar_grad_im ) ) then
      this%ar_grad_im(m)%f1 = 0.0
      this%ai_grad_im(m)%f1 = 0.0
    endif

  enddo

  if ( only_f1 ) return

  ! clean up f2
  do m = 0, this%max_mode

    if ( this%cfr_re(m)%has2d() ) then
      this%cfr_re(m)%f2 = 0.0
      this%cfi_re(m)%f2 = 0.0
    endif
    if ( associated( this%sr_re ) .and. this%sr_re(m)%has2d() ) then
      this%sr_re(m)%f2 = 0.0
      this%si_re(m)%f2 = 0.0
    endif
    if ( associated( this%axir_re ) .and. this%axir_re(m)%has2d() ) then
      this%axir_re(m)%f2 = 0.0
      this%axii_re(m)%f2 = 0.0
    endif
    if ( associated( this%ar_grad_re ) .and. this%ar_grad_re(m)%has2d() ) then
      this%ar_grad_re(m)%f2 = 0.0
      this%ai_grad_re(m)%f2 = 0.0
    endif

    if ( m == 0 ) cycle

    if ( this%cfr_im(m)%has2d() ) then
      this%cfr_im(m)%f2 = 0.0
      this%cfi_im(m)%f2 = 0.0
    endif
    if ( associated( this%sr_im ) ) then
      this%sr_im(m)%f2 = 0.0
      this%si_im(m)%f2 = 0.0
    endif
    if ( associated( this%axir_im ) .and. this%axir_im(m)%has2d() ) then
      this%axir_im(m)%f2 = 0.0
      this%axii_im(m)%f2 = 0.0
    endif
    if ( associated( this%ar_grad_im ) .and. this%ar_grad_im(m)%has2d() ) then
      this%ar_grad_im(m)%f2 = 0.0
      this%ai_grad_im(m)%f2 = 0.0
    endif

  enddo

end subroutine zero_field_laser

end module field_laser_class