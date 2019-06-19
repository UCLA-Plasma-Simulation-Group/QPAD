module field_psi_class

use parallel_pipe_class
use grid_class
use field_class
use field_solver_class
use field_src_class
use ufield_class
use param
use sys
use mpi

implicit none

private

character(len=20), parameter :: cls_name = "field_psi"
integer, parameter :: cls_level = 3

public :: field_psi

type, extends( field ) :: field_psi

  ! private

  class( field_solver ), dimension(:), pointer :: solver => null()
  real, dimension(:), pointer :: buf_re => null(), buf_im => null() ! buffer for source term

  contains

  generic :: new => init_field_psi
  procedure :: del => end_field_psi
  ! generic :: read_input => read_input_field_psi
  generic :: solve => solve_field_psi

  procedure, private :: init_field_psi
  procedure, private :: end_field_psi
  procedure, private :: set_source
  procedure, private :: get_solution
  procedure, private :: solve_field_psi

end type field_psi

contains

subroutine init_field_psi( this, pp, gp, num_modes, part_shape, boundary, iter_tol )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape, boundary
  real, intent(in) :: iter_tol

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  real :: dr
  character(len=20), save :: sname = "init_field_psi"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = gp%get_ndp(1)
  dr = gp%get_dr()

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/1, 1/)
    gc_num(:,2) = (/2, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 1
  ! call initialization routine of the parent class
  call this%field%new( pp, gp, dim, num_modes, gc_num, entity=p_entity_plasma )

  ! initialize solver
  allocate( this%solver( 0:num_modes ) )
  do i = 0, num_modes
    call this%solver(i)%new( pp, gp, i, dr, &
      kind=p_fk_psi, bnd=boundary, stype=p_hypre_cycred, tol=iter_tol )
  enddo

  allocate( this%buf_re(nrp), this%buf_im(nrp) )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_psi

subroutine end_field_psi( this )

  implicit none

  class( field_psi ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = 'end_field_psi'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%solver(i)%del()
  enddo
  deallocate( this%solver )

  if ( associated( this%buf_re ) ) deallocate( this%buf_re )
  if ( associated( this%buf_im ) ) deallocate( this%buf_im )

  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_psi

subroutine set_source( this, mode, q_re, q_im )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( ufield ), intent(in) :: q_re
  class( ufield ), intent(in), optional :: q_im
  integer, intent(in) :: mode

  integer :: i, nrp, noff, dtype, ierr, comm, idproc
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, save :: local_sum, global_sum
  real :: dr, dr2, rmax
  character(len=20), save :: sname = 'set_source'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve psi' )

  nrp   = q_re%get_ndp(1)
  noff  = q_re%get_noff(1)
  dr    = this%dr
  dr2   = dr**2
  dtype = q_re%pp%getmreal()
  comm  = q_re%pp%getlgrp()
  idproc= q_re%pp%getlidproc()
  rmax  = q_re%get_nd(1) * dr

  f1_re => q_re%get_f1()
  if ( present(q_im) ) then
    f1_im => q_im%get_f1()
  endif

  this%buf_re = 0.0
  if ( present(q_im) ) this%buf_im = 0.0
  if ( mode == 0 ) then

    select case ( this%solver(0)%bnd )

    case ( p_bnd_zero, p_bnd_open )
      do i = 1, nrp
        this%buf_re(i) = -1.0 * f1_re(1,i)
      enddo

      ! if ( idproc == 0 ) this%q_ax = -1.0*dr2*f1_re(1,0)
      ! call MPI_BCAST( this%q_ax, 1, dtype, 0, comm, ierr )

    end select

  elseif ( mode > 0 .and. present(q_im) ) then

    do i = 1, nrp
      this%buf_re(i) = -1.0 * f1_re(1,i)
      this%buf_im(i) = -1.0 * f1_im(1,i)
    enddo

  else
    call write_err( 'Invalid input arguments!' )
  endif

  call stop_tprof( 'solve psi' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source

subroutine get_solution( this, mode, q_ax )

  implicit none

  class( field_psi ), intent(inout) :: this
  integer, intent(in) :: mode
  real, intent(in), optional :: q_ax

  integer :: i, nrp, noff
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: r, r2, rmax, rmax2
  character(len=20), save :: sname = 'get_solution'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve psi' )

  nrp   = this%rf_re(mode)%get_ndp(1)
  noff  = this%rf_re(mode)%get_noff(1)
  rmax  = this%rf_re(mode)%get_nd(1) * this%dr
  rmax2 = rmax**2

  if ( mode == 0 .and. present(q_ax) ) then
    select case ( this%solver(0)%bnd )

    case ( p_bnd_zero, p_bnd_open )

      ! add contribution of the source terms on axis
      do i = 1, nrp
        r = (i+noff-0.5) * this%dr
        this%buf_re(i) = this%buf_re(i) + q_ax * log(r/rmax)
      enddo

    end select
  endif

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nrp
    f1_re(1,i) = this%buf_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nrp
      f1_im(1,i) = this%buf_im(i)
    enddo
  endif

  call stop_tprof( 'solve psi' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution

subroutine solve_field_psi( this, q )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( field_rho ), intent(in) :: q

  type( ufield ), dimension(:), pointer :: q_re => null(), q_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_psi'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  q_re => q%get_rf_re()
  q_im => q%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source( i, q_re(i) )
      call this%solver(i)%solve( this%buf_re )
      call this%get_solution( i, q%q_ax )
      cycle
    endif

    call this%set_source( i, q_re(i), q_im(i) )
    call this%solver(i)%solve( this%buf_re )
    call this%solver(i)%solve( this%buf_im )
    call this%get_solution(i)

  enddo

  call this%copy_gc_f1( bnd_ax = .true. )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_psi

end module field_psi_class