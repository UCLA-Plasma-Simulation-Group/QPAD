module field_psi_class

use parallel_module
use options_class
use field_class
use field_solver_class
use field_src_class
use ufield_class
use param
use sysutil_module

implicit none

private

character(len=32), parameter :: cls_name = "field_psi"
integer, parameter :: cls_level = 3

public :: field_psi

type, extends( field ) :: field_psi

  ! private

  class( field_solver ), dimension(:), pointer :: solver => null()
  real, dimension(:), pointer :: buf_re => null(), buf_im => null() ! buffer for source term

  contains

  generic :: new => init_field_psi
  procedure :: init_field_psi
  procedure :: del => end_field_psi
  procedure :: solve => solve_field_psi
  procedure, private :: set_source
  procedure, private :: get_solution

end type field_psi

contains

subroutine init_field_psi( this, opts, num_modes, part_shape, boundary )

  implicit none

  class( field_psi ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: num_modes, part_shape, boundary

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  real :: dr
  character(len=32), save :: sname = "init_field_psi"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = opts%get_ndp(1)
  dr  = opts%get_dr()

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
  call this%field%new( opts, dim, num_modes, gc_num, entity=p_entity_plasma )

  ! initialize solver
  allocate( this%solver( 0:num_modes ) )
  do i = 0, num_modes
    call this%solver(i)%new( opts, i, dr, &
      kind=p_fk_psi, bnd=boundary, stype=p_hypre_cycred )
  enddo

  allocate( this%buf_re(nrp), this%buf_im(nrp) )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_psi

subroutine end_field_psi( this )

  implicit none

  class( field_psi ), intent(inout) :: this

  integer :: i
  character(len=32), save :: sname = 'end_field_psi'

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

  integer :: i, nrp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=32), save :: sname = 'set_source'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve psi' )

  nrp = q_re%get_ndp(1)

  f1_re => q_re%get_f1()
  this%buf_re = 0.0
  if ( present(q_im) ) then
    f1_im => q_im%get_f1()
    this%buf_im = 0.0
  endif

  if ( mode == 0 ) then

    do i = 1, nrp
      this%buf_re(i) = -1.0 * f1_re(1,i)
    enddo

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

subroutine get_solution( this, mode )

  implicit none

  class( field_psi ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp, idproc
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=32), save :: sname = 'get_solution'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve psi' )

  nrp    = this%rf_re(mode)%get_ndp(1)
  idproc = id_proc_loc()

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nrp
    f1_re(1,i) = this%buf_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nrp
      f1_im(1,i) = this%buf_im(i)
    enddo

    ! psi(m>0) vanishes on axis
    if ( idproc == 0 ) then
      f1_re(1,1) = 0.0
      f1_im(1,1) = 0.0
    endif
  endif

  call stop_tprof( 'solve psi' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution

subroutine solve_field_psi( this, q )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( field_rho ), intent(inout) :: q

  type( ufield ), dimension(:), pointer :: q_re => null(), q_im => null()
  integer :: mode
  character(len=32), save :: sname = 'solve_field_psi'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  q_re => q%get_rf_re()
  q_im => q%get_rf_im()

  do mode = 0, this%num_modes

    if ( mode == 0 ) then
      call this%set_source( mode, q_re(mode) )
      call this%solver(mode)%solve( this%buf_re )
      call this%get_solution(mode)
      cycle
    endif

    call this%set_source( mode, q_re(mode), q_im(mode) )
    call this%solver(mode)%solve( this%buf_re )
    call this%solver(mode)%solve( this%buf_im )
    call this%get_solution(mode)

  enddo

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_psi

end module field_psi_class