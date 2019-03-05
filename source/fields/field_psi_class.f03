module field_psi_class

use parallel_pipe_class
use grid_class
use field_class
use field_solver_class
use field_src_class
use ufield_class
use param
use system

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

subroutine init_field_psi( this, pp, gp, dr, dxi, num_modes, part_shape )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape
  real, intent(in) :: dr, dxi

  integer, dimension(2,2) :: gc_num
  integer :: dim, i
  integer, dimension(2) :: ndp, noff, nd
  character(len=20), save :: sname = "init_field_psi"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd = gp%get_nd()
  ndp = gp%get_ndp()
  noff = gp%get_noff()

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/1, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 1
  ! call initialization routine of the parent class
  call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num, entity=p_entity_plasma )

  ! initialize solver
  allocate( this%solver( 0:num_modes ) )
  do i = 0, num_modes
    call this%solver(i)%new( pp, gp, i, dr, &
      kind=p_fk_psi, stype=p_hypre_cycred, tol=1.0d-6 )
  enddo

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

  integer :: i, nrp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'set_source'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = q_re%get_ndp(1)

  f1_re => q_re%get_f1()
  if ( .not. associated( this%buf_re ) ) then
    allocate( this%buf_re( nrp ) )
  elseif ( size(this%buf_re) < nrp ) then
    deallocate( this%buf_re )
    allocate( this%buf_re( nrp ) )
  endif

  if ( present(q_im) ) then
    f1_im => q_im%get_f1()
    if ( .not. associated( this%buf_im ) ) then
      allocate( this%buf_im( nrp ) )
    elseif ( size(this%buf_im) < nrp ) then
      deallocate( this%buf_im )
      allocate( this%buf_im( nrp ) )
    endif
  endif

  this%buf_re = 0.0
  if ( present(q_im) ) this%buf_im = 0.0
  if ( mode == 0 ) then
    do i = 1, nrp
      this%buf_re(i) = f1_re(1,i)
    enddo
  elseif ( mode > 0 .and. present(q_im) ) then
    do i = 1, nrp
      this%buf_re(i) = f1_re(1,i)
      this%buf_im(i) = f1_im(1,i)
    enddo
  else
    call write_err( 'Invalid input arguments!' )
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source

subroutine get_solution( this, mode )

  implicit none

  class( field_psi ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = this%rf_re(mode)%get_ndp(1)


  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nd1p
    f1_re(1,i) = this%buf_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1_im(1,i) = this%buf_im(i)
    enddo
  endif

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
      call this%get_solution(i)
      cycle
    endif

    call this%set_source( i, q_re(i), q_im(i) )
    call this%solver(i)%solve( this%buf_re )
    call this%solver(i)%solve( this%buf_im )
    call this%get_solution(i)

  enddo

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_psi

end module field_psi_class