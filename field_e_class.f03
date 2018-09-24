module field_e_class

use field_class
use field_solver_class
use ufield_class
use param
use system

implicit none

private

character(len=20), parameter :: cls_name = "field_e"
integer, parameter :: cls_level = 1

public :: field_e

type, extends( field ) :: field_e

  ! private

  class( field_solver ), dimension(:), pointer :: solver => null()
  real, dimension(:), pointer :: buf => null() ! buffer for source term

  contains

  generic :: new => init_field_e
  procedure :: del => end_field_e
  ! generic :: read_input => read_input_field_e
  generic :: solve => solve_field_ez, solve_field_eperp

  procedure, private :: init_field_e
  procedure, private :: end_field_e
  procedure, private :: sort_src
  procedure, private :: sort_sol
  procedure, private :: solve_field_ez
  procedure, private :: solve_field_eperp

end type field_e

contains

subroutine init_field_e( this, num_modes, dr, dxi, nd, nvp, order, part_shape )

  implicit none

  class( field_e ), intent(inout) :: this
  integer, intent(in) :: num_modes, order, part_shape
  real, intent(in) :: dr, dxi
  integer, intent(in), dimension(2) :: nd, nvp

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, solver_type, kind
  integer, dimension(2) :: ndp, noff
  real :: tol
  character(len=20), save :: sname = "init_field_e"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ndp = nd / nvp
  noff = (/0,0/)
  solver_type = p_hypre_cycred
  kind = p_fk_e
  tol = 1.0d-6

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/0, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    print *, "Quadratic particle shape not implemented."
    stop

  case default

    print *, "Invalid particle shape."
    stop

  end select

  dim = 3
  ! call initialization routine of the parent class
  call this%field%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )

  ! initialize solver
  allocate( this%solver( 0:num_modes ) )
  do i = 0, num_modes
    call this%solver(i)%new( ndp, noff, order, kind, i, dr, solver_type, tol )
  enddo 

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_e

subroutine end_field_e( this )

  implicit none

  class( field_e ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = 'end_field_e'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%solver(i)%del()
  enddo
  deallocate( this%solver )

  if ( associated( this%buf ) ) deallocate( this%buf )

  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_e

subroutine sort_src( this, jay )

  implicit none

  class( field_e ), intent(inout) :: this
  class( ufield ), intent(in) :: jay

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1 => null()
  character(len=20), save :: sname = 'sort_src'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = jay%get_ndp(1)
  if ( .not. associated( this%buf ) ) then
    allocate( this%buf( nd1p ) )
  endif

  f1 => q%get_f1()
  do i = 1, nd1p
    this%buf(i) = f1(3,i)
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine sort_src

subroutine sort_sol( this, mode, part )

  implicit none

  class( field_e ), intent(inout) :: this
  integer, intent(in) :: mode, part

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1 => null()
  character(len=20), save :: sname = 'sort_sol'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = this%rf_re(mode)%get_ndp(1)
  if ( part == p_real ) then

    f1 => this%rf_re(mode)%get_f1()
    do i = 1, nd1p
      f1(3,i) = this%buf(i)
    enddo

  else

    f1 => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1(3,i) = this%buf(i)
    enddo

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine sort_sol

subroutine solve_field_ez( this, jay )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field ), intent(in) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_e'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%num_modes

    call this%sort_src( jay_re(i) )
    call this%solver(i)%solve( this%buf )
    call this%sort_sol( i, part=p_real )

    if ( i == 0 ) cycle

    call this%sort_src( jay_im(i) )
    call this%solver(i)%solve( this%buf )
    call this%sort_sol( i, part=p_imag )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_ez

subroutine solve_field_eperp( this, psi, b )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_psi ), intent(in) :: psi
  class( field_b ), intent(in) :: b



end subroutine solve_field_eperp

end module field_e_class