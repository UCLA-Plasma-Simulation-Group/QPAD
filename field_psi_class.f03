module field_psi_class

use field_class
use field_solver_class
use ufield_class
use param

implicit none

private

type, extends( field ) :: field_psi

  private

  class( field_solver ), dimension(:), pointer :: solver => null()

  contains

  generic :: new => init_field_psi
  procedure :: del => end_field_psi
  ! generic :: read_input => read_input_field_psi
  generic :: solve => solve_field_psi

  procedure, private :: init_field_psi
  procedure, private :: end_field_psi
  procedure, private :: sort_src
  procedure, private :: sort_sol
  procedure, private :: solve_field_psi

end type field_psi

contains

subroutine init_field_psi( this, num_modes, dr, dxi, nd, nvp, order, part_shape )

  implicit none

  class( field_psi ), intent(inout) :: this
  integer, intent(in) :: num_modes, order, part_shape
  real, intent(in) :: dr, dxi
  integer, intent(in), dimension(2) :: nd, nvp

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, solver_type, precond_type
  integer, dimension(2) :: ndp, noff
  real :: tol

  ndp = nd / nvp
  noff = (/0,0/)
  solver_type = p_hypre_pcg
  precond_type = p_hypre_smg
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

  dim = 1
  ! call initialization routine of the parent class
  call this%field%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )

  ! initialize solver
  allocate( this%solver( 0:num_modes ) )
  do i = 0, num_modes
    call this%solver(i)%new( ndp, noff, order, i, dr, solver_type, &
      precond_type, tol )
  enddo 

end subroutine init_field_psi

subroutine end_field_psi( this )

  implicit none

  class( field_psi ), intent(inout) :: this

  integer :: i

  do i = 0, this%num_modes
    call this%solver(i)%del()
  enddo
  deallocate( this%solver )

  call this%field%del()

end subroutine end_field_psi

subroutine sort_src( this, q )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( ufield ), intent(in) :: q

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1 => null()

  if ( .not. associated( this%buf ) ) then
    nd1p = q%get_ndp(1)
    allocate( this%buf( nd1p ) )
  endif

  f1 => q%get_f1()
  do i = 1, nd1p
    this%buf(i) = f1(1,i)
  enddo

end subroutine sort_src

subroutine sort_sol( this, mode, part )

  implicit none

  class( field_psi ), intent(inout) :: this
  integer, intent(in) :: mode, part

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1 => null()

  nd1p = this%rf_re(mode)%get_ndp(1)
  if ( part == p_real ) then

    f1 => this%rf_re(mode)%get_f1()
    do i = 1, nd1p
      f1(1,i) = this%buf(i)
    enddo

  else

    f1 => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1(1,i) = this%buf(i)
    enddo

  endif

end subroutine sort_sol

subroutine solve_field_psi( this, q )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( field ), intent(in) :: q

  class( ufield ), dimension(:), pointer :: q_re => null(), q_im => null()
  integer :: i

  q_re => q%get_rf_re()
  q_im => q%get_rf_im()

  do i = 0, this%num_modes

    call this%sort_src( q_re(i) )
    call this%solver(i)%solve( this%buf )
    call this%sort_sol( i, part=p_real )

    if ( i == 0 ) cycle

    call this%sort_src( q_im(i) )
    call this%solver(i)%solve( this%buf )
    call this%sort_sol( i, part=p_imag )

  enddo

end subroutine solve_field_psi

end module field_psi_class