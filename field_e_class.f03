module field_e_class

use parallel_pipe_class
use grid_class
use field_class
use field_b_class
use field_psi_class
use field_src_class
use field_solver_class
use ufield_class
use param
use system

implicit none

private

character(len=20), parameter :: cls_name = "field_e"
integer, parameter :: cls_level = 3

public :: field_e

type, extends( field ) :: field_e

  ! private

  class( field_solver ), dimension(:), pointer :: solver_ez => null()
  real, dimension(:), pointer :: buf_re => null(), buf_im => null()

  contains

  generic :: new => init_field_e
  procedure :: del => end_field_e
  ! generic :: read_input => read_input_field_e
  generic :: solve => solve_field_ez, solve_field_eperp, solve_field_eperp_beam

  procedure, private :: init_field_e
  procedure, private :: end_field_e
  procedure, private :: set_source_ez
  procedure, private :: get_solution_ez
  procedure, private :: solve_field_ez
  procedure, private :: solve_field_eperp
  procedure, private :: solve_field_eperp_beam

end type field_e

contains

subroutine init_field_e( this, pp, gp, dr, dxi, num_modes, part_shape, entity )

  implicit none

  class( field_e ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape, entity
  real, intent(in) :: dr, dxi

  integer, dimension(2,2) :: gc_num
  integer :: dim, i
  integer, dimension(2) :: ndp, noff, nd
  real :: tol
  character(len=20), save :: sname = "init_field_e"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd = gp%get_nd()
  ndp = gp%get_ndp()
  noff = gp%get_noff()
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
  call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num, entity )

  ! initialize solver
  select case ( entity )
  case ( p_entity_plasma )
    allocate( this%solver_ez( 0:num_modes ) )
    do i = 0, num_modes
      call this%solver_ez(i)%new( nd, ndp, noff, p_fk_ez, i, dr, p_hypre_cycred, tol )
    enddo 
  case ( p_entity_beam )
    ! do nothing
  case default
    call write_err( 'Invalid field entity type.' )
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_e

subroutine end_field_e( this )

  implicit none

  class( field_e ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = 'end_field_e'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( this%entity == p_entity_plasma ) then
    do i = 0, this%num_modes
      call this%solver_ez(i)%del()
    enddo
    deallocate( this%solver_ez )
  endif

  if ( associated( this%buf_re ) ) deallocate( this%buf_re )
  if ( associated( this%buf_im ) ) deallocate( this%buf_im )

  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_e

subroutine set_source_ez( this, mode, jay_re, jay_im )

  implicit none

  class( field_e ), intent(inout) :: this
  integer, intent(in) :: mode
  class( ufield ), intent(in) :: jay_re
  class( ufield ), intent(in), optional :: jay_im

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: idrh, idr, a1, a2, a3, b
  character(len=20), save :: sname = 'set_source_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = jay_re%get_ndp(1)
  idr = 1.0 / this%dr
  idrh = 0.5 * idr

  f1_re => jay_re%get_f1()
  if ( .not. associated( this%buf_re ) ) then
    allocate( this%buf_re( nd1p ) )
  elseif ( size(this%buf_re) < nd1p ) then
    deallocate( this%buf_re )
    allocate( this%buf_re( nd1p ) )
  endif

  if ( present(jay_im) ) then
    f1_im => jay_im%get_f1()
    if ( .not. associated( this%buf_im ) ) then
      allocate( this%buf_im( nd1p ) )
    elseif ( size(this%buf_im) < nd1p ) then
      deallocate( this%buf_im )
      allocate( this%buf_im( nd1p ) )
    endif
  endif

  if ( mode == 0 ) then
    
    do i = 2, nd1p-1

      a1 = -idrh * (i-1.0) / (i-0.5)
      a2 =  idrh / (i-0.5)
      a3 =  idrh * i / (i-0.5)

      this%buf_re(i) = a1 * f1_re(1,i-1) + a2 * f1_re(1,i) + a3 * f1_re(1,i+1)

    enddo

    ! calculate the derivatives at the boundary and axis
    this%buf_re(1) = idr * ( f1_re(1,1) + f1_re(1,2) )
    a2 = idr * (nd1p+0.5) / (nd1p-0.5)
    this%buf_re(nd1p) = -idr * f1_re(1,nd1p-1) + a2 * f1_re(1,nd1p)

  elseif ( mode > 0 .and. present( jay_im ) ) then
    
    do i = 2, nd1p-1

      a1 = -idrh * (i-1.0) / (i-0.5)
      a2 = idrh / (i-0.5)
      a3 = idrh * i / (i-0.5)
      b  = idr * real(mode) / (i-0.5)

      this%buf_re(i) = a1 * f1_re(1,i-1) + a2 * f1_re(1,i) + a3 * f1_re(1,i+1) - &
                        b * f1_im(2,i)
      this%buf_im(i) = a1 * f1_im(1,i-1) + a2 * f1_im(1,i) + a3 * f1_im(1,i+1) + &
                        b * f1_re(2,i)

    enddo

    ! calculate the derivatives at the boundary and axis????????????????????????????????????
    this%buf_re(1) = idr * ( f1_re(2,1) + f1_re(2,2) - 2.0 * real(mode) * f1_im(1,1) )
    this%buf_im(1) = idr * ( f1_im(2,1) + f1_im(2,2) + 2.0 * real(mode) * f1_re(1,1) )
    a2 = idr * (nd1p+0.5) / (nd1p-0.5)
    b  = idr * real(mode) / (nd1p-0.5)
    this%buf_re(nd1p) = -idr * f1_re(2,nd1p-1) + a2 * f1_re(2,nd1p) - b * f1_im(1,nd1p)
    this%buf_im(nd1p) = -idr * f1_im(2,nd1p-1) + a2 * f1_im(2,nd1p) + b * f1_re(1,nd1p)

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_ez

subroutine get_solution_ez( this, mode )

  implicit none

  class( field_e ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = this%rf_re(mode)%get_ndp(1)
  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nd1p
    f1_re(3,i) = this%buf_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1_im(3,i) = this%buf_im(i)
    enddo
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_ez

subroutine solve_field_ez( this, jay )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_jay ), intent(in) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_ez( i, jay_re(i) )
      call this%solver_ez(i)%solve( this%buf_re )
      call this%get_solution_ez(i)
      cycle
    endif

    call this%set_source_ez( i, jay_re(i), jay_im(i) )
    call this%solver_ez(i)%solve( this%buf_re )
    call this%solver_ez(i)%solve( this%buf_im )
    call this%get_solution_ez(i)

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_ez

subroutine solve_field_eperp( this, b, psi )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_b ), intent(in) :: b
  class( field_psi ), intent(in) :: psi

  type( ufield ), dimension(:), pointer :: b_re => null(), b_im => null()
  type( ufield ), dimension(:), pointer :: psi_re => null(), psi_im => null()
  real, dimension(:,:), pointer :: ub_re => null(), ub_im => null()
  real, dimension(:,:), pointer :: upsi_re => null(), upsi_im => null()
  real, dimension(:,:), pointer :: ue_re => null(), ue_im => null()
  integer :: mode, i, nd1p
  real :: idr, idrh, ir
  character(len=20), save :: sname = 'solve_field_eperp'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  idr = 1.0 / this%dr
  idrh = idr * 0.5
  nd1p = this%rf_re(0)%get_ndp(1)

  b_re => b%get_rf_re()
  b_im => b%get_rf_im()
  psi_re => psi%get_rf_re()
  psi_im => psi%get_rf_im()

  do mode = 0, this%num_modes

    ub_re => b_re(mode)%get_f1()
    upsi_re => psi_re(mode)%get_f1()
    ue_re => this%rf_re(mode)%get_f1()
    if ( mode == 0 ) then
      do i = 2, nd1p-1
        ue_re(1,i) = ub_re(2,i) - idrh * ( upsi_re(1,i+1) - upsi_re(1,i-1) )
        ue_re(2,i) = -ub_re(1,i)
      enddo
      ue_re(1,1) = ub_re(2,1) - idr * ( upsi_re(1,2) - upsi_re(1,1) )
      ue_re(2,1) = -ub_re(1,1)
      ue_re(1,nd1p) = ub_re(2,nd1p) - idr * ( upsi_re(1,nd1p) - upsi_re(1,nd1p-1) )
      ue_re(2,nd1p) = -ub_re(1,nd1p)
      cycle
    endif

    ub_im => b_im(mode)%get_f1()
    upsi_im => psi_im(mode)%get_f1()
    ue_im => this%rf_im(mode)%get_f1()

    do i = 2, nd1p-1
      ir = idr / (i-0.5)
      ue_re(1,i) = ub_re(2,i) - idrh * ( upsi_re(1,i+1) - upsi_re(1,i-1) )
      ue_re(2,i) = -ub_re(1,i) + ir * mode * upsi_im(1,i)

      ue_im(1,i) = ub_im(2,i) - idrh * ( upsi_im(1,i+1) - upsi_im(1,i-1) )
      ue_im(2,i) = -ub_im(1,i) - ir * mode * upsi_re(1,i)
    enddo
    ir = 2.0 * idr
    ue_re(1,1) = ub_re(2,1) - idr * ( upsi_re(1,2) - upsi_re(1,1) )
    ue_re(2,1) = -ub_re(1,1) + ir * mode * upsi_im(1,1)
    ue_im(1,1) = ub_im(2,1) - idr * ( upsi_im(1,2) - upsi_im(1,1) )
    ue_im(2,1) = -ub_im(1,1) - ir * mode * upsi_re(1,1)
    ir = idr / (nd1p-0.5)
    ue_re(1,nd1p) = ub_re(2,nd1p) - idr * ( upsi_re(1,nd1p) - upsi_re(1,nd1p-1) )
    ue_re(2,nd1p) = -ub_re(1,nd1p) + ir * mode * upsi_im(1,nd1p)
    ue_im(1,nd1p) = ub_im(2,nd1p) - idr * ( upsi_im(1,nd1p) - upsi_im(1,nd1p-1) )
    ue_im(2,nd1p) = -ub_im(1,nd1p) - ir * mode * upsi_re(1,nd1p)

  enddo

end subroutine solve_field_eperp

subroutine solve_field_eperp_beam( this, b )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_b ), intent(in) :: b

  type( ufield ), dimension(:), pointer :: b_re => null(), b_im => null()
  real, dimension(:,:), pointer :: ub_re => null(), ub_im => null()
  real, dimension(:,:), pointer :: ue_re => null(), ue_im => null()
  integer :: mode, i, nd1p
  character(len=20), save :: sname = 'solve_field_eperp_beam'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  b_re => b%get_rf_re()
  b_im => b%get_rf_im()
  nd1p = this%rf_re(0)%get_ndp(1)

  do mode = 0, this%num_modes

    ub_re => b_re(mode)%get_f1()
    ue_re => this%rf_re(mode)%get_f1()

    do i = 1, nd1p
      ue_re(1,i) = ub_re(2,i)
      ue_re(2,i) = -ub_re(1,i)
    enddo

    if ( mode == 0 ) cycle

    ub_im => b_im(mode)%get_f1()
    ue_im => this%rf_im(mode)%get_f1()

    do i = 1, nd1p
      ue_im(1,i) = ub_im(2,i)
      ue_im(2,i) = -ub_im(1,i)
    enddo

  enddo

end subroutine solve_field_eperp_beam

end module field_e_class