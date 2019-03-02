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
  character(len=20), save :: sname = "init_field_e"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd = gp%get_nd()
  ndp = gp%get_ndp()
  noff = gp%get_noff()

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/0, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 3
  ! call initialization routine of the parent class
  call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num, entity )

  ! initialize solver
  select case ( entity )
  case ( p_entity_plasma )
    allocate( this%solver_ez( 0:num_modes ) )
    do i = 0, num_modes
      call this%solver_ez(i)%new( pp, gp, i, dr, kind=p_fk_ez, &
        stype=p_hypre_cycred, tol=1.0d-6 )
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

  integer :: i, nrp, noff, idproc, nvp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: idrh, idr, k0, a1, a2, a3, b
  character(len=20), save :: sname = 'set_source_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = jay_re%get_ndp(1)
  idr = 1.0 / this%dr
  idrh = 0.5 * idr
  noff = jay_re%get_noff(1)
  nvp = jay_re%pp%getlnvp()
  idproc = jay_re%pp%getlidproc()

  f1_re => jay_re%get_f1()
  if ( .not. associated( this%buf_re ) ) then
    allocate( this%buf_re( nrp ) )
  elseif ( size(this%buf_re) < nrp ) then
    deallocate( this%buf_re )
    allocate( this%buf_re( nrp ) )
  endif

  if ( present(jay_im) ) then
    f1_im => jay_im%get_f1()
    if ( .not. associated( this%buf_im ) ) then
      allocate( this%buf_im( nrp ) )
    elseif ( size(this%buf_im) < nrp ) then
      deallocate( this%buf_im )
      allocate( this%buf_im( nrp ) )
    endif
  endif

  this%buf_re = 0.0
  if ( present(jay_im) ) this%buf_im = 0.0
  if ( mode == 0 ) then
    
    do i = 1, nrp

      k0 = real(i+noff) - 0.5
      a1 = -idrh * (k0-0.5) / k0
      a2 =  idrh / k0
      a3 =  idrh * (k0+0.5) / k0

      this%buf_re(i) = a1 * f1_re(1,i-1) + a2 * f1_re(1,i) + a3 * f1_re(1,i+1)

    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf_re(1) = idr * ( f1_re(1,1) + f1_re(1,2) )
    endif
    if ( idproc == nvp-1 ) then
      a2 = idr * (nrp+noff+0.5) / (nrp+noff-0.5)
      this%buf_re(nrp) = -idr * f1_re(1,nrp-1) + a2 * f1_re(1,nrp)
    endif

  elseif ( mode > 0 .and. present( jay_im ) ) then
    
    do i = 1, nrp

      k0 = real(i+noff) - 0.5
      a1 = -idrh * (k0-0.5) / k0
      a2 = idrh / k0
      a3 = idrh * (k0+0.5) / k0
      b  = idr * real(mode) / k0

      this%buf_re(i) = a1 * f1_re(1,i-1) + a2 * f1_re(1,i) + a3 * f1_re(1,i+1) - &
                        b * f1_im(2,i)
      this%buf_im(i) = a1 * f1_im(1,i-1) + a2 * f1_im(1,i) + a3 * f1_im(1,i+1) + &
                        b * f1_re(2,i)

    enddo

    ! calculate the derivatives at the boundary and axis????????????????????????????????????
    if ( idproc == 0 ) then
      this%buf_re(1) = idr * ( f1_re(2,1) + f1_re(2,2) - 2.0 * real(mode) * f1_im(1,1) )
      this%buf_im(1) = idr * ( f1_im(2,1) + f1_im(2,2) + 2.0 * real(mode) * f1_re(1,1) )
    endif
    if ( idproc == nvp-1 ) then
      k0 = real(nrp+noff) - 0.5
      a2 = idr * (k0+1.0) / k0
      b  = idr * real(mode) / k0
      this%buf_re(nrp) = -idr * f1_re(2,nrp-1) + a2 * f1_re(2,nrp) - b * f1_im(1,nrp)
      this%buf_im(nrp) = -idr * f1_im(2,nrp-1) + a2 * f1_im(2,nrp) + b * f1_re(1,nrp)
    endif

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
  class( field_jay ), intent(inout) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call jay%copy_gc_f1()

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
  class( field_psi ), intent(inout) :: psi

  type( ufield ), dimension(:), pointer :: b_re => null(), b_im => null()
  type( ufield ), dimension(:), pointer :: psi_re => null(), psi_im => null()
  real, dimension(:,:), pointer :: ub_re => null(), ub_im => null()
  real, dimension(:,:), pointer :: upsi_re => null(), upsi_im => null()
  real, dimension(:,:), pointer :: ue_re => null(), ue_im => null()
  integer :: mode, i, nrp, noff, idproc, nvp
  real :: idr, idrh, ir, k0
  character(len=20), save :: sname = 'solve_field_eperp'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call psi%copy_gc_f1()

  idr = 1.0 / this%dr
  idrh = idr * 0.5
  nrp = this%rf_re(0)%get_ndp(1)

  b_re => b%get_rf_re()
  b_im => b%get_rf_im()
  psi_re => psi%get_rf_re()
  psi_im => psi%get_rf_im()

  noff = this%rf_re(0)%get_noff(1)
  nvp = this%rf_re(0)%pp%getlnvp()
  idproc = this%rf_re(0)%pp%getlidproc()

  do mode = 0, this%num_modes

    ub_re => b_re(mode)%get_f1()
    upsi_re => psi_re(mode)%get_f1()
    ue_re => this%rf_re(mode)%get_f1()
    if ( mode == 0 ) then
      do i = 1, nrp
        ue_re(1,i) = ub_re(2,i) - idrh * ( upsi_re(1,i+1) - upsi_re(1,i-1) )
        ue_re(2,i) = -ub_re(1,i)
      enddo
      if ( idproc == 0 ) then
        ue_re(1,1) = ub_re(2,1) - idrh * ( 4.0 * upsi_re(1,2) - upsi_re(1,3) - 3.0 * upsi_re(1,1) )
      endif
      if ( idproc == nvp-1 ) then
        ue_re(1,nrp) = ub_re(2,nrp) + idrh * ( 4.0 * upsi_re(1,nrp-1) - upsi_re(1,nrp-2) - 3.0 * upsi_re(1,nrp) )
      endif
      cycle
    endif

    ub_im => b_im(mode)%get_f1()
    upsi_im => psi_im(mode)%get_f1()
    ue_im => this%rf_im(mode)%get_f1()

    do i = 1, nrp
      k0 = real(i+noff) - 0.5
      ir = idr / k0
      ue_re(1,i) = ub_re(2,i) - idrh * ( upsi_re(1,i+1) - upsi_re(1,i-1) )
      ue_re(2,i) = -ub_re(1,i) + ir * mode * upsi_im(1,i)

      ue_im(1,i) = ub_im(2,i) - idrh * ( upsi_im(1,i+1) - upsi_im(1,i-1) )
      ue_im(2,i) = -ub_im(1,i) - ir * mode * upsi_re(1,i)
    enddo
    if ( idproc == 0 ) then
      ir = 2.0 * idr
      ue_re(1,1) = ub_re(2,1) - idrh * ( 4.0 * upsi_re(1,2) - upsi_re(1,3) - 3.0 * upsi_re(1,1) )
      ue_im(1,1) = ub_im(2,1) - idrh * ( 4.0 * upsi_im(1,2) - upsi_im(1,3) - 3.0 * upsi_im(1,1) )
    endif
    if ( idproc == nvp-1 ) then
      ir = idr / (nrp+noff-0.5)
      ue_re(1,nrp) = ub_re(2,nrp) + idrh * ( 4.0 * upsi_re(1,nrp-1) - upsi_re(1,nrp-2) - 3.0 * upsi_re(1,nrp) )
      ue_im(1,nrp) = ub_im(2,nrp) + idrh * ( 4.0 * upsi_im(1,nrp-1) - upsi_im(1,nrp-2) - 3.0 * upsi_im(1,nrp) )
    endif

  enddo

end subroutine solve_field_eperp

subroutine solve_field_eperp_beam( this, b )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_b ), intent(in) :: b

  type( ufield ), dimension(:), pointer :: b_re => null(), b_im => null()
  real, dimension(:,:), pointer :: ub_re => null(), ub_im => null()
  real, dimension(:,:), pointer :: ue_re => null(), ue_im => null()
  integer :: mode, i, nrp
  character(len=20), save :: sname = 'solve_field_eperp_beam'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  b_re => b%get_rf_re()
  b_im => b%get_rf_im()
  nrp = this%rf_re(0)%get_ndp(1)

  do mode = 0, this%num_modes

    ub_re => b_re(mode)%get_f1()
    ue_re => this%rf_re(mode)%get_f1()

    do i = 1, nrp
      ue_re(1,i) = ub_re(2,i)
      ue_re(2,i) = -ub_re(1,i)
    enddo

    if ( mode == 0 ) cycle

    ub_im => b_im(mode)%get_f1()
    ue_im => this%rf_im(mode)%get_f1()

    do i = 1, nrp
      ue_im(1,i) = ub_im(2,i)
      ue_im(2,i) = -ub_im(1,i)
    enddo

  enddo

end subroutine solve_field_eperp_beam

end module field_e_class