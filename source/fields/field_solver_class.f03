module field_solver_class

use grid_class
use parallel_pipe_class

use mpi
use param
use sys
use debug_tool

implicit none

private

public :: field_solver
public :: HYPRE_BUF

character(len=20), parameter :: cls_name = "field_solver"
integer, parameter :: cls_level = 4

integer, dimension(4), save :: itime
double precision, save :: dtime

real, dimension(:), pointer, save :: HYPRE_BUF => null()

type :: field_solver ! class for HYPRE solver

  ! HYPRE parameters
  integer, dimension(:), pointer :: offsets => null()
  integer, dimension(:), pointer :: stencil_idx => null()
  integer :: num_stencil
  integer :: stype
  integer :: kind
  integer :: mode
  integer :: bnd
  real :: tol

  integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver, precond
  integer(HYPRE_TYPE) :: par_A, par_b, par_x
  integer :: iupper, ilower

  contains

  generic :: new => init_field_solver
  generic :: solve => solve_equation
  generic :: del => end_field_solver

  procedure, private :: init_field_solver, end_field_solver
  procedure, private :: set_struct_solver
  procedure, private :: solve_equation
  procedure, private :: set_struct_grid
  procedure, private :: set_struct_stencil
  procedure, private :: set_struct_matrix
  ! procedure, private :: set_ij_matrix
  ! procedure, private :: set_ij_solver

end type field_solver


contains

! =====================================================================
! Class field_solver implementation
! =====================================================================

subroutine init_field_solver( this, pp, gp, mode, dr, kind, bnd, stype, tol )

  implicit none

  class( field_solver ), intent(inout) :: this
  class( parallel_pipe ), intent(in) :: pp
  class( grid ), intent(in) :: gp
  integer, intent(in) :: kind, stype, mode, bnd
  real, intent(in) :: dr, tol

  integer :: i, j, ierr, comm
  character(len=20), save :: sname = "init_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%stype = stype
  this%tol   = tol
  this%mode  = mode
  this%kind  = kind
  this%bnd   = bnd

  ! setup HYPRE grid
  comm = pp%getlgrp()

  ! select case ( this%kind )

  ! case ( p_fk_psi, p_fk_ez, p_fk_bz, p_fk_bt, p_fk_bplus, p_fk_bminus )

  !   call this%set_struct_grid( pp, gp )
  !   call this%set_struct_stencil()
  !   call this%set_struct_matrix( pp, gp, dr )

  !   call HYPRE_StructVectorCreate( comm, this%grid, this%b, ierr )
  !   call HYPRE_StructVectorInitialize( this%b, ierr )

  !   call HYPRE_StructVectorCreate( comm, this%grid, this%x, ierr )
  !   call HYPRE_StructVectorInitialize( this%x, ierr )

  !   call this%set_struct_solver( pp )

  ! case ( p_fk_bt_old, p_fk_bt_iter )

  !   call this%set_ij_matrix( pp, gp, dr )

  !   call HYPRE_IJVectorCreate( comm, this%ilower, this%iupper, this%b, ierr )
  !   call HYPRE_IJVectorSetObjectType( this%b, HYPRE_PARCSR, ierr )
  !   call HYPRE_IJVectorInitialize( this%b, ierr )

  !   call HYPRE_IJVectorCreate( comm, this%ilower, this%iupper, this%x, ierr )
  !   call HYPRE_IJVectorSetObjectType( this%x, HYPRE_PARCSR, ierr )
  !   call HYPRE_IJVectorInitialize( this%x, ierr )

  !   call this%set_ij_solver( pp )

  ! end select

  call this%set_struct_grid( pp, gp )
  call this%set_struct_stencil()
  call this%set_struct_matrix( pp, gp, dr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%b, ierr )
  call HYPRE_StructVectorInitialize( this%b, ierr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%x, ierr )
  call HYPRE_StructVectorInitialize( this%x, ierr )

  call this%set_struct_solver( pp )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_solver

subroutine end_field_solver( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: ierr
  character(len=20), save :: sname = "end_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( associated(HYPRE_BUF) ) deallocate( HYPRE_BUF )

  ! select case ( this%kind )
  ! case ( p_fk_psi, p_fk_bz, p_fk_ez, p_fk_bt, p_fk_bplus, p_fk_bminus )
  !   call HYPRE_StructGridDestroy( this%grid, ierr )
  !   call HYPRE_StructStencilDestroy( this%stencil, ierr )
  !   call HYPRE_StructVectorDestroy( this%b, ierr )
  !   call HYPRE_StructVectorDestroy( this%x, ierr )
  !   call HYPRE_StructMatrixDestroy( this%A, ierr )
  ! case ( p_fk_bt_iter, p_fk_bt_old )
  !   call HYPRE_IJMatrixDestroy( this%A, ierr )
  !   call HYPRE_IJVectorDestroy( this%b, ierr )
  !   call HYPRE_IJVectorDestroy( this%x, ierr )
  ! end select

  call HYPRE_StructGridDestroy( this%grid, ierr )
  call HYPRE_StructStencilDestroy( this%stencil, ierr )
  call HYPRE_StructVectorDestroy( this%b, ierr )
  call HYPRE_StructVectorDestroy( this%x, ierr )
  call HYPRE_StructMatrixDestroy( this%A, ierr )

  select case ( this%stype )
  case ( p_hypre_cycred )
    call HYPRE_StructCycRedDestroy( this%solver, ierr )
  case ( p_hypre_smg )
    call HYPRE_StructSMGDestroy( this%solver, ierr )
  case ( p_hypre_pcg )
    call HYPRE_StructPCGDestroy( this%solver, ierr )
    call HYPRE_StructSMGDestroy( this%precond, ierr )
  case ( p_hypre_gmres )
    call HYPRE_StructGMRESDestroy( this%solver, ierr )
    call HYPRE_StructJacobiDestroy( this%precond, ierr )
  case ( p_hypre_amg )
    call HYPRE_BoomerAMGDestroy( this%solver, ierr )
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_solver

subroutine set_struct_solver( this, pp )

  implicit none

  class( field_solver ), intent(inout) :: this
  class( parallel_pipe ), intent(in) :: pp

  integer :: n_post = 1, n_pre = 1, maxiter = 1000, maxiter_pre = 10
  integer :: i, ierr, precond_id, comm
  character(len=20), save :: sname = "set_struct_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = pp%getlgrp()

  select case ( this%stype )

  case ( p_hypre_cycred )

    call write_stdout( 'mode '//num2str(this%mode)//': Using Cyclic Reduction solver' )
    call HYPRE_StructCycRedCreate( comm, this%solver, ierr )
    call HYPRE_StructCycRedSetup( this%solver, this%A, this%b, this%x, ierr )

  case ( p_hypre_smg )

    call write_stdout( 'mode '//num2str(this%mode)//': Using SMG solver' )
    call HYPRE_StructSMGCreate( comm, this%solver, ierr )
    call HYPRE_StructSMGSetMemoryUse( this%solver, 0, ierr )
    call HYPRE_StructSMGSetMaxIter( this%solver, maxiter, ierr )
    call HYPRE_StructSMGSetTol( this%solver, this%tol, ierr )
    call HYPRE_StructSMGSetPrintLevel( this%solver, 1, ierr )
    call HYPRE_StructSMGSetRelChange( this%solver, 0, ierr )
    call HYPRE_StructSMGSetNumPreRelax( this%solver, n_pre, ierr )
    call HYPRE_StructSMGSetNumPostRelax( this%solver, n_post, ierr )
    call HYPRE_StructSMGSetLogging( this%solver, 1, ierr )
    call HYPRE_StructSMGSetup( this%solver, this%A, this%b, this%x, ierr )

  case ( p_hypre_pcg )

    call write_stdout( 'mode '//num2str(this%mode)//': Using PCG solver' )
    call HYPRE_StructPCGCreate( comm, this%solver, ierr )
    call HYPRE_StructPCGSetMaxIter( this%solver, maxiter, ierr )
    call HYPRE_StructPCGSetTol( this%solver, this%tol, ierr )
    call HYPRE_StructPCGSetPrintLevel( this%solver, 2, ierr )
    call HYPRE_StructPCGSetTwoNorm( this%solver, 1, ierr )
    call HYPRE_StructPCGSetRelChange( this%solver, 0, ierr )
    call HYPRE_StructPCGSetLogging( this%solver, 1, ierr )

    ! set up preconditioner
    precond_id = 7
    call HYPRE_StructJacobiCreate( comm, this%precond, ierr )
    ! call HYPRE_StructSMGSetMemoryUse( this%precond, 0, ierr )
    call HYPRE_StructJacobiSetMaxIter( this%precond, 20, ierr )
    call HYPRE_StructJacobiSetTol( this%precond, this%tol, ierr )
    ! call HYPRE_StructSMGSetNumPreRelax( this%precond, n_pre, ierr )
    ! call HYPRE_StructSMGSetNumPostRelax( this%precond, n_post, ierr )
    ! call HYPRE_StructSMGSetLogging( this%precond, 0, ierr )

    call HYPRE_StructPCGSetPrecond( this%solver, precond_id, this%precond, ierr )
    call HYPRE_StructPCGSetup( this%solver, this%A, this%b, this%x, ierr )

  case ( p_hypre_gmres )

    call write_stdout( 'mode '//num2str(this%mode)//': Using GMRES solver' )
    call HYPRE_StructGMRESCreate( comm, this%solver, ierr )
    call HYPRE_StructGMRESSetMaxIter( this%solver, maxiter, ierr )
    call HYPRE_StructGMRESSetTol( this%solver, this%tol, ierr )
    call HYPRE_StructGMRESSetPrintLevel( this%solver, 2, ierr )
    call HYPRE_StructGMRESSetLogging( this%solver, 1, ierr )

    call HYPRE_StructJacobiCreate( comm, this%precond, ierr )
    call HYPRE_StructJacobiSetMaxIter( this%precond, 20, ierr )
    call HYPRE_StructJacobiSetTol( this%precond, 0.0, ierr )
    call HYPRE_StructJacobiSetZeroGuess( this%precond, ierr )

    precond_id = 7
    call HYPRE_StructGMRESSetPrecond( this%solver, precond_id, this%precond, ierr )
    call HYPRE_StructGMRESSetup( this%solver, this%A, this%b, this%x, ierr )

  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_solver

subroutine solve_equation( this, src_sol )

  implicit none

  class( field_solver ), intent(inout) :: this
  real, intent(inout), dimension(:), pointer :: src_sol

  integer :: local_size, i, ierr
  integer, dimension(:), pointer, save :: rows => null()
  character(len=20), save :: sname = "solve_equation"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( this%kind )
  case ( p_fk_psi )
    call start_tprof( 'solve psi' )
  case ( p_fk_ez )
    call start_tprof( 'solve ez' )
  case ( p_fk_bz )
    call start_tprof( 'solve bz' )
  case ( p_fk_bt )
    call start_tprof( 'solve beam bt' )
  case ( p_fk_bplus, p_fk_bminus )
    call start_tprof( 'solve plasma bt' )
  end select

  ! select case ( this%kind )
  ! case ( p_fk_psi, p_fk_ez, p_fk_bz, p_fk_bt, p_fk_bplus, p_fk_bminus )
  !   call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
  !   call HYPRE_StructVectorAssemble( this%b, ierr )
  ! case ( p_fk_bt_old, p_fk_bt_iter )
  !   local_size = this%iupper - this%ilower + 1
  !   if ( .not. associated(rows) ) then
  !     allocate( rows(local_size) )
  !     do i = this%ilower, this%iupper
  !       rows(i-this%ilower+1) = i
  !     enddo
  !   endif
  !   call HYPRE_IJVectorSetValues( this%b, local_size, rows, src_sol, ierr )
  !   call HYPRE_IJVectorAssemble( this%b, ierr )
  !   call HYPRE_IJVectorGetObject( this%b, this%par_b, ierr )
  !   call HYPRE_IJVectorAssemble( this%x, ierr )
  !   call HYPRE_IJVectorGetObject( this%x, this%par_x, ierr )
  ! end select

  call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
  call HYPRE_StructVectorAssemble( this%b, ierr )

  select case ( this%stype )
  case ( p_hypre_cycred )
    call HYPRE_StructCycRedSolve( this%solver, this%A, this%b, this%x, ierr )
  case ( p_hypre_smg )
    call HYPRE_StructSMGSolve( this%solver, this%A, this%b, this%x, ierr )
  case ( p_hypre_pcg )
    call HYPRE_StructPCGSolve( this%solver, this%A, this%b, this%x, ierr )
  case ( p_hypre_gmres )
    call HYPRE_StructGMRESSolve( this%solver, this%A, this%b, this%x, ierr )
  case ( p_hypre_amg )
    call HYPRE_BoomerAMGSolve( this%solver, this%par_A, this%par_b, this%par_x, ierr )
  end select

  ! select case ( this%kind )
  ! case ( p_fk_psi, p_fk_ez, p_fk_bz, p_fk_bt, p_fk_bplus, p_fk_bminus )
  !   call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )
  ! case ( p_fk_bt_old, p_fk_bt_iter )
  !   call HYPRE_IJVectorGetValues( this%x, local_size, rows, src_sol, ierr )
  ! end select

  call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )

  select case ( this%kind )
  case ( p_fk_psi )
    call stop_tprof( 'solve psi' )
  case ( p_fk_ez )
    call stop_tprof( 'solve ez' )
  case ( p_fk_bz )
    call stop_tprof( 'solve bz' )
  case ( p_fk_bt )
    call stop_tprof( 'solve beam bt' )
  case ( p_fk_bplus, p_fk_bminus )
    call stop_tprof( 'solve plasma bt' )
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_equation

subroutine set_struct_grid( this, pp, gp )

  implicit none

  class( field_solver ), intent(inout) :: this
  class( parallel_pipe ), intent(in) :: pp
  class( grid ), intent(in) :: gp

  integer :: comm, ierr, dim = 1
  character(len=20), save :: sname = "set_struct_grid"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = pp%getlgrp()
  this%ilower = gp%get_noff(1) + 1
  this%iupper = gp%get_noff(1) + gp%get_ndp(1)

  call HYPRE_StructGridCreate( comm, 1, this%grid, ierr )
  call HYPRE_StructGridSetExtents( this%grid, this%ilower, this%iupper, ierr )
  call HYPRE_StructGridAssemble( this%grid, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_grid

subroutine set_struct_stencil( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: i, ierr
  character(len=20), save :: sname = "set_struct_stencil"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%num_stencil = 3
  if ( .not. associated( this%offsets ) ) then
    allocate( this%offsets( this%num_stencil ) )
  endif

  if ( .not. associated( this%stencil_idx ) ) then
    allocate( this%stencil_idx( this%num_stencil ) )
    do i = 1, this%num_stencil
      this%stencil_idx(i) = i-1
    enddo
  endif

  this%offsets = (/-1, 0, 1/)

  call HYPRE_StructStencilCreate( 1, this%num_stencil, this%stencil, ierr )
  do i = 1, this%num_stencil
    call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(i), ierr )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_stencil

subroutine set_struct_matrix( this, pp, gp, dr )

  implicit none

  class( field_solver ), intent(inout) :: this
  class( parallel_pipe ), intent(in) :: pp
  class( grid ), intent(in) :: gp
  real, intent(in) :: dr

  integer :: i, ierr, local_vol, nr, noff, m
  integer :: comm, lidproc, lnvp
  real :: dr2, m2, j, jmax
  character(len=20), save :: sname = "set_struct_matrix"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = pp%getlgrp()
  lidproc = pp%getlidproc()
  lnvp = pp%getlnvp()
  noff = gp%get_noff(1)

  dr2 = dr*dr
  m = this%mode
  m2 = real(m*m)
  nr = this%iupper - this%ilower + 1
  local_vol = nr * this%num_stencil

  if ( .not. associated( HYPRE_BUF ) ) then
    allocate( HYPRE_BUF( local_vol ) )
  elseif ( size(HYPRE_BUF) < local_vol ) then
    deallocate( HYPRE_BUF )
    allocate( HYPRE_BUF( local_vol ) )
  endif

  call HYPRE_StructMatrixCreate( comm, this%grid, this%stencil, this%A, ierr )
  call HYPRE_StructMatrixInitialize( this%A, ierr )

  ! set the matrix element and lower boundary
  select case ( this%kind )

  case ( p_fk_psi, p_fk_bt, p_fk_ez, p_fk_bz )

    ! set from the second grid point of each partition
    j = real(noff)
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
      HYPRE_BUF(i+1) = -2.0 - m2 / j**2
      HYPRE_BUF(i+2) = 1.0 + 0.5 / j
    enddo

    ! set the first grid point of each partition
    if (lidproc == 0) then

      if (m == 0) then
        HYPRE_BUF(1) = 0.0
        HYPRE_BUF(2) = -4.0
        HYPRE_BUF(3) = 4.0
      else
        ! matrix elements 1 to 3 are given arbitrarily to make sure the matrix
        ! is not singular. The vanishing of element 4 indicates the on-axis field
        ! value is zero.
        HYPRE_BUF(1) = 0.0
        HYPRE_BUF(2) = 1.0
        HYPRE_BUF(3) = 0.0
        HYPRE_BUF(4) = 0.0
      endif

    else

      j = real(noff)
      HYPRE_BUF(1) = 1.0 - 0.5 / j
      HYPRE_BUF(2) = -2.0 - m2 / j**2
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  case ( p_fk_bplus )

    ! set from the second grid point of each partition
    j = real(noff)
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
      HYPRE_BUF(i+1) = -2.0 - ((m+1)/j)**2 ) - dr2
      HYPRE_BUF(i+2) = 1.0 + 0.5 / j
    enddo

    ! set the first grid point of each partition
    if (lidproc == 0) then

      ! matrix elements 1 to 3 are given arbitrarily to make sure the matrix
      ! is not singular. The vanishing of element 4 indicates the on-axis field
      ! value is zero.
      HYPRE_BUF(1) = 0.0
      HYPRE_BUF(2) = 1.0
      HYPRE_BUF(3) = 0.0
      HYPRE_BUF(4) = 0.0

    else

      j = real(noff)
      HYPRE_BUF(1) = 1.0 - 0.5 / j
      HYPRE_BUF(2) = -2.0 - ((m+1)/j)**2 ) - dr2
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  case ( p_fk_bminus )

    ! set from the second grid point of each partition
    j = real(noff)
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
      HYPRE_BUF(i+1) = -2.0 - ((m-1)/j)**2 ) - dr2
      HYPRE_BUF(i+2) = 1.0 + 0.5 / j
    enddo

    ! set the first grid point of each partition
    if (lidproc == 0) then

      if (m == 1) then
        HYPRE_BUF(1) = 0.0
        HYPRE_BUF(2) = -4.0 - dr2
        HYPRE_BUF(3) = 4.0
      else
        ! matrix elements 1 to 3 are given arbitrarily to make sure the matrix
        ! is not singular. The vanishing of element 4 indicates the on-axis field
        ! value is zero.
        HYPRE_BUF(1) = 0.0
        HYPRE_BUF(2) = 1.0
        HYPRE_BUF(3) = 0.0
        HYPRE_BUF(4) = 0.0
      endif

    else

      j = real(noff)
      HYPRE_BUF(1) = 1.0 - 0.5 / j
      HYPRE_BUF(2) = -2.0 - ((m-1)/j)**2 ) - dr*dr
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  end select

  ! set the upper boundary
  if ( lidproc == lnvp-1 ) then

    select case ( this%bnd )

    case ( p_bnd_zero )

      HYPRE_BUF(local_vol) = 0.0

    case ( p_bnd_open )

      jmax = real(noff + nr)

      select case ( this%kind )

      case ( p_fk_psi )

        if ( this%mode == 0 ) then
          HYPRE_BUF(local_vol) = 0.0
        else
          HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-m/jmax) * HYPRE_BUF(local_vol)
          HYPRE_BUF(local_vol) = 0.0
        endif

      case ( p_fk_bt )

        if ( this%mode == 0 ) then
          HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0+1.0/(jmax*log(jmax*dr))) * HYPRE_BUF(local_vol)
          HYPRE_BUF(local_vol) = 0.0
        else
          HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-m/jmax) * HYPRE_BUF(local_vol)
          HYPRE_BUF(local_vol) = 0.0
        endif

      case ( p_fk_ez, p_fk_bz )

        if ( this%mode == 0 ) then
          HYPRE_BUF(local_vol) = 0.0
        else
          HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-m/jmax) * HYPRE_BUF(local_vol)
          HYPRE_BUF(local_vol) = 0.0
        endif

      case ( p_fk_bplus )

        HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-(m+1)/jmax) * HYPRE_BUF(local_vol)
        HYPRE_BUF(local_vol) = 0.0

      case ( p_fk_bminus )

        HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-(m+1)/jmax) * HYPRE_BUF(local_vol)
        HYPRE_BUF(local_vol) = 0.0

      case default
        call write_err( 'Invalid field type!' )
      end select

    case default
      call write_err( 'Invalid boundary condition!' )
    end select

  endif

  call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
    this%stencil_idx, HYPRE_BUF, ierr )

  call HYPRE_StructMatrixAssemble( this%A, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_matrix

! subroutine set_ij_solver( this, pp )

!   implicit none

!   class( field_solver ), intent(inout) :: this
!   class( parallel_pipe ), intent(in) :: pp

!   integer :: n_post = 1, n_pre = 1, maxiter = 1000, maxiter_pre = 10
!   integer :: i, ierr, precond_id, comm
!   character(len=20), save :: sname = "set_ij_solver"

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   comm = pp%getlgrp()

!   select case ( this%stype )

!   case ( p_hypre_amg )

!     call write_stdout( 'mode '//num2str(this%mode)//": Using Boomer AMG solver" )
!     call HYPRE_BoomerAMGCreate( this%solver, ierr )
!     call HYPRE_BoomerAMGSetMaxLevels( this%solver, 20, ierr )
!     call HYPRE_BoomerAMGSetMaxIter( this%solver, 100, ierr )
!     call HYPRE_BoomerAMGSetTol( this%solver, this%tol, ierr )
!     ! call HYPRE_BoomerAMGSetPrintLevel( this%solver, 3, ierr )
!     ! call HYPRE_BoomerAMGSetOldDefault( this%solver, ierr )
!     call HYPRE_BoomerAMGSetRelaxType( this%solver, 6, ierr )
!     call HYPRE_BoomerAMGSetRelaxOrder( this%solver, 1, ierr )
!     call HYPRE_BoomerAMGSetNumSweeps( this%solver, 10, ierr )
!     ! call HYPRE_BoomerAMGSetCycleNumSweeps( this%solver, 10, 3, ierr )
!     ! call HYPRE_BoomerAMGSetMaxCoarseSize( this%solver, 20 )
!     ! call HYPRE_BoomerAMGSetMinCoarseSize( this%solver, 5 )
!     call HYPRE_BoomerAMGSetup( this%solver, this%par_A, this%par_b, this%par_x, ierr )

!   case ( p_hypre_smg )

!     call write_stdout( 'mode '//num2str(this%mode)//': Using SMG solver' )
!     call HYPRE_StructSMGCreate( comm, this%solver, ierr )
!     call HYPRE_StructSMGSetMemoryUse( this%solver, 0, ierr )
!     call HYPRE_StructSMGSetMaxIter( this%solver, maxiter, ierr )
!     call HYPRE_StructSMGSetTol( this%solver, this%tol, ierr )
!     call HYPRE_StructSMGSetPrintLevel( this%solver, 1, ierr )
!     call HYPRE_StructSMGSetRelChange( this%solver, 0, ierr )
!     call HYPRE_StructSMGSetNumPreRelax( this%solver, n_pre, ierr )
!     call HYPRE_StructSMGSetNumPostRelax( this%solver, n_post, ierr )
!     call HYPRE_StructSMGSetLogging( this%solver, 1, ierr )
!     call HYPRE_StructSMGSetup( this%solver, this%A, this%b, this%x, ierr )

!   case ( p_hypre_parpcg )

!     call write_stdout( 'mode '//num2str(this%mode)//': Using ParCSR PCG solver' )
!     call HYPRE_ParCSRPCGCreate( comm, this%solver, ierr )
!     call HYPRE_ParCSRPCGSetMaxIter( this%solver, maxiter, ierr )
!     call HYPRE_ParCSRPCGSetTol( this%solver, this%tol, ierr)
!     call HYPRE_ParCSRPCGSetTwoNorm( this%solver, 1, ierr)
!     ! call HYPRE_ParCSRPCGSetPrintLevel(this%solver, 2, ierr)
!     ! call HYPRE_ParCSRPCGSetLogging(this%solver, 1, ierr)

!     call HYPRE_BoomerAMGCreate( this%precond, ierr )
!     ! call HYPRE_BoomerAMGSetPrintLevel(this%precond, 1, ierr)
!     call HYPRE_BoomerAMGSetCoarsenType( this%precond, 6, ierr )
!     ! call HYPRE_BoomerAMGSetRelaxType(this%precond, 6, ierr)
!     ! call HYPRE_BoomerAMGSetNumSweeps(this%precond, 1, ierr)
!     call HYPRE_BoomerAMGSetTol( this%precond, 0.0d0, ierr )
!     call HYPRE_BoomerAMGSetMaxIter( this%precond, 1, ierr )

!     precond_id = 2
!     call HYPRE_ParCSRPCGSetPrecond( this%solver, precond_id, this%precond, ierr )

!     call HYPRE_ParCSRPCGSetup( this%solver, this%par_A, this%par_b, this%par_x, ierr )

!   end select

!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine set_ij_solver

! subroutine set_ij_matrix( this, pp, gp, dr )

!   implicit none

!   class( field_solver ), intent(inout) :: this
!   class( parallel_pipe ), intent(in) :: pp
!   class( grid ), intent(in) :: gp
!   real, intent(in) :: dr

!   integer :: local_size, ierr, m, m2, i, nr, noff
!   integer :: comm
!   integer, dimension(:), pointer :: cols
!   real :: idr2, k0, k_minus, k_plus, rmax
!   character(len=20), save :: sname = "set_struct_matrix"

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   comm = pp%getlgrp()
!   nr = gp%get_nd(1)
!   noff = gp%get_noff(1)
!   this%ilower = 4*noff + 1
!   this%iupper = 4*noff + 4*gp%get_ndp(1)
!   rmax = ( noff + nr - 0.5 ) * dr

!   local_size = this%iupper - this%ilower + 1

!   call HYPRE_IJMatrixCreate( comm, this%ilower, this%iupper, &
!       this%ilower, this%iupper, this%A, ierr )
!   call HYPRE_IJMatrixSetObjectType( this%A, HYPRE_PARCSR, ierr )
!   call HYPRE_IJMatrixInitialize( this%A, ierr )

!   if ( .not. associated( HYPRE_BUF ) ) then
!     allocate( HYPRE_BUF(4) )
!   elseif ( size(HYPRE_BUF) < 4 ) then
!     deallocate( HYPRE_BUF )
!     allocate( HYPRE_BUF(4) )
!   endif
!   allocate( cols(4) )

!   m = this%mode
!   m2 = m*m
!   idr2 = 1.0 / dr**2
!   HYPRE_BUF = 0.0


!   do i = this%ilower, this%iupper, 4

!     k0 = i/4 + 0.5
!     k_minus = k0 - 0.5
!     k_plus  = k0 + 0.5

!     ! set Re(Br)
!     cols = (/i-4, i, i+3, i+4/)
!     HYPRE_BUF(1) = k_minus / k0 * idr2
!     HYPRE_BUF(2) = -idr2 * (2.0 + (m2+1.0) / k0**2)
!     HYPRE_BUF(3) = 2.0*m / k0**2 * idr2
!     HYPRE_BUF(4) = k_plus  / k0 * idr2
!     if ( this%kind == p_fk_bt_iter ) HYPRE_BUF(2) = HYPRE_BUF(2) - 1.0

!     if ( i == 1 ) then
!       call HYPRE_IJMatrixSetValues( this%A, 1, 3, i, cols(2), HYPRE_BUF(2), ierr )
!     elseif ( i == 4*nr-3 ) then
!       select case ( this%bnd )
!       case ( p_bnd_zero )
!         ! do nothing
!       case ( p_bnd_open )
!         HYPRE_BUF(2) = HYPRE_BUF(2) + (1-(m+1)*dr/rmax) * HYPRE_BUF(4)
!       end select
!       call HYPRE_IJMatrixSetValues( this%A, 1, 3, i, cols, HYPRE_BUF, ierr )
!     else
!       call HYPRE_IJMatrixSetValues( this%A, 1, 4, i, cols, HYPRE_BUF, ierr )
!     endif

!     ! set Im(Br)
!     cols = (/i-3, i+1, i+2, i+5/)

!     HYPRE_BUF(1) = k_minus / k0 * idr2
!     HYPRE_BUF(2) = -idr2 * (2.0 + (m2+1.0) / k0**2)
!     HYPRE_BUF(3) = -2.0*m / k0**2 * idr2
!     HYPRE_BUF(4) = k_plus  / k0 * idr2
!     if ( this%kind == p_fk_bt_iter ) HYPRE_BUF(2) = HYPRE_BUF(2) - 1.0

!     if ( i == 1 ) then
!       call HYPRE_IJMatrixSetValues( this%A, 1, 3, i+1, cols(2), HYPRE_BUF(2), ierr )
!     elseif ( i == 4*nr-3 ) then
!       select case ( this%bnd )
!       case ( p_bnd_zero )
!         ! do nothing
!       case ( p_bnd_open )
!         HYPRE_BUF(2) = HYPRE_BUF(2) + (1-(m+1)*dr/rmax) * HYPRE_BUF(4)
!       end select
!       call HYPRE_IJMatrixSetValues( this%A, 1, 3, i+1, cols, HYPRE_BUF, ierr )
!     else
!       call HYPRE_IJMatrixSetValues( this%A, 1, 4, i+1, cols, HYPRE_BUF, ierr )
!     endif

!     ! set Re(Bphi)
!     cols = (/i-2, i+1, i+2, i+6/)

!     HYPRE_BUF(1) = k_minus / k0 * idr2
!     HYPRE_BUF(2) = -2.0*m / k0**2 * idr2
!     HYPRE_BUF(3) = -idr2 * (2.0 + (m2+1.0) / k0**2)
!     HYPRE_BUF(4) = k_plus  / k0 * idr2
!     if ( this%kind == p_fk_bt_iter ) HYPRE_BUF(3) = HYPRE_BUF(3) - 1.0

!     if ( i == 1 ) then
!       call HYPRE_IJMatrixSetValues( this%A, 1, 3, i+2, cols(2), HYPRE_BUF(2), ierr )
!     elseif ( i == 4*nr-3 ) then
!       select case ( this%bnd )
!       case ( p_bnd_zero )
!         ! do nothing
!       case ( p_bnd_open )
!         HYPRE_BUF(3) = HYPRE_BUF(3) + (1-(m+1)*dr/rmax) * HYPRE_BUF(4)
!       end select
!       call HYPRE_IJMatrixSetValues( this%A, 1, 3, i+2, cols, HYPRE_BUF, ierr )
!     else
!       call HYPRE_IJMatrixSetValues( this%A, 1, 4, i+2, cols, HYPRE_BUF, ierr )
!     endif

!     ! set Im(Bphi)
!     cols = (/i-1, i, i+3, i+7/)

!     HYPRE_BUF(1) = k_minus / k0 * idr2
!     HYPRE_BUF(2) = 2.0*m / k0**2 * idr2
!     HYPRE_BUF(3) = -idr2 * (2.0 + (m2+1.0) / k0**2)
!     HYPRE_BUF(4) = k_plus  / k0 * idr2
!     if ( this%kind == p_fk_bt_iter ) HYPRE_BUF(3) = HYPRE_BUF(3) - 1.0

!     if ( i == 1 ) then
!       call HYPRE_IJMatrixSetValues( this%A, 1, 3, i+3, cols(2), HYPRE_BUF(2), ierr )
!     elseif ( i == 4*nr-3 ) then
!       select case ( this%bnd )
!       case ( p_bnd_zero )
!         ! do nothing
!       case ( p_bnd_open )
!         HYPRE_BUF(3) = HYPRE_BUF(3) + (1-(m+1)*dr/rmax) * HYPRE_BUF(4)
!       end select
!       call HYPRE_IJMatrixSetValues( this%A, 1, 3, i+3, cols, HYPRE_BUF, ierr )
!     else
!       call HYPRE_IJMatrixSetValues( this%A, 1, 4, i+3, cols, HYPRE_BUF, ierr )
!     endif

!   enddo

!   call HYPRE_IJMatrixAssemble( this%A, ierr )
!   call HYPRE_IJMatrixGetObject( this%A, this%par_A, ierr )

!   deallocate( cols )

!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine set_ij_matrix


end module field_solver_class