module field_solver_class

use mpi
use param

implicit none

private

public :: field_solver
public :: HYPRE_BUF

real, dimension(:), pointer :: HYPRE_BUF => null()

type :: field_solver ! class for HYPRE solver

  integer :: order

  ! HYPRE parameters
  integer, dimension(:), pointer :: offsets => null()
  integer, dimension(:), pointer :: stencil_idx => null()
  integer :: num_stencil
  integer :: solver_type
  integer :: precond_type
  integer :: kind
  integer :: mode
  real :: tol

  integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver, precond
  integer :: iupper, ilower
  ! real, dimension(:), pointer :: hypre_buf => null()

  contains

  generic :: new => init_field_solver
  generic :: solve => solve_hypre_equation
  generic :: del => end_field_solver

  procedure, private :: init_field_solver, end_field_solver
  procedure, private :: set_hypre_solver
  procedure, private :: solve_hypre_equation
  ! procedure(set_hypre_grid), private, deferred :: set_hypre_grid
  ! procedure(set_hypre_stencil), private, deferred :: set_hypre_stencil
  ! procedure(set_hypre_matrix), private, deferred :: set_hypre_matrix
  procedure, private :: set_hypre_grid
  procedure, private :: set_hypre_stencil
  procedure, private :: set_hypre_matrix

end type field_solver

! abstract interface
! subroutine set_hypre_grid( this, comm, noff, ndp )
!   import field_solver
!   implicit none
!   class( field_solver ), intent(inout) :: this
!   integer, intent(in) :: comm
!   integer, intent(in), dimension(2) :: noff, ndp
! end subroutine set_hypre_grid

! subroutine set_hypre_stencil( this )
!   import field_solver
!   implicit none
!   class( field_solver ), intent(inout) :: this
! end subroutine set_hypre_stencil

! subroutine set_hypre_matrix( this, comm, dr )
!   import field_solver
!   implicit none
!   class( field_solver ), intent(inout) :: this
!   integer, intent(in) :: comm
!   real, intent(in) :: dr
! end subroutine set_hypre_matrix
! end interface

! type, extends( field_solver ) :: field_solver_A

!   private

!   contains

!   procedure, private :: set_hypre_grid => set_hypre_grid_AB
!   procedure, private :: set_hypre_stencil => set_hypre_stencil_AB
!   procedure, private :: set_hypre_matrix => set_hypre_matrix_A

! end type field_solver_A

! type, extends( field_solver ) :: field_solver_B

!   private

!   contains

!   procedure, private :: set_hypre_grid => set_hypre_grid_AB
!   procedure, private :: set_hypre_stencil => set_hypre_stencil_AB
!   procedure, private :: set_hypre_matrix => set_hypre_matrix_B

! end type field_solver_B

contains

! =====================================================================
! Class field_solver implementation
! =====================================================================

subroutine init_field_solver( this, ndp, noff, order, mode, dr, solver_type, precond_type, tol )

  implicit none

  class( field_solver ), intent(inout) :: this
  integer, intent(in), dimension(2) :: ndp, noff
  integer, intent(in) :: order, solver_type, precond_type, mode
  real, intent(in) :: dr, tol

  integer :: i, j, ierr

  this%order = order
  this%solver_type = solver_type
  this%precond_type = precond_type
  this%tol = tol
  this%mode = mode

  ! setup HYPRE grid
  call this%set_hypre_grid( MPI_COMM_WORLD, noff, ndp )
  
  call this%set_hypre_stencil()
  
  call this%set_hypre_matrix( MPI_COMM_WORLD, dr )

  call HYPRE_StructVectorCreate( MPI_COMM_WORLD, this%grid, this%b, ierr )
  call HYPRE_StructVectorInitialize( this%b, ierr )

  call HYPRE_StructVectorCreate( MPI_COMM_WORLD, this%grid, this%x, ierr )
  call HYPRE_StructVectorInitialize( this%x, ierr )

  call this%set_hypre_solver( MPI_COMM_WORLD )

end subroutine init_field_solver

subroutine end_field_solver( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: ierr

  if ( associated(HYPRE_BUF) ) deallocate( HYPRE_BUF )

  call HYPRE_StructGridDestroy( this%grid, ierr )
  call HYPRE_StructStencilDestroy( this%stencil, ierr )
  call HYPRE_StructVectorDestroy( this%b, ierr )
  call HYPRE_StructVectorDestroy( this%x, ierr )
  call HYPRE_StructMatrixDestroy( this%A, ierr )

  select case ( this%solver_type )

  case ( p_hypre_smg )

    call HYPRE_StructSMGDestroy( this%solver, ierr )

  case ( p_hypre_pfmg )

    call HYPRE_StructPFMGDestroy( this%solver, ierr )

  case ( p_hypre_pcg )

    call HYPRE_StructPCGDestroy( this%solver, ierr )

    select case ( this%precond_type )
    case ( p_hypre_smg )
      call HYPRE_StructSMGDestroy( this%precond, ierr )
    case ( p_hypre_pfmg )
      call HYPRE_StructPFMGDestroy( this%precond, ierr )
    end select

  end select

end subroutine end_field_solver

subroutine set_hypre_solver( this, comm )

  implicit none

  class( field_solver ), intent(inout) :: this
  integer, intent(in) :: comm

  integer(HYPRE_TYPE) :: solver, precond, A, b, x
  integer :: n_post = 1, n_pre = 1, maxiter = 1000
  integer :: i, ierr, precond_id
  real :: tol

  b = this%b
  x = this%x
  tol = this%tol

  A = this%A
  solver = this%solver
  select case ( this%solver_type )

  case ( p_hypre_smg )

    print *, "Using SMG solver"
    call HYPRE_StructSMGCreate( comm, solver, ierr )
    call HYPRE_StructSMGSetMemoryUse( solver, 0, ierr )
    call HYPRE_StructSMGSetMaxIter( solver, maxiter, ierr )
    call HYPRE_StructSMGSetTol( solver, tol, ierr )
    ! call HYPRE_StructSMGSetPrintLevel( solver, 2, ierr )
    call HYPRE_StructSMGSetRelChange( solver, 0, ierr )
    call HYPRE_StructSMGSetNumPreRelax( solver, n_pre, ierr )
    call HYPRE_StructSMGSetNumPostRelax( solver, n_post, ierr )
    call HYPRE_StructSMGSetLogging( solver, 1, ierr )
    call HYPRE_StructSMGSetup( solver, A, b, x, ierr )

  case ( p_hypre_pfmg )

    print *, 'Using PFMG solver'
    call HYPRE_StructPFMGCreate( comm, solver, ierr )
    call HYPRE_StructPFMGSetMaxIter( solver, maxiter, ierr )
    call HYPRE_StructPFMGSetTol( solver, tol, ierr )
    ! call HYPRE_StructPFMGSetPrintLevel( solver, 2, ierr )
    call HYPRE_StructPFMGSetRelChange( solver, 0, ierr )
    ! weighted Jacobi = 1; red-black GS = 2
    call HYPRE_StructPFMGSetRelaxType( solver, 1, ierr )
    call HYPRE_StructPFMGSetNumPreRelax( solver, n_pre, ierr )
    call HYPRE_StructPFMGSetNumPostRelax( solver, n_post, ierr )
    call HYPRE_StructPFMGSetLogging( solver, 1, ierr )
    call HYPRE_StructPFMGSetup( solver, A, b, x, ierr )

  case ( p_hypre_pcg )

    print *, 'Using PCG solver'
    call HYPRE_StructPCGCreate( comm, solver, ierr )
    call HYPRE_StructPCGSetMaxIter( solver, maxiter, ierr )
    call HYPRE_StructPCGSetTol( solver, tol, ierr )
    ! call HYPRE_StructPCGSetPrintLevel( solver, 2, ierr )
    call HYPRE_StructPCGSetTwoNorm( solver, 1, ierr )
    call HYPRE_StructPCGSetRelChange( solver, 0, ierr )
    call HYPRE_StructPCGSetLogging( solver, 1, ierr )

    ! set up preconditioner
    precond = this%precond
    select case ( this%precond_type )

    case ( p_hypre_smg )

      precond_id = 0
      call HYPRE_StructSMGCreate( comm, precond, ierr )
      call HYPRE_StructSMGSetMemoryUse( precond, 0, ierr )
      call HYPRE_StructSMGSetMaxIter( precond, maxiter, ierr )
      call HYPRE_StructSMGSetTol( precond, tol, ierr )
      call HYPRE_StructSMGSetNumPreRelax( precond, n_pre, ierr )
      call HYPRE_StructSMGSetNumPostRelax( precond, n_post, ierr )
      call HYPRE_StructSMGSetLogging( precond, 0, ierr )

    case ( p_hypre_pfmg )

      precond_id = 1
      call HYPRE_StructPFMGCreate( comm, precond, ierr )
      call HYPRE_StructPFMGSetMaxIter( precond, maxiter, ierr )
      call HYPRE_StructPFMGSetTol( precond, tol, ierr )
      ! 1 - Jacobi, 2 - red-black GS
      call HYPRE_StructPFMGSetRelaxType( precond, 1, ierr )
      call HYPRE_StructPFMGSetNumPreRelax( precond, n_pre, ierr )
      call HYPRE_StructPFMGSetNumPostRelax( precond, n_post, ierr )
      call HYPRE_StructPFMGSetLogging( precond, 0, ierr )

    case default

      print *, "Preconditioner must be set up for PCG solver"
      stop

    end select
    call HYPRE_StructPCGSetPrecond( solver, precond_id, precond, ierr )

    call HYPRE_StructPCGSetup( solver, A, b, x, ierr )

  end select

end subroutine set_hypre_solver

subroutine solve_hypre_equation( this, src_sol )

  implicit none

  class( field_solver ), intent(inout) :: this
  real, intent(inout), dimension(:), pointer :: src_sol

  integer :: ierr

  call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )

  select case ( this%solver_type )
  case ( p_hypre_smg )
    call HYPRE_StructSMGSolve( this%solver, this%A, this%b, this%x, ierr )
  case ( p_hypre_pfmg )
    call HYPRE_StructPFMGSolve( this%solver, this%A, this%b, this%x, ierr )
  case ( p_hypre_pcg )
    call HYPRE_StructPCGSolve( this%solver, this%A, this%b, this%x, ierr )
  end select

  call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )

end subroutine solve_hypre_equation

subroutine set_hypre_grid( this, comm, noff, ndp )

  implicit none

  class( field_solver ), intent(inout) :: this
  integer, intent(in) :: comm
  integer, intent(in), dimension(2) :: noff, ndp

  integer :: ierr, dim = 1

  call HYPRE_StructGridCreate( comm, dim, this%grid, ierr )

  select case ( this%kind )

  case ( p_fk_psi, p_fk_ez, p_fk_bz, p_fk_br_iter, p_fk_bphi_iter )
  
    this%ilower = noff(1) + 1
    this%iupper = noff(1) + ndp(1)

  case ( p_fk_bperp )

    this%ilower = 4 * noff(1) + 1
    this%iupper = 4 * ( noff(1) + ndp(1) )

  case default

    print *, "Invalid field kind for solver."
    stop

  end select

  call HYPRE_StructGridSetExtents( this%grid, this%ilower, this%iupper, ierr )
  call HYPRE_StructGridAssemble( this%grid, ierr )

end subroutine set_hypre_grid

subroutine set_hypre_stencil( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: i, ierr

  select case ( this%order )
  
  case ( p_fs_2order )

    select case ( this%kind )

    case ( p_fk_psi, p_fk_ez, p_fk_bz, p_fk_br_iter, p_fk_bphi_iter )

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

    case ( p_fk_bperp )

      this%num_stencil = 7
      if ( .not. associated( this%offsets ) ) then
        allocate( this%offsets( this%num_stencil ) )
      endif

      if ( .not. associated( this%stencil_idx ) ) then
        allocate( this%stencil_idx( this%num_stencil ) )
        do i = 1, this%num_stencil
          this%stencil_idx(i) = i-1
        enddo
      endif

      this%offsets = (/-4, -3, -1, 0, 1, 3, 4/)

    case default

      print *, "Invalid field kind for solver."
      stop

    end select

  case ( p_fs_4order )

    print *, "HYPRE solver not implemented for 4order algorithm"
    stop

  case default

    print *, "HYPRE solver not implemented for specified algorithm"
    stop

  end select

  call HYPRE_StructStencilCreate( 1, this%num_stencil, this%stencil, ierr )
  do i = 1, this%num_stencil
    call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(i), ierr )
  enddo

end subroutine set_hypre_stencil

subroutine set_hypre_matrix( this, comm, dr )

  implicit none

  class( field_solver ), intent(inout) :: this
  integer, intent(in) :: comm
  real, intent(in) :: dr

  integer :: i, k, ierr, local_vol
  real :: idr2, m, m2, k1, k2

  idr2 = 1.0 / (dr*dr)
  m = this%mode
  m2 = m*m
  local_vol = (this%iupper - this%ilower + 1) * this%num_stencil

  if ( .not. associated( HYPRE_BUF ) ) then
    allocate( HYPRE_BUF( local_vol ) )
  elseif ( size(HYPRE_BUF) < local_vol ) then
    deallocate( HYPRE_BUF )
    allocate( HYPRE_BUF( local_vol ) )
  endif

  call HYPRE_StructMatrixCreate( comm, this%grid, this%stencil, this%A, ierr )
  call HYPRE_StructMatrixInitialize( this%A, ierr )

  k = 0
  select case ( this%order )

  case ( p_fs_2order )

    select case ( this%kind )

    case ( p_fk_psi, p_fk_ez, p_fk_bz )
      
      do i = 1, local_vol, this%num_stencil
        k = k + 1
        k1 = k - 0.5
        k2 = k - 1.0
        HYPRE_BUF(i) = k2 / k1 * idr2
        HYPRE_BUF(i+1) = -1.0 * (2.0 + m2 / k1**2) * idr2
        HYPRE_BUF(i+2) = k / k1 * idr2
      enddo

      HYPRE_BUF(1) = 0.0
      HYPRE_BUF(local_vol) = 0.0
      if ( this%kind == p_fk_bz ) then
        HYPRE_BUF(local_vol-2) = 2.0 * idr2
      endif

    case ( p_fk_br_iter, p_fk_bphi_iter )

      do i = 1, local_vol, this%num_stencil
        k = k + 1
        k1 = k - 0.5
        k2 = k - 1.0
        HYPRE_BUF(i) = k2 / k1 * idr2
        HYPRE_BUF(i+1) = -1.0 * (2.0 + (m2+1.0) / k1**2) * idr2 - 1.0
        HYPRE_BUF(i+2) = k / k1 * idr2
      enddo

      HYPRE_BUF(1) = 0.0
      HYPRE_BUF(local_vol) = 0.0
      if ( this%kind == p_fk_bphi_iter ) then
        HYPRE_BUF(local_vol-2) = 2.0 * idr2
      endif

    case ( p_fk_bperp )

      do i = 1, local_vol, this%num_stencil*4
        k = k + 1
        k1 = k - 0.5
        k2 = k - 1.0
        ! stencil for Re(Br)
        HYPRE_BUF(i)    = k2 / k1 * idr2
        HYPRE_BUF(i+1)  = 0.0
        HYPRE_BUF(i+2)  = 0.0
        HYPRE_BUF(i+3)  = -1.0 * (2.0 + (m2+1.0) / k1**2) * idr2
        HYPRE_BUF(i+4)  = 0.0
        HYPRE_BUF(i+5)  = 2.0*m / k1**2 * idr2
        HYPRE_BUF(i+6)  = k / k1 * idr2

        ! stencil for Im(Br)
        HYPRE_BUF(i+7)   = k2 / k1 * idr2
        HYPRE_BUF(i+8)   = 0.0
        HYPRE_BUF(i+9)   = 0.0
        HYPRE_BUF(i+10)  = -1.0 * (2.0 + (m2+1.0) / k1**2) * idr2
        HYPRE_BUF(i+11)  = -2.0*m / k1**2 * idr2
        HYPRE_BUF(i+12)  = 0.0
        HYPRE_BUF(i+13)  = k / k1 * idr2

        ! stencil for Re(Bphi)
        HYPRE_BUF(i+14)  = k2 / k1 * idr2
        HYPRE_BUF(i+15)  = 0.0
        HYPRE_BUF(i+16)  = -2.0*m / k1**2 * idr2
        HYPRE_BUF(i+17)  = -1.0 * (2.0 + (m2+1.0) / k1**2) * idr2
        HYPRE_BUF(i+18)  = 0.0
        HYPRE_BUF(i+19)  = 0.0
        HYPRE_BUF(i+20)  = k / k1 * idr2

        ! stencil for Im(Bphi)
        HYPRE_BUF(i+21)  = k2 / k1 * idr2
        HYPRE_BUF(i+22)  = 2.0*m / k1**2 * idr2
        HYPRE_BUF(i+23)  = 0.0
        HYPRE_BUF(i+24)  = -1.0 * (2.0 + (m2+1.0) / k1**2) * idr2
        HYPRE_BUF(i+25)  = 0.0
        HYPRE_BUF(i+26)  = 0.0
        HYPRE_BUF(i+27)  = k / k1 * idr2
      enddo

    end select

  case ( p_fs_4order )

    print *, "HYPRE solver not implemented for 4order algorithm"
    stop

  case default

    print *, "HYPRE solver not implemented for specified algorithm"
    stop

  end select

  call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
    this%stencil_idx, HYPRE_BUF, ierr )
  call HYPRE_StructMatrixAssemble( this%A, ierr )

end subroutine set_hypre_matrix

end module field_solver_class