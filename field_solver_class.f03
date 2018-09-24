module field_solver_class

use mpi
use param
use system

implicit none

private

public :: field_solver
public :: HYPRE_BUF

character(len=20), parameter :: cls_name = "field_solver"
integer, parameter :: cls_level = 2

real, dimension(:), pointer :: HYPRE_BUF => null()

type :: field_solver ! class for HYPRE solver

  integer :: order

  ! HYPRE parameters
  integer, dimension(:), pointer :: offsets => null()
  integer, dimension(:), pointer :: stencil_idx => null()
  integer :: num_stencil
  integer :: solver_type
  integer :: kind
  integer :: mode
  real :: tol

  integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver
  integer :: iupper, ilower
  ! real, dimension(:), pointer :: hypre_buf => null()

  contains

  generic :: new => init_field_solver
  generic :: solve => solve_hypre_equation
  generic :: del => end_field_solver

  procedure, private :: init_field_solver, end_field_solver
  procedure, private :: set_hypre_solver
  procedure, private :: solve_hypre_equation
  procedure, private :: set_hypre_grid
  procedure, private :: set_hypre_stencil
  procedure, private :: set_hypre_matrix

end type field_solver


contains

! =====================================================================
! Class field_solver implementation
! =====================================================================

subroutine init_field_solver( this, ndp, noff, order, kind, mode, dr, solver_type, tol )

  implicit none

  class( field_solver ), intent(inout) :: this
  integer, intent(in), dimension(2) :: ndp, noff
  integer, intent(in) :: order, kind, solver_type, mode
  real, intent(in) :: dr, tol

  integer :: i, j, ierr, comm
  character(len=20), save :: sname = "init_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%order = order
  this%solver_type = solver_type
  this%tol = tol
  this%mode = mode
  this%kind = kind

  ! setup HYPRE grid
  comm = MPI_COMM_WORLD
  call this%set_hypre_grid( comm, noff, ndp )
  
  call this%set_hypre_stencil()
  
  call this%set_hypre_matrix( comm, dr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%b, ierr )
  call HYPRE_StructVectorInitialize( this%b, ierr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%x, ierr )
  call HYPRE_StructVectorInitialize( this%x, ierr )

  call this%set_hypre_solver( comm )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_solver

subroutine end_field_solver( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: ierr
  character(len=20), save :: sname = "end_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( associated(HYPRE_BUF) ) deallocate( HYPRE_BUF )

  call HYPRE_StructGridDestroy( this%grid, ierr )
  call HYPRE_StructStencilDestroy( this%stencil, ierr )
  call HYPRE_StructVectorDestroy( this%b, ierr )
  call HYPRE_StructVectorDestroy( this%x, ierr )
  call HYPRE_StructMatrixDestroy( this%A, ierr )

  select case ( this%solver_type )
  case ( p_hypre_cycred )
    call HYPRE_StructCycRedDestroy( this%solver, ierr )
  case ( p_hypre_pcg )
    call HYPRE_StructPCGDestroy( this%solver, ierr )
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_solver

subroutine set_hypre_solver( this, comm )

  implicit none

  class( field_solver ), intent(inout) :: this
  integer, intent(in) :: comm

  integer :: n_post = 1, n_pre = 1, maxiter = 1000
  integer :: i, ierr
  character(len=20), save :: sname = "set_hypre_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( this%solver_type )

  case ( p_hypre_cycred )

    print *, "Using Cyclic Reduction solver"
    call HYPRE_StructCycRedCreate( comm, this%solver, ierr )
    call HYPRE_StructCycRedSetup( this%solver, this%A, this%b, this%x, ierr )

  case ( p_hypre_pcg )

    print *, 'Using PCG solver'
    call HYPRE_StructPCGCreate( comm, this%solver, ierr )
    call HYPRE_StructPCGSetMaxIter( this%solver, maxiter, ierr )
    call HYPRE_StructPCGSetTol( this%solver, this%tol, ierr )
    ! call HYPRE_StructPCGSetPrintLevel( solver, 2, ierr )
    call HYPRE_StructPCGSetTwoNorm( this%solver, 1, ierr )
    call HYPRE_StructPCGSetRelChange( this%solver, 0, ierr )
    call HYPRE_StructPCGSetLogging( this%solver, 1, ierr )
    call HYPRE_StructPCGSetup( this%solver, this%A, this%b, this%x, ierr )

  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_hypre_solver

subroutine solve_hypre_equation( this, src_sol )

  implicit none

  class( field_solver ), intent(inout) :: this
  real, intent(inout), dimension(:), pointer :: src_sol

  integer :: ierr
  character(len=20), save :: sname = "set_hypre_equation"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )

  select case ( this%solver_type )
  case ( p_hypre_cycred )
    call HYPRE_StructCycRedSolve( this%solver, this%A, this%b, this%x, ierr )
  case ( p_hypre_pcg )
    call HYPRE_StructPCGSolve( this%solver, this%A, this%b, this%x, ierr )
  end select

  call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )

  ! print *, "vector x = ", src_sol(this%ilower:this%iupper)

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_hypre_equation

subroutine set_hypre_grid( this, comm, noff, ndp )

  implicit none

  class( field_solver ), intent(inout) :: this
  integer, intent(in) :: comm
  integer, intent(in), dimension(2) :: noff, ndp

  integer :: ierr, dim = 1
  character(len=20), save :: sname = "set_hypre_grid"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

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

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_hypre_grid

subroutine set_hypre_stencil( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: i, ierr
  character(len=20), save :: sname = "set_hypre_stencil"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

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

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_hypre_stencil

subroutine set_hypre_matrix( this, comm, dr )

  implicit none

  class( field_solver ), intent(inout) :: this
  integer, intent(in) :: comm
  real, intent(in) :: dr

  integer :: i, k, ierr, local_vol
  real :: idr2, m, m2, k1, k2
  character(len=20), save :: sname = "set_hypre_matrix"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

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

      ! lower boundary
      HYPRE_BUF(1) = 0.0

      ! upper boundary
      if ( this%kind == p_fk_ez .or. &
        ( this%mode > 0 .and. this%kind == p_fk_psi ) ) then
        HYPRE_BUF(local_vol) = 0.0
      endif
      if ( this%kind == p_fk_bz .or. &
        ( this%mode == 0 .and. this%kind == p_fk_psi ) ) then

        HYPRE_BUF(local_vol) = 0.0
        HYPRE_BUF(local_vol-1) = 0.0 ! this is for eliminating the ambiguity of pure Neumann boundary
        HYPRE_BUF(local_vol-2) = 2.0 * idr2
      endif

    ! case ( p_fk_br_iter, p_fk_bphi_iter )

    !   do i = 1, local_vol, this%num_stencil
    !     k = k + 1
    !     k1 = k - 0.5
    !     k2 = k - 1.0
    !     HYPRE_BUF(i) = k2 / k1 * idr2
    !     HYPRE_BUF(i+1) = -1.0 * (2.0 + (m2+1.0) / k1**2) * idr2 - 1.0
    !     HYPRE_BUF(i+2) = k / k1 * idr2
    !   enddo

    !   HYPRE_BUF(1) = 0.0
    !   HYPRE_BUF(local_vol) = 0.0
    !   if ( this%kind == p_fk_bphi_iter ) then
    !     HYPRE_BUF(local_vol-2) = 2.0 * idr2
    !   endif

    case ( p_fk_bperp )

      do i = 1, local_vol, this%num_stencil*4
        k = k + 1
        k1 = k - 0.5
        k2 = k - 1.0
        ! stencil for Re(Br)
        HYPRE_BUF(i)    = k2 / k1 * idr2
        HYPRE_BUF(i+1)  = 0.0
        HYPRE_BUF(i+2)  = 0.0
        HYPRE_BUF(i+3)  = -idr2 * (2.0 + (m2+1.0) / k1**2)
        HYPRE_BUF(i+4)  = 0.0
        HYPRE_BUF(i+5)  = 2.0*m / k1**2 * idr2
        HYPRE_BUF(i+6)  = k / k1 * idr2

        ! stencil for Im(Br)
        HYPRE_BUF(i+7)   = k2 / k1 * idr2
        HYPRE_BUF(i+8)   = 0.0
        HYPRE_BUF(i+9)   = 0.0
        HYPRE_BUF(i+10)  = -idr2 * (2.0 + (m2+1.0) / k1**2)
        HYPRE_BUF(i+11)  = -2.0*m / k1**2 * idr2
        HYPRE_BUF(i+12)  = 0.0
        HYPRE_BUF(i+13)  = k / k1 * idr2

        ! stencil for Re(Bphi)
        HYPRE_BUF(i+14)  = k2 / k1 * idr2
        HYPRE_BUF(i+15)  = 0.0
        HYPRE_BUF(i+16)  = -2.0*m / k1**2 * idr2
        HYPRE_BUF(i+17)  = -idr2 * (2.0 + (m2+1.0) / k1**2)
        HYPRE_BUF(i+18)  = 0.0
        HYPRE_BUF(i+19)  = 0.0
        HYPRE_BUF(i+20)  = k / k1 * idr2

        ! stencil for Im(Bphi)
        HYPRE_BUF(i+21)  = k2 / k1 * idr2
        HYPRE_BUF(i+22)  = 2.0*m / k1**2 * idr2
        HYPRE_BUF(i+23)  = 0.0
        HYPRE_BUF(i+24)  = -idr2 * (2.0 + (m2+1.0) / k1**2)
        HYPRE_BUF(i+25)  = 0.0
        HYPRE_BUF(i+26)  = 0.0
        HYPRE_BUF(i+27)  = k / k1 * idr2
      enddo

      ! boundary condition to be set

    case ( p_fk_bperp_iter )

      do i = 1, local_vol, this%num_stencil*4
        k = k + 1
        k1 = k - 0.5
        k2 = k - 1.0
        ! stencil for Re(Br)
        HYPRE_BUF(i)    = k2 / k1 * idr2
        HYPRE_BUF(i+1)  = 0.0
        HYPRE_BUF(i+2)  = 0.0
        HYPRE_BUF(i+3)  = -1.0 - idr2 * (2.0 + (m2+1.0) / k1**2)
        HYPRE_BUF(i+4)  = 0.0
        HYPRE_BUF(i+5)  = 2.0*m / k1**2 * idr2
        HYPRE_BUF(i+6)  = k / k1 * idr2

        ! stencil for Im(Br)
        HYPRE_BUF(i+7)   = k2 / k1 * idr2
        HYPRE_BUF(i+8)   = 0.0
        HYPRE_BUF(i+9)   = 0.0
        HYPRE_BUF(i+10)  = -1.0 - idr2 * (2.0 + (m2+1.0) / k1**2)
        HYPRE_BUF(i+11)  = -2.0*m / k1**2 * idr2
        HYPRE_BUF(i+12)  = 0.0
        HYPRE_BUF(i+13)  = k / k1 * idr2

        ! stencil for Re(Bphi)
        HYPRE_BUF(i+14)  = k2 / k1 * idr2
        HYPRE_BUF(i+15)  = 0.0
        HYPRE_BUF(i+16)  = -2.0*m / k1**2 * idr2
        HYPRE_BUF(i+17)  = -1.0 - idr2 * (2.0 + (m2+1.0) / k1**2)
        HYPRE_BUF(i+18)  = 0.0
        HYPRE_BUF(i+19)  = 0.0
        HYPRE_BUF(i+20)  = k / k1 * idr2

        ! stencil for Im(Bphi)
        HYPRE_BUF(i+21)  = k2 / k1 * idr2
        HYPRE_BUF(i+22)  = 2.0*m / k1**2 * idr2
        HYPRE_BUF(i+23)  = 0.0
        HYPRE_BUF(i+24)  = -1.0 - idr2 * (2.0 + (m2+1.0) / k1**2)
        HYPRE_BUF(i+25)  = 0.0
        HYPRE_BUF(i+26)  = 0.0
        HYPRE_BUF(i+27)  = k / k1 * idr2
      enddo

      HYPRE_BUF(1) = 0.0
      HYPRE_BUF(8) = 0.0
      HYPRE_BUF(15) = 0.0
      HYPRE_BUF(22) = 0.0

      HYPRE_BUF(local_vol-21) = 0.0
      HYPRE_BUF(local_vol-14) = 0.0
      HYPRE_BUF(local_vol-7) = 0.0
      HYPRE_BUF(local_vol) = 0.0

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

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_hypre_matrix

end module field_solver_class