module field_solver_class

use options_class 
use parallel_module
use mpi
use param
use sysutil_module
use debug_tool
use ufield_class

implicit none

private

public :: field_solver
public :: HYPRE_BUF

character(len=32), parameter :: cls_name = "field_solver"
integer, parameter :: cls_level = 4

real, dimension(:), pointer, save :: HYPRE_BUF => null()

type :: field_solver ! class for HYPRE solver

  ! HYPRE parameters
  integer, dimension(:), pointer :: offsets => null()
  integer, dimension(:), pointer :: stencil_idx => null()
  integer :: num_stencil
  integer :: stype ! currently there's only cyclic reduction
  integer :: kind
  integer :: mode
  integer :: bnd

  integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver
  integer :: iupper, ilower

  contains

  procedure :: new => init_field_solver
  procedure :: del => end_field_solver 
  procedure :: solve => solve_equation
  procedure, private :: set_struct_solver
  procedure, private :: set_struct_grid
  procedure, private :: set_struct_stencil
  procedure, private :: set_struct_matrix 
  procedure, private :: set_struct_matrixv1 
 
end type field_solver

integer, save :: nofff
real, save :: drr
 
contains

! =====================================================================
! Class field_solver implementation
! =====================================================================

subroutine init_field_solver( this, opts, mode, dr, kind, bnd, stype )

  implicit none

  class( field_solver ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: kind, stype, mode, bnd
  real, intent(in) :: dr

  integer :: ierr, comm
  character(len=32), save :: sname = "init_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%stype = stype
  this%mode  = mode
  this%kind  = kind
  this%bnd   = bnd

  ! setup HYPRE grid
  comm = comm_loc()

  call this%set_struct_grid( opts )
  call this%set_struct_stencil()
  call this%set_struct_matrix( opts, dr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%b, ierr )
  call HYPRE_StructVectorInitialize( this%b, ierr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%x, ierr )
  call HYPRE_StructVectorInitialize( this%x, ierr )

  call this%set_struct_solver()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_solver

subroutine end_field_solver( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: ierr
  character(len=32), save :: sname = "end_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( associated(HYPRE_BUF) ) deallocate( HYPRE_BUF )

  call HYPRE_StructGridDestroy( this%grid, ierr )
  call HYPRE_StructStencilDestroy( this%stencil, ierr )
  call HYPRE_StructVectorDestroy( this%b, ierr )
  call HYPRE_StructVectorDestroy( this%x, ierr )
  call HYPRE_StructMatrixDestroy( this%A, ierr )

  call HYPRE_StructCycRedDestroy( this%solver, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_solver

subroutine set_struct_solver( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: ierr, comm
  character(len=32), save :: sname = "set_struct_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()

  ! call write_stdout( 'mode '//num2str(this%mode)//': Using Cyclic Reduction solver' )
  call HYPRE_StructCycRedCreate( comm, this%solver, ierr )
  call HYPRE_StructCycRedSetup( this%solver, this%A, this%b, this%x, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_solver

subroutine solve_equation( this, src_sol, psi_re, q_re)

  implicit none

  class( field_solver ), intent(inout) :: this
  real, intent(inout), dimension(:), pointer :: src_sol
  type(ufield), intent(inout), pointer, optional :: psi_re
  type(ufield), intent(inout), pointer, optional :: q_re
!   real, intent(inout), optional:: psisum
!   real, intent(inout), optional :: qsum

  integer :: ierr
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  character(len=32), save :: sname = "solve_equation"

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
    write(2,*) 'solve_equation'
    f1_re => psi_re%get_f1()
    f2_re => q_re%get_f1()
!     write(2,*) psisum,'solve_equation psi'
!     write(2,*) qsum,'solve_equation q'
    call start_tprof( 'solve plasma bt' )
    write(2,*) 'satrt set_struct_matrixv1'
    call set_struct_matrixv1( this, f1_re, f2_re)
    call this%set_struct_solver()
    write(2,*) 'end set_struct_matrixv1'
  case ( p_fk_vpotz, p_fk_vpotp, p_fk_vpotm )
    call start_tprof( 'solve plasma A' )
  end select 

  call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
  call HYPRE_StructVectorAssemble( this%b, ierr )

  call HYPRE_StructCycRedSolve( this%solver, this%A, this%b, this%x, ierr )

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
  case ( p_fk_vpotz, p_fk_vpotp, p_fk_vpotm )
    call stop_tprof( 'solve plasma A' )
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_equation

subroutine set_struct_grid( this, opts )

  implicit none

  class( field_solver ), intent(inout) :: this
  type( options ), intent(in) :: opts

  integer :: comm, ierr
  character(len=32), save :: sname = "set_struct_grid"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()
  this%ilower = opts%get_noff(1) + 1
  this%iupper = opts%get_noff(1) + opts%get_ndp(1)

  call HYPRE_StructGridCreate( comm, 1, this%grid, ierr )
  call HYPRE_StructGridSetExtents( this%grid, this%ilower, this%iupper, ierr )
  call HYPRE_StructGridAssemble( this%grid, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_grid

subroutine set_struct_stencil( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: i, ierr
  character(len=32), save :: sname = "set_struct_stencil"

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

subroutine set_struct_matrix( this, opts, dr )

  implicit none

  class( field_solver ), intent(inout) :: this
  type( options ), intent(in) :: opts
  real, intent(in) :: dr

  integer :: i, ierr, local_vol, nr, noff, m
  integer :: comm, lidproc, lnvp
  real :: dr2, m2, j, jmax
  character(len=32), save :: sname = "set_struct_matrix"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm    = comm_loc()
  lidproc = id_proc_loc()
  lnvp    = num_procs_loc()
  noff    = opts%get_noff(1)
  nofff   = noff
  drr     = dr

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

  case ( p_fk_psi, p_fk_bt, p_fk_ez, p_fk_bz, p_fk_vpotz )

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
      HYPRE_BUF(i+1) = -2.0 - ((m+1)/j)**2 - dr2
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
      HYPRE_BUF(2) = -2.0 - ((m+1)/j)**2 - dr2
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  case ( p_fk_bminus )

    ! set from the second grid point of each partition
    j = real(noff)
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
      HYPRE_BUF(i+1) = -2.0 - ((m-1)/j)**2 - dr2
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
      HYPRE_BUF(2) = -2.0 - ((m-1)/j)**2 - dr2
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  case ( p_fk_vpotp )

    ! set from the second grid point of each partition
    j = real(noff)
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
      HYPRE_BUF(i+1) = -2.0 - ((m+1)/j)**2
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
      HYPRE_BUF(2) = -2.0 - ((m+1)/j)**2
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  case ( p_fk_vpotm )

    ! set from the second grid point of each partition
    j = real(noff)
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
      HYPRE_BUF(i+1) = -2.0 - ((m-1)/j)**2
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
      HYPRE_BUF(2) = -2.0 - ((m-1)/j)**2
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  end select

  ! set the upper boundary
  if ( lidproc == lnvp-1 ) then

    select case ( this%bnd )

    ! to be deleted
    case ( p_bnd_zero )

      HYPRE_BUF(local_vol) = 0.0

    case ( p_bnd_open )

      jmax = real(noff + nr)

      select case ( this%kind )

      case ( p_fk_psi, p_fk_vpotz )

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

      case ( p_fk_bplus, p_fk_vpotp )

        HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-(m+1)/jmax) * HYPRE_BUF(local_vol)
        HYPRE_BUF(local_vol) = 0.0

      case ( p_fk_bminus, p_fk_vpotm )

        HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-(m+1)/jmax) * HYPRE_BUF(local_vol)
        HYPRE_BUF(local_vol) = 0.0

      case default
        call write_err( 'Invalid field type!' )
      end select

    case default
      call write_err( 'Invalid boundary condition!' )
    end select

  endif

  HYPRE_BUF = HYPRE_BUF / dr2

  call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
    this%stencil_idx, HYPRE_BUF, ierr )

  call HYPRE_StructMatrixAssemble( this%A, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_matrix

subroutine set_struct_matrixv1( this, psi_re, q_re )
! subroutine set_struct_matrixv1( this, psisum, qsum )

  implicit none

  class( field_solver ), intent(inout) :: this
!   type( options ), intent(in) :: opts
!   real, intent(in) :: dr
  real, intent(in), dimension(:,:), optional, pointer :: psi_re
  real, intent(in), dimension(:,:), optional, pointer :: q_re
!   real, intent(in) :: psisum
!   real, intent(in) :: qsum

  integer :: i, ierr, local_vol, nr, noff, m
  integer :: comm, lidproc, lnvp
  real :: dr2, m2, j, jmax, dr
  character(len=32), save :: sname = "set_struct_matrix"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm    = comm_loc()
  lidproc = id_proc_loc()
  lnvp    = num_procs_loc()
  noff    = nofff
  dr      = drr
!   noff    = this%opts%get_noff(1)

  dr2 = drr*drr
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

  case ( p_fk_psi, p_fk_bt, p_fk_ez, p_fk_bz, p_fk_vpotz )

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
!     write(2,*) j, 'j number'
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
!       write(2,*) i, 'i number'
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
!       HYPRE_BUF(i+1) = -2.0 - ((m+1)/j)**2 
      HYPRE_BUF(i+1) = -2.0 - ((m+1)/j)**2 + dr2 * (q_re(1,j)-1)/(1+psi_re(1,j)) 
!       + dr2 /(1836.5-psi_re(1,j))
!         dr2 * 65.000000000000014/(1-0.00054451*psisum)
!       write(2,*) qsum/(1+psisum), 'coefficient ratio'
!       write(2,*) q_re(1,j), 'q_re in field_solver'
!       write(2,*) psi_re(1,j), 'psi_re'
!       write(2,*) 'now'
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
!       HYPRE_BUF(2) = -2.0 - ((m+1)/j)**2 
      HYPRE_BUF(2) = -2.0 - ((m+1)/j)**2 + dr2 * (q_re(1,j)-1)/(1+psi_re(1,j)) 
!       + dr2 /(1836.5-psi_re(1,j))      
!         dr2 * 65.000000000000014/(1-0.00054451*psisum)
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  case ( p_fk_bminus )

    ! set from the second grid point of each partition
    j = real(noff)
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
!       HYPRE_BUF(i+1) = -2.0 - ((m-1)/j)**2 
!       HYPRE_BUF(i+1) = -2.0 - ((m-1)/j)**2 + dr2 * (qsum)/(1+psisum) 
!       HYPRE_BUF(i+1) = -2.0 - ((m-1)/j)**2 + dr2 * (qsum - 65.000000000000014)/(1+psisum) + &
!        dr2 * 65.000000000000014/(1-0.00054451*psisum)
      HYPRE_BUF(i+1) = -2.0 - ((m-1)/j)**2 + dr2 * (q_re(1,j)-1)/(1+psi_re(1,j)) 
!       + dr2 /(1836.5-psi_re(1,j))
!       write(2,*) qsum/(1+psisum), 'coefficient ratio'
!       write(2,*) qsum, 'q_re in field_solver'
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
!       HYPRE_BUF(2) = -2.0 - ((m-1)/j)**2 
      HYPRE_BUF(2) = -2.0 - ((m-1)/j)**2 + dr2 * (q_re(1,j)-1)/(1+psi_re(1,j)) 
!       + dr2 /(1836.5-psi_re(1,j))
!        dr2 * 65.000000000000014/(1-0.00054451*psisum)
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  case ( p_fk_vpotp )

    ! set from the second grid point of each partition
    j = real(noff)
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
      HYPRE_BUF(i+1) = -2.0 - ((m+1)/j)**2
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
      HYPRE_BUF(2) = -2.0 - ((m+1)/j)**2
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  case ( p_fk_vpotm )

    ! set from the second grid point of each partition
    j = real(noff)
    do i = 4, local_vol, this%num_stencil
      j = j + 1.0
      HYPRE_BUF(i)   = 1.0 - 0.5 / j
      HYPRE_BUF(i+1) = -2.0 - ((m-1)/j)**2
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
      HYPRE_BUF(2) = -2.0 - ((m-1)/j)**2
      HYPRE_BUF(3) = 1.0 + 0.5 / j

    endif

  end select

  ! set the upper boundary
  if ( lidproc == lnvp-1 ) then

    select case ( this%bnd )

    ! to be deleted
    case ( p_bnd_zero )

      HYPRE_BUF(local_vol) = 0.0

    case ( p_bnd_open )

      jmax = real(noff + nr)

      select case ( this%kind )

      case ( p_fk_psi, p_fk_vpotz )

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

      case ( p_fk_bplus, p_fk_vpotp )

        HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-(m+1)/jmax) * HYPRE_BUF(local_vol)
        HYPRE_BUF(local_vol) = 0.0

      case ( p_fk_bminus, p_fk_vpotm )

        HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-(m+1)/jmax) * HYPRE_BUF(local_vol)
        HYPRE_BUF(local_vol) = 0.0

      case default
        call write_err( 'Invalid field type!' )
      end select

    case default
      call write_err( 'Invalid boundary condition!' )
    end select

  endif

  HYPRE_BUF = HYPRE_BUF / dr2

  call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
    this%stencil_idx, HYPRE_BUF, ierr )

  call HYPRE_StructMatrixAssemble( this%A, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_matrixv1 

end module field_solver_class
