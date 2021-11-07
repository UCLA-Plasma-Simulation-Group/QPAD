module laser_solver_class

use options_class
use parallel_module
use sysutil_module
use pcr_solver_class
use debug_tool

implicit none

private

public :: laser_solver

character(len=32), parameter :: cls_name = "laser_solver"
integer, parameter :: cls_level = 4

type, extends(pcr_solver) :: laser_solver

  contains

  procedure, private :: set_struct_matrix => set_struct_matrix_laser_solver

end type laser_solver

contains

subroutine set_struct_matrix_laser_solver( this, opts, dr )

  implicit none

  class( laser_solver ), intent(inout) :: this
  type( options ), intent(in) :: opts
  real, intent(in) :: dr

  integer :: i, ierr, local_vol, nr, noff, m
  integer :: comm, lidproc, lnvp
  real :: dr2, m2, j, jmax
  character(len=32), save :: sname = "set_struct_matrix_laser_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm    = comm_loc()
  lidproc = id_proc_loc()
  lnvp    = num_procs_loc()
  noff    = opts%get_noff(1)

  dr2 = dr*dr
  m = this%mode
  m2 = real(m*m)
  nr = this%iupper - this%ilower + 1
  local_vol = nr * this%num_stencil

  if ( .not. associated( this%HYPRE_BUF ) ) then
    allocate( this%HYPRE_BUF( local_vol ) )
  elseif ( size(this%HYPRE_BUF) < local_vol ) then
    deallocate( this%HYPRE_BUF )
    allocate( this%HYPRE_BUF( local_vol ) )
  endif

  call HYPRE_StructMatrixCreate( comm, this%grid, this%stencil, this%A, ierr )
  call HYPRE_StructMatrixInitialize( this%A, ierr )

  ! set the matrix elements of each partition
  j = real(noff)
  do i = 1, local_vol, this%num_stencil
    j = j + 1.0
    this%HYPRE_BUF(i)   = 1.0 - 0.5 / j
    this%HYPRE_BUF(i+1) = -2.0 - m2 / j**2
    this%HYPRE_BUF(i+2) = 1.0 + 0.5 / j
  enddo

  ! set the inner boundary
  if (lidproc == 0) then

    if (m == 0) then
      this%HYPRE_BUF(1) = 0.0
      this%HYPRE_BUF(2) = -4.0
      this%HYPRE_BUF(3) = 4.0
    else
      ! matrix elements 1 to 3 are given arbitrarily to make sure the matrix
      ! is not singular. The vanishing of element 4 indicates the on-axis field
      ! value is zero.
      this%HYPRE_BUF(1) = 0.0
      this%HYPRE_BUF(2) = 1.0
      this%HYPRE_BUF(3) = 0.0
      this%HYPRE_BUF(4) = 0.0
    endif

  endif

  ! set the outter boundary
  if ( lidproc == lnvp-1 ) then

    jmax = real(noff + nr)
    if ( this%mode == 0 ) then
      this%HYPRE_BUF(local_vol) = 0.0
    else
      this%HYPRE_BUF(local_vol-1) = this%HYPRE_BUF(local_vol-1) + (1.0-m/jmax) * this%HYPRE_BUF(local_vol)
      this%HYPRE_BUF(local_vol) = 0.0
    endif

  endif

  ! normalize the matrix elements
  this%HYPRE_BUF = this%HYPRE_BUF / dr2

  call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
    this%stencil_idx, this%HYPRE_BUF, ierr )
  call HYPRE_StructMatrixAssemble( this%A, ierr )

  ! release the buffer to save more memory
  deallocate(this%HYPRE_BUF)
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_matrix_laser_solver

end module laser_solver_class
