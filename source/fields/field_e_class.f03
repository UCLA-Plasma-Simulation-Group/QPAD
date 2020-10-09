module field_e_class

use parallel_module
use options_class
use field_class
use field_b_class
use field_psi_class
use field_src_class
use field_solver_class
use ufield_class
use param
use sysutil
use mpi

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

  generic :: solve => solve_field_ez_fast, solve_field_ez, &
                      solve_field_et, solve_field_et_beam
  generic :: new => init_field_e

  procedure :: init_field_e
  procedure :: del => end_field_e
  procedure, private :: set_source_ez
  procedure, private :: get_solution_ez
  procedure, private :: solve_field_ez
  procedure, private :: solve_field_ez_fast ! to be deleted in future
  procedure, private :: solve_field_et
  procedure, private :: solve_field_et_beam

end type field_e

contains

subroutine init_field_e( this, opts, num_modes, part_shape, boundary, entity )

  implicit none

  class( field_e ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: num_modes, part_shape, entity, boundary

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  real :: dr
  character(len=20), save :: sname = "init_field_e"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = opts%get_ndp(1)
  dr  = opts%get_dr()

  select case ( part_shape )

  case ( p_ps_linear )

    gc_num(:,1) = (/1, 1/)
    gc_num(:,2) = (/2, 1/)

  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 3
  ! call initialization routine of the parent class
  call this%field%new( opts, dim, num_modes, gc_num, entity )

  ! initialize solver
  select case ( entity )
  case ( p_entity_plasma )
    allocate( this%solver_ez( 0:num_modes ) )
    do i = 0, num_modes
      call this%solver_ez(i)%new( opts, i, dr, kind=p_fk_ez, &
        bnd=boundary, stype=p_hypre_cycred )
    enddo
    allocate( this%buf_re(nrp), this%buf_im(nrp) )
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
  class( ufield ), intent(in) :: jay_re
  class( ufield ), intent(in), optional :: jay_im
  integer, intent(in) :: mode

  integer :: i, nrp, noff, idproc, nvp, ierr, i1, i2
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: idrh, idr, dr, dr2, ir, div
  character(len=20), save :: sname = 'set_source_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve ez' )

  nrp    = jay_re%get_ndp(1)
  idr    = 1.0 / this%dr
  idrh   = 0.5 * idr
  noff   = jay_re%get_noff(1)
  nvp    = num_procs_loc()
  idproc = id_proc_loc()
  dr     = this%dr
  dr2    = dr*dr

  f1_re => jay_re%get_f1()
  if ( present(jay_im) ) then
    f1_im => jay_im%get_f1()
  endif

  this%buf_re = 0.0
  if ( present(jay_im) ) this%buf_im = 0.0

  if ( mode == 0 ) then

    i1 = 1
    i2 = nrp
    if ( idproc == 0 ) i1 = 2
    if ( idproc == nvp-1 ) i2 = nrp-2

    div = 0.0
    do i = i1, i2
      ir = idr / real(i+noff-1)
      this%buf_re(i) = idrh * ( f1_re(1,i+1) - f1_re(1,i-1) ) + ir * f1_re(1,i)
      div = div + this%buf_re(i) * real(i+noff-1)
    enddo

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-2)
      this%buf_re(nrp-1) = idrh * ( f1_re(1,nrp) - f1_re(1,nrp-2) ) + ir * f1_re(1,nrp-1)
      ir = idr / real(nrp+noff-1)
      this%buf_re(nrp) = idr * ( f1_re(1,nrp) - f1_re(1,nrp-1) ) + ir * f1_re(1,nrp)
      div = div - idrh * (f1_re(1,nrp-2) + f1_re(1,nrp-1)) * (real(nrp+noff) - 2.5)
    ! else
    !   ir = idr / real(nrp+noff-1)
    !   this%buf_re(nrp) = idrh * ( f1_re(1,nrp+1) - f1_re(1,nrp-1) ) + ir * f1_re(1,nrp)
    !   div = div + this%buf_re(nrp) * real(nrp+noff-1)
    endif

    call MPI_REDUCE( div, this%buf_re(1), 1, p_dtype_real, MPI_SUM, 0, comm_loc(), ierr )
    if ( idproc == 0 ) this%buf_re(1) = -8.0 * this%buf_re(1)

  elseif ( mode > 0 .and. present( jay_im ) ) then

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      this%buf_re(i) = idrh * ( f1_re(1,i+1) - f1_re(1,i-1) ) + ir * f1_re(1,i) - mode * ir * f1_im(2,i)
      this%buf_im(i) = idrh * ( f1_im(1,i+1) - f1_im(1,i-1) ) + ir * f1_im(1,i) + mode * ir * f1_re(2,i)
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      if ( mod(mode,2) == 0 ) then
        this%buf_re(1) = 2.0 * idr * f1_re(1,2) - mode * idr * f1_im(2,2)
        this%buf_im(1) = 2.0 * idr * f1_im(1,2) + mode * idr * f1_re(2,2)
      else
        this%buf_re(1) = 0.0
        this%buf_im(1) = 0.0

        ! since j_r(m=1) is multiplied by factor 8 on axis, the derivative on index=2 is
        ! calculated using forward difference
        if ( mode == 1 ) then
          ir = idr
          this%buf_re(2) = idr * ( f1_re(1,3) - f1_re(1,2) ) + ir * f1_re(1,2) - mode * ir * f1_im(2,2)
          this%buf_im(2) = idr * ( f1_im(1,3) - f1_im(1,2) ) + ir * f1_im(1,2) + mode * ir * f1_re(2,2)
        endif
      endif
    else
      ir = idr / real(noff)
      this%buf_re(1) = idrh * ( f1_re(1,2) - f1_re(1,0) ) + ir * f1_re(1,1) - mode * ir * f1_im(2,1)
      this%buf_im(1) = idrh * ( f1_im(1,2) - f1_im(1,0) ) + ir * f1_im(1,1) + mode * ir * f1_re(2,1)
    endif

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-1)
      this%buf_re(nrp) = idrh * ( 3.0 * f1_re(1,nrp) - 4.0 * f1_re(1,nrp-1) + f1_re(1,nrp-2) ) + ir * f1_re(1,nrp) &
                         -mode * ir * f1_im(2,nrp)
      this%buf_im(nrp) = idrh * ( 3.0 * f1_im(1,nrp) - 4.0 * f1_im(1,nrp-1) + f1_im(1,nrp-2) ) + ir * f1_im(1,nrp) &
                         +mode * ir * f1_re(2,nrp)
    endif

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call stop_tprof( 'solve ez' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_ez

subroutine get_solution_ez( this, mode )

  implicit none

  class( field_e ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve ez' )

  nrp = this%rf_re(mode)%get_ndp(1)

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nrp
    f1_re(3,i) = this%buf_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nrp
      f1_im(3,i) = this%buf_im(i)
    enddo

    ! Ez(m>0) vanishes on axis
    f1_re(3,1) = 0.0
    f1_im(3,1) = 0.0
  endif

  call stop_tprof( 'solve ez' )
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

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_ez

subroutine solve_field_ez_fast( this, psi, idx )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_psi ), intent(in) :: psi
  integer, intent(in) :: idx

  type( ufield ), dimension(:), pointer :: psi_re => null(), psi_im => null()
  type( ufield ), dimension(:), pointer :: e_re => null(), e_im => null()
  real, dimension(:,:), pointer :: e_f1_re => null(), e_f1_im => null()
  real, dimension(:,:,:), pointer :: e_f2_re => null(), e_f2_im => null()
  real, dimension(:,:), pointer :: psi_f1_re => null(), psi_f1_im => null()
  real, dimension(:,:,:), pointer :: psi_f2_re => null(), psi_f2_im => null()
  integer :: i, j, nrp
  real :: idxih, idxi
  character(len=20), save :: sname = 'solve_field_ez_fast'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve ez' )

  psi_re => psi%get_rf_re()
  psi_im => psi%get_rf_im()
  e_re => this%get_rf_re()
  e_im => this%get_rf_im()

  nrp = e_re(0)%get_ndp(1)
  idxih = 0.5 / this%dxi
  idxi = 1.0 / this%dxi

  do i = 0, this%num_modes

    e_f1_re => e_re(i)%get_f1()
    e_f2_re => e_re(i)%get_f2()
    psi_f1_re => psi_re(i)%get_f1()
    psi_f2_re => psi_re(i)%get_f2()

    do j = 1, nrp
      ! quadratic extrapolation
      e_f1_re(3,j) = idxih * ( 3.0 * psi_f1_re(1,j) - 4.0 * psi_f2_re(1,j,idx-1) + psi_f2_re(1,j,idx-2) )

      ! linear extrapolation
      ! e_f1_re(3,j) = idxi * ( psi_f1_re(1,j) - psi_f2_re(1,j,idx-1) )

      ! smooth
      ! e_f1_re(3,j) = ( e_f1_re(3,j) + e_f2_re(3,j,idx-1) + e_f2_re(3,j,idx-2) ) / 3
    enddo

    if ( i == 0 ) cycle

    e_f1_im => e_im(i)%get_f1()
    e_f2_im => e_im(i)%get_f2()
    psi_f1_im => psi_im(i)%get_f1()
    psi_f2_im => psi_im(i)%get_f2()

    do j = 1, nrp
      ! quadratic extrapolation
      e_f1_im(3,j) = idxih * ( 3.0 * psi_f1_im(1,j) - 4.0 * psi_f2_im(1,j,idx-1) + psi_f2_im(1,j,idx-2) )

      ! linear extrapolation
      ! e_f1_im(3,j) = idxi * ( psi_f1_im(1,j) - psi_f2_im(1,j,idx-1) )

      ! smooth
      ! e_f1_im(3,j) = ( e_f1_im(3,j) + e_f2_im(3,j,idx-1) + e_f2_im(3,j,idx-2) ) / 3
    enddo

  enddo

  call stop_tprof( 'solve ez' )

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_ez_fast

subroutine solve_field_et( this, b, psi )

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
  real :: idr, idrh, ir
  character(len=20), save :: sname = 'solve_field_et'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve plasma et' )

  idr = 1.0 / this%dr
  idrh = idr * 0.5
  nrp = this%rf_re(0)%get_ndp(1)

  b_re   => b%get_rf_re()
  b_im   => b%get_rf_im()
  psi_re => psi%get_rf_re()
  psi_im => psi%get_rf_im()

  noff   = this%rf_re(0)%get_noff(1)
  nvp    = num_procs_loc()
  idproc = id_proc_loc()

  ! m=0 mode
  ub_re   => b_re(0)%get_f1()
  upsi_re => psi_re(0)%get_f1()
  ue_re   => this%rf_re(0)%get_f1()

  do i = 1, nrp
    ue_re(1,i) = ub_re(2,i) - idrh * ( upsi_re(1,i+1) - upsi_re(1,i-1) )
    ue_re(2,i) = -ub_re(1,i)
  enddo
  if ( idproc == 0 ) then
    ue_re(1,1) = 0.0
    ue_re(2,1) = 0.0
  endif
  if ( idproc == nvp-1 ) then
    ue_re(1,nrp) = ub_re(2,nrp) + idrh * ( 4.0 * upsi_re(1,nrp-1) - upsi_re(1,nrp-2) - 3.0 * upsi_re(1,nrp) )
  endif

  ! m>0 mode
  do mode = 1, this%num_modes

    ub_re   => b_re(mode)%get_f1()
    ub_im   => b_im(mode)%get_f1()
    upsi_re => psi_re(mode)%get_f1()
    upsi_im => psi_im(mode)%get_f1()
    ue_re   => this%rf_re(mode)%get_f1()
    ue_im   => this%rf_im(mode)%get_f1()

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      ue_re(1,i) = ub_re(2,i) - idrh * ( upsi_re(1,i+1) - upsi_re(1,i-1) )
      ue_im(1,i) = ub_im(2,i) - idrh * ( upsi_im(1,i+1) - upsi_im(1,i-1) )
      ue_re(2,i) = -ub_re(1,i) + ir * mode * upsi_im(1,i)
      ue_im(2,i) = -ub_im(1,i) - ir * mode * upsi_re(1,i)
    enddo

    if ( idproc == 0 ) then
      if ( mode == 1 ) then
        ue_re(1,1) = ub_re(2,1) - idr * upsi_re(1,2)
        ue_im(1,1) = ub_im(2,1) - idr * upsi_im(1,2)
        ue_re(2,1) = -ub_re(1,1) + idr * upsi_im(1,2)
        ue_im(2,1) = -ub_im(1,1) - idr * upsi_re(1,2)
      else
        ue_re(1,1) = 0.0
        ue_im(1,1) = 0.0
        ue_re(2,1) = 0.0
        ue_im(2,1) = 0.0
      endif
    else
      ir = idr / real(noff)
      ue_re(1,1) = ub_re(2,1) - idrh * ( upsi_re(1,2) - upsi_re(1,0) )
      ue_im(1,1) = ub_im(2,1) - idrh * ( upsi_im(1,2) - upsi_im(1,0) )
      ue_re(2,1) = -ub_re(1,1) + ir * mode * upsi_im(1,1)
      ue_im(2,1) = -ub_im(1,1) - ir * mode * upsi_re(1,1)
    endif

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-1)
      ue_re(1,nrp) = ub_re(2,nrp) + idrh * ( 4.0 * upsi_re(1,nrp-1) - upsi_re(1,nrp-2) - 3.0 * upsi_re(1,nrp) )
      ue_im(1,nrp) = ub_im(2,nrp) + idrh * ( 4.0 * upsi_im(1,nrp-1) - upsi_im(1,nrp-2) - 3.0 * upsi_im(1,nrp) )
    endif

  enddo

  call stop_tprof( 'solve plasma et' )

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_et

subroutine solve_field_et_beam( this, b )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_b ), intent(in) :: b

  type( ufield ), dimension(:), pointer :: b_re => null(), b_im => null()
  real, dimension(:,:), pointer :: ub_re => null(), ub_im => null()
  real, dimension(:,:), pointer :: ue_re => null(), ue_im => null()
  integer :: mode, i, nrp
  character(len=20), save :: sname = 'solve_field_et_beam'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve beam et' )

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

  call stop_tprof( 'solve beam et' )

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_et_beam

end module field_e_class