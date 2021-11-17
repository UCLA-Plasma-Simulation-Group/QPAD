module field_laser_class

use parallel_module
use options_class
use field_class
use field_solver_class
use ufield_class
use param
use sysutil_module
use mpi
use fpcr_penta_class

implicit none

private

character(len=20), parameter :: cls_name = "field_laser"
integer, parameter :: cls_level = 3

public :: field_laser

type, extends( field ) :: field_laser

  type( fpcr_penta ), dimension(:), pointer :: pgc_solver => null()
  real, dimension(:), pointer :: buf_re => null(), buf_im => null()

  contains

  generic   :: new   => init_field_laser
  procedure :: alloc => alloc_field_laser
  procedure :: del   => end_field_laser
  procedure :: solve => solve_field_laser
  procedure, private :: set_rhs
  procedure, private :: get_solution
  procedure, private :: init_solver

end type field_laser

contains

subroutine alloc_field_laser( this, max_mode )

  implicit none

  class( field_laser ), intent(inout) :: this
  integer, intent(in) :: max_mode

  if ( .not. associated( this%pgc_solver ) ) then
    allocate( this%pgc_solver(0:max_mode) )
  endif

end subroutine alloc_field_laser

subroutine init_field_laser( this, opts, max_mode, iter )

  implicit none

  class( field_laser ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: max_mode, iter

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  real :: dr
  character(len=20), save :: sname = "init_field_laser"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = opts%get_ndp(1)
  dr  = opts%get_dr()

  gc_num(:,1) = (/1, 1/)
  gc_num(:,2) = (/iter, iter/)

  dim = 1
  ! call initialization routine of the parent class
  call this%field%new( opts, dim, max_mode, gc_num )

  ! initialize solver
  call this%init_solver()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_laser

subroutine end_field_laser( this )

  implicit none

  class( field_laser ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = 'end_field_laser'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%pgc_solver(i)%destroy()
  enddo
  deallocate( this%pgc_solver )
  if ( associated( this%buf_re ) ) deallocate( this%buf_re )
  if ( associated( this%buf_im ) ) deallocate( this%buf_im )
  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_laser

subroutine init_solver( this )

  implicit none
  class( field_laser ), intent(inout) :: this
end subroutine init_solver

subroutine set_source_ez( this, mode, jay_re, jay_im )

  implicit none

  class( field_laser ), intent(inout) :: this
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

  class( field_laser ), intent(inout) :: this
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

subroutine solve_field_laserz( this, jay )

  implicit none

  class( field_laser ), intent(inout) :: this
  class( field_jay ), intent(inout) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_laserz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%max_mode

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

end subroutine solve_field_laserz

subroutine solve_field_laserz_fast( this, psi, idx )

  implicit none

  class( field_laser ), intent(inout) :: this
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
  character(len=20), save :: sname = 'solve_field_laserz_fast'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve ez' )

  psi_re => psi%get_rf_re()
  psi_im => psi%get_rf_im()
  e_re => this%get_rf_re()
  e_im => this%get_rf_im()

  nrp = e_re(0)%get_ndp(1)
  idxih = 0.5 / this%dxi
  idxi = 1.0 / this%dxi

  do i = 0, this%max_mode

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

end subroutine solve_field_laserz_fast

subroutine solve_field_lasert( this, b, psi )

  implicit none

  class( field_laser ), intent(inout) :: this
  class( field_b ), intent(in) :: b
  class( field_psi ), intent(inout) :: psi

  type( ufield ), dimension(:), pointer :: b_re => null(), b_im => null()
  type( ufield ), dimension(:), pointer :: psi_re => null(), psi_im => null()
  real, dimension(:,:), pointer :: ub_re => null(), ub_im => null()
  real, dimension(:,:), pointer :: upsi_re => null(), upsi_im => null()
  real, dimension(:,:), pointer :: ue_re => null(), ue_im => null()
  integer :: mode, i, nrp, noff, idproc, nvp
  real :: idr, idrh, ir
  character(len=20), save :: sname = 'solve_field_lasert'

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
  do mode = 1, this%max_mode

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

end subroutine solve_field_lasert

subroutine solve_field_lasert_beam( this, b )

  implicit none

  class( field_laser ), intent(inout) :: this
  class( field_b ), intent(in) :: b

  type( ufield ), dimension(:), pointer :: b_re => null(), b_im => null()
  real, dimension(:,:), pointer :: ub_re => null(), ub_im => null()
  real, dimension(:,:), pointer :: ue_re => null(), ue_im => null()
  integer :: mode, i, nrp
  character(len=20), save :: sname = 'solve_field_lasert_beam'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve beam et' )

  b_re => b%get_rf_re()
  b_im => b%get_rf_im()
  nrp = this%rf_re(0)%get_ndp(1)

  do mode = 0, this%max_mode

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

end subroutine solve_field_lasert_beam

end module laser_pgc_class