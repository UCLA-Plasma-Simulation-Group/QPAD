module field_vpot_class

use parallel_module
use options_class
use field_class
use field_b_class
use field_psi_class
use field_src_class
use field_solver_class
use ufield_class
use param
use sysutil_module
use mpi

implicit none

private

character(len=20), parameter :: cls_name = "field_vpot"
integer, parameter :: cls_level = 3

public :: field_vpot

type, extends( field ) :: field_vpot

  ! private

  class( field_solver ), dimension(:), pointer :: solver_vpotz => null()
  class( field_solver ), dimension(:), pointer :: solver_vpotp => null()
  class( field_solver ), dimension(:), pointer :: solver_vpotm => null()
  real, dimension(:), pointer :: buf1_re => null(), buf1_im => null()
  real, dimension(:), pointer :: buf2_re => null(), buf2_im => null()

  contains

  generic :: new => init_field_vpot

  procedure :: init_field_vpot
  procedure :: alloc => alloc_field_vpot
  procedure :: del => end_field_vpot
  procedure, private :: set_source_vpotz
  procedure, private :: get_solution_vpotz
  procedure, private :: set_source_vpott
  procedure, private :: get_solution_vpott
  procedure :: solve_vpotz => solve_field_vpotz
  procedure :: solve_vpott => solve_field_vpott

end type field_vpot

contains

subroutine alloc_field_vpot( this, max_mode )

  implicit none

  class( field_vpot ), intent(inout) :: this
  integer, intent(in) :: max_mode

  if ( .not. associated( this%solver_vpotz ) ) then
    allocate( field_solver :: this%solver_vpotz(0:max_mode) )
  endif

  if ( .not. associated( this%solver_vpotp ) ) then
    allocate( field_solver :: this%solver_vpotp(0:max_mode) )
  endif

  if ( .not. associated( this%solver_vpotm ) ) then
    allocate( field_solver :: this%solver_vpotm(0:max_mode) )
  endif

end subroutine alloc_field_vpot

subroutine init_field_vpot( this, opts, max_mode, part_shape, boundary )

  implicit none

  class( field_vpot ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: max_mode, part_shape, boundary

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  real :: dr
  character(len=20), save :: sname = "init_field_vpot"

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
  call this%field%new( opts, dim, max_mode, gc_num, entity=p_entity_plasma )

  ! initialize solver
  do i = 0, max_mode
    call this%solver_vpotz(i)%new( opts, i, dr, kind=p_fk_vpotz, &
      bnd=boundary, stype=p_hypre_cycred )
    call this%solver_vpotp(i)%new( opts, i, dr, kind=p_fk_vpotp, &
      bnd=boundary, stype=p_hypre_cycred )
    call this%solver_vpotm(i)%new( opts, i, dr, kind=p_fk_vpotm, &
      bnd=boundary, stype=p_hypre_cycred )
  enddo
  allocate( this%buf1_re(nrp), this%buf1_im(nrp) )
  allocate( this%buf2_re(nrp), this%buf2_im(nrp) )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_vpot

subroutine end_field_vpot( this )

  implicit none

  class( field_vpot ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = 'end_field_vpot'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%max_mode
    call this%solver_vpotz(i)%del()
    call this%solver_vpotp(i)%del()
    call this%solver_vpotm(i)%del()
  enddo
  deallocate( this%solver_vpotz )
  deallocate( this%solver_vpotp )
  deallocate( this%solver_vpotm )

  if ( associated( this%buf1_re ) ) deallocate( this%buf1_re )
  if ( associated( this%buf2_re ) ) deallocate( this%buf2_re )
  if ( associated( this%buf1_im ) ) deallocate( this%buf1_im )
  if ( associated( this%buf2_im ) ) deallocate( this%buf2_im )

  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_vpot

subroutine set_source_vpotz( this, mode, jay_re, jay_im )

  implicit none

  class( field_vpot ), intent(inout) :: this
  class( ufield ), intent(in) :: jay_re
  class( ufield ), intent(in), optional :: jay_im
  integer, intent(in) :: mode

  integer :: i, nrp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'set_source_vpotz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve vpotz' )

  nrp = jay_re%get_ndp(1)

  f1_re => jay_re%get_f1()
  if ( present(jay_im) ) then
    f1_im => jay_im%get_f1()
  endif

  this%buf1_re = 0.0
  if ( present(jay_im) ) this%buf1_im = 0.0

  if ( mode == 0 ) then

    do i = 1, nrp
      this%buf1_re(i) = -1.0 * f1_re(3,i)
    enddo

  elseif ( mode > 0 .and. present( jay_im ) ) then

    do i = 1, nrp
      this%buf1_re(i) = -1.0 * f1_re(3,i)
      this%buf1_im(i) = -1.0 * f1_im(3,i)
    enddo

  else
    call write_err( 'Invalid input arguments!' )
  endif

  call stop_tprof( 'solve vpotz' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_vpotz

subroutine set_source_vpott( this, mode, jay_re, jay_im )

  implicit none

  class( field_vpot ), intent(inout) :: this
  class( ufield ), intent(in) :: jay_re
  class( ufield ), intent(in), optional :: jay_im
  integer, intent(in) :: mode

  integer :: i, nrp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'set_source_vpott'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve vpott' )

  nrp    = jay_re%get_ndp(1)

  f1_re => jay_re%get_f1()
  this%buf1_re = 0.0
  this%buf2_re = 0.0

  if ( present(jay_im) ) then
    f1_im => jay_im%get_f1()
    this%buf1_im = 0.0
    this%buf2_im = 0.0
  endif

  if ( mode == 0 ) then

    do i = 1, nrp
      this%buf1_re(i) = -f1_re(1,i)
      this%buf2_re(i) = -f1_re(2,i)
    enddo

  elseif ( mode > 0 .and. present(jay_im) ) then

    do i = 1, nrp
      this%buf1_re(i) = -f1_re(1,i) + f1_im(2,i)
      this%buf1_im(i) = -f1_im(1,i) - f1_re(2,i)
      this%buf2_re(i) = -f1_re(1,i) - f1_im(2,i)
      this%buf2_im(i) = -f1_im(1,i) + f1_re(2,i)
    enddo

  else
    call write_err( 'Invalid input arguments!' )
  endif

  call stop_tprof( 'solve vpott' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_vpott

subroutine get_solution_vpotz( this, mode )

  implicit none

  class( field_vpot ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_vpotz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve vpotz' )

  nrp = this%rf_re(mode)%get_ndp(1)

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nrp
    f1_re(3,i) = this%buf1_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nrp
      f1_im(3,i) = this%buf1_im(i)
    enddo

    ! Az (m>0) vanishes on axis
    f1_re(3,1) = 0.0
    f1_im(3,1) = 0.0
  endif

  call stop_tprof( 'solve vpotz' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_vpotz

subroutine get_solution_vpott( this, mode )

  implicit none

  class( field_vpot ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp, idproc
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_vpott'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve vpott' )

  nrp    = this%rf_re(mode)%get_ndp(1)
  idproc = id_proc_loc()

  if ( mode == 0 ) then

    f1_re => this%rf_re(mode)%get_f1()
    do i = 1, nrp
      f1_re(1,i) = this%buf1_re(i) ! Re(Ar)
      f1_re(2,i) = this%buf2_re(i) ! Re(Aphi)
    enddo

    if ( idproc == 0 ) then
      f1_re(1,1) = 0.0
      f1_re(2,1) = 0.0
    endif

  else

    f1_re => this%rf_re(mode)%get_f1()
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nrp
      f1_re(1,i) = 0.5 * ( this%buf1_re(i) + this%buf2_re(i) ) ! Re(Ar)
      f1_im(1,i) = 0.5 * ( this%buf1_im(i) + this%buf2_im(i) ) ! Im(Ar)
      f1_re(2,i) = 0.5 * ( this%buf1_im(i) - this%buf2_im(i) ) ! Re(Aphi)
      f1_im(2,i) = 0.5 * (-this%buf1_re(i) + this%buf2_re(i) ) ! Im(Aphi)
    enddo

    if ( idproc == 0 ) then
      if ( mode /= 1 ) then
        f1_re(1,1) = 0.0
        f1_im(1,1) = 0.0
        f1_re(2,1) = 0.0
        f1_im(2,1) = 0.0
      endif
    endif

  endif

  call stop_tprof( 'solve vpott' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_vpott

subroutine solve_field_vpotz( this, jay )

  implicit none

  class( field_vpot ), intent(inout) :: this
  class( field_jay ), intent(inout) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_vpotz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%max_mode

    if ( i == 0 ) then
      call this%set_source_vpotz( i, jay_re(i) )
      call this%solver_vpotz(i)%solve( this%buf1_re )
      call this%get_solution_vpotz(i)
      cycle
    endif

    call this%set_source_vpotz( i, jay_re(i), jay_im(i) )
    call this%solver_vpotz(i)%solve( this%buf1_re )
    call this%solver_vpotz(i)%solve( this%buf1_im )
    call this%get_solution_vpotz(i)

  enddo

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_vpotz

subroutine solve_field_vpott( this, jay )

  implicit none

  class( field_vpot ), intent(inout) :: this
  class( field_jay ), intent(inout) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_vpott'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%max_mode

    if ( i == 0 ) then
      call this%set_source_vpott( i, jay_re(i) )
      call this%solver_vpotp(i)%solve( this%buf1_re )
      call this%solver_vpotm(i)%solve( this%buf2_re )
      call this%get_solution_vpott(i)
      cycle
    endif

    call this%set_source_vpott( i, jay_re(i), jay_im(i) )
    call this%solver_vpotp(i)%solve( this%buf1_re )
    call this%solver_vpotp(i)%solve( this%buf1_im )
    call this%solver_vpotm(i)%solve( this%buf2_re )
    call this%solver_vpotm(i)%solve( this%buf2_im )
    call this%get_solution_vpott(i)

  enddo

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_vpott

end module field_vpot_class