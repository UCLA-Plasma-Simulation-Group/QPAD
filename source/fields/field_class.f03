module field_class

use parallel_pipe_class
use grid_class
use ufield_class
use hdf5io_class
use param
use system

implicit none

private

public :: field

character(len=20), parameter :: cls_name = "field"
integer, parameter :: cls_level = 3

type :: field

  ! private

  class( ufield ), dimension(:), pointer :: rf_re => null()
  class( ufield ), dimension(:), pointer :: rf_im => null()

  ! class( emf_solver ) :: solver

  real :: dr, dxi
  integer :: num_modes
  integer :: entity

  contains

  generic :: new => init_field
  ! generic :: read_input => read_input_field
  procedure :: del => end_field
  generic :: get_rf_re => get_rf_re_all, get_rf_re_mode
  generic :: get_rf_im => get_rf_im_all, get_rf_im_mode
  generic :: copy_gc => copy_gc_local, copy_gc_stage
  generic :: write_hdf5 => write_hdf5_single, write_hdf5_pipe
  ! generic :: solve => solve_field
  ! generic :: smooth
  procedure :: copy_slice
  procedure :: get_dr, get_dxi, get_num_modes

  procedure, private :: init_field, end_field 
  procedure, private :: get_rf_re_all, get_rf_re_mode, get_rf_im_all, get_rf_im_mode
  procedure, private :: copy_gc_local, copy_gc_stage
  procedure, private :: write_hdf5_single, write_hdf5_pipe

end type field


contains

! =====================================================================
! Class field implementation
! =====================================================================
subroutine init_field( this, pp, gp, dim, dr, dxi, num_modes, gc_num, entity )

  implicit none

  class( field ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, dim
  real, intent(in) :: dr, dxi
  integer, intent(in), dimension(2,2) :: gc_num
  integer, intent(in), optional :: entity

  integer :: i
  character(len=20), save :: sname = "init_field"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%num_modes = num_modes
  this%dr = dr
  this%dxi = dxi

  if ( present(entity) ) then
    this%entity = entity
  else
    this%entity = p_entity_none
  endif

  allocate( this%rf_re(0:num_modes) )
  allocate( this%rf_im(num_modes) )
  call this%rf_re(0)%new( pp, gp, dim, gc_num, has_2d=.true. )
  do i = 1, this%num_modes
    call this%rf_re(i)%new( pp, gp, dim, gc_num, has_2d=.true. )
    call this%rf_im(i)%new( pp, gp, dim, gc_num, has_2d=.true. )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field

subroutine end_field( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "end_field"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%rf_re(0)%del()
  do i = 1, this%num_modes
    call this%rf_re(i)%del()
    call this%rf_im(i)%del()
  enddo
  deallocate( this%rf_re, this%rf_im )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field

subroutine copy_slice( this, idx, dir )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: idx, dir

  integer :: i
  character(len=20), save :: sname = "copy_slice"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_slice( idx, dir )
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_slice( idx, dir )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_slice

subroutine copy_gc_local( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "copy_gc_local"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_gc()
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_gc()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_local

subroutine copy_gc_stage( this, dir )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: dir

  integer :: i
  character(len=20), save :: sname = "copy_gc_stage"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_gc( dir )
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_gc( dir )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_stage

subroutine write_hdf5_single( this, files, dim )

  implicit none

  class( field ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim

  integer :: i
  character(len=32), save :: sname = 'write_hdf5_single'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%rf_re(i)%write_hdf5( files(1), dim )
      cycle
    endif

    call this%rf_re(i)%write_hdf5( files(2*i), dim )
    call this%rf_im(i)%write_hdf5( files(2*i+1), dim )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_single

subroutine write_hdf5_pipe( this, files, dim, rtag, stag, id )

  implicit none

  class( field ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim, rtag, stag
  integer, intent(inout) :: id

  integer :: i
  character(len=32), save :: sname = 'write_hdf5_pipe'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%rf_re(i)%write_hdf5( files(1), dim, rtag, stag, id )
      cycle
    endif

    call this%rf_re(i)%write_hdf5( files(2*i), dim, rtag, stag, id )
    call this%rf_im(i)%write_hdf5( files(2*i+1), dim, rtag, stag, id )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_pipe

function get_dr( this )

  implicit none

  class( field ), intent(in) :: this
  real :: get_dr

  get_dr = this%dr

end function get_dr

function get_dxi( this )

  implicit none

  class( field ), intent(in) :: this
  real :: get_dxi

  get_dxi = this%dxi
  
end function get_dxi

function get_num_modes( this )

  implicit none

  class( field ), intent(in) :: this
  integer :: get_num_modes

  get_num_modes = this%num_modes
  
end function get_num_modes

function get_rf_re_all( this )

  implicit none

  class( field ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_rf_re_all

  get_rf_re_all => this%rf_re

end function get_rf_re_all

function get_rf_re_mode( this, mode )

  implicit none

  class( field ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_rf_re_mode

  get_rf_re_mode => this%rf_re(mode)

end function get_rf_re_mode

function get_rf_im_all( this )

  implicit none

  class( field ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_rf_im_all

  get_rf_im_all => this%rf_im

end function get_rf_im_all

function get_rf_im_mode( this, mode )

  implicit none

  class( field ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_rf_im_mode

  get_rf_im_mode => this%rf_im(mode)

end function get_rf_im_mode

end module field_class