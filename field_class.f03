module field_class

use ufield_class
use param
use system

implicit none

private

public :: field

character(len=20), parameter :: cls_name = "field"
integer, parameter :: cls_level = 1

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
  ! generic :: solve => solve_field
  ! generic :: write_hdf5
  ! generic :: smooth
  procedure :: get_dr, get_dxi, get_num_modes

  procedure, private :: init_field, end_field 
  procedure, private :: get_rf_re_all, get_rf_re_mode, get_rf_im_all, get_rf_im_mode
  ! procedure :: read_input_field

end type field


contains

! =====================================================================
! Class field implementation
! =====================================================================
subroutine init_field( this, num_modes, dim, dr, dxi, nd, nvp, gc_num, entity )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: num_modes, dim
  real, intent(in) :: dr, dxi
  integer, intent(in), dimension(2) :: nd, nvp
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
  call this%rf_re(0)%new( dim, nd, nvp, gc_num, has_2d=.true. )
  do i = 1, this%num_modes
    call this%rf_re(i)%new( dim, nd, nvp, gc_num, has_2d=.true. )
    call this%rf_im(i)%new( dim, nd, nvp, gc_num, has_2d=.true. )
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