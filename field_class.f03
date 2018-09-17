module field_class

use ufield_class
use param

implicit none

private

public :: field

type :: field

  ! private

  class( ufield ), dimension(:), pointer :: rf_re => null()
  class( ufield ), dimension(:), pointer :: rf_im => null()

  real, dimension(:), pointer :: buf => null() ! buffer to store source term and solution for padding

  ! class( emf_solver ) :: solver

  real :: dr, dxi
  integer :: num_modes
  ! integer :: solver_order ! order of finite difference solver

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
subroutine init_field( this, num_modes, dim, dr, dxi, nd, nvp, gc_num )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: num_modes, dim
  real, intent(in) :: dr, dxi
  integer, intent(in), dimension(2) :: nd, nvp
  integer, intent(in), dimension(2,2) :: gc_num

  integer :: i

  this%num_modes = num_modes
  this%dr = dr
  this%dxi = dxi

  allocate( this%rf_re(0:num_modes) )
  allocate( this%rf_im(num_modes) )
  call this%rf_re(0)%new( dim, nd, nvp, gc_num, has_2d=.true. )
  do i = 1, this%num_modes
    call this%rf_re(i)%new( dim, nd, nvp, gc_num, has_2d=.true. )
    call this%rf_im(i)%new( dim, nd, nvp, gc_num, has_2d=.true. )
  enddo

end subroutine init_field

subroutine end_field( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i

  call this%rf_re(0)%del()
  do i = 1, this%num_modes
    call this%rf_re(i)%del()
    call this%rf_im(i)%del()
  enddo
  deallocate( this%rf_re, this%rf_im )

  if ( associated( this%buf ) ) deallocate( this%buf )

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
  class( ufield ), dimension(:), pointer :: get_rf_re_all

  get_rf_re_all => this%rf_re

end function get_rf_re_all

function get_rf_re_mode( this, mode )

  implicit none

  class( field ), intent(in) :: this
  integer, intent(in) :: mode
  class( ufield ), pointer :: get_rf_re_mode

  get_rf_re_mode => this%rf_re(mode)

end function get_rf_re_mode

function get_rf_im_all( this )

  implicit none

  class( field ), intent(in) :: this
  class( ufield ), dimension(:), pointer :: get_rf_im_all

  get_rf_im_all => this%rf_im

end function get_rf_im_all

function get_rf_im_mode( this, mode )

  implicit none

  class( field ), intent(in) :: this
  integer, intent(in) :: mode
  class( ufield ), pointer :: get_rf_im_mode

  get_rf_im_mode => this%rf_im(mode)

end function get_rf_im_mode

end module field_class