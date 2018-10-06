module field_src_class

use field_class
use param
use system

implicit none

private

character(len=20), parameter :: cls_name = "field_src"
integer, parameter :: cls_level = 1

public :: field_rho, field_jay, field_djdxi

type, extends( field ) :: field_rho

  contains

  generic :: new => init_field_rho
  procedure, private :: init_field_rho

end type field_rho

type, extends( field ) :: field_jay

  contains

  generic :: new => init_field_jay
  procedure, private :: init_field_jay

end type field_jay

type, extends( field ) :: field_djdxi

  contains

  generic :: new => init_field_djdxi
  procedure, private :: init_field_djdxi

end type field_djdxi

contains

subroutine init_field_rho( this, num_modes, dr, dxi, nd, nvp, part_shape )

  implicit none

  class( field_rho ), intent(inout) :: this
  integer, intent(in) :: num_modes, part_shape
  real, intent(in) :: dr, dxi
  integer, intent(in), dimension(2) :: nd, nvp

  integer, dimension(2,2) :: gc_num
  integer :: dim, i
  integer, dimension(2) :: ndp, noff
  character(len=20), save :: sname = "init_field_rho"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ndp = nd / nvp
  noff = (/0,0/)

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/0, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 1
  ! call initialization routine of the parent class
  call this%field%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_rho

subroutine init_field_jay( this, num_modes, dr, dxi, nd, nvp, part_shape )

  implicit none

  class( field_jay ), intent(inout) :: this
  integer, intent(in) :: num_modes, part_shape
  real, intent(in) :: dr, dxi
  integer, intent(in), dimension(2) :: nd, nvp

  integer, dimension(2,2) :: gc_num
  integer :: dim, i
  integer, dimension(2) :: ndp, noff
  character(len=20), save :: sname = "init_field_jay"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ndp = nd / nvp
  noff = (/0,0/)

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/0, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 3
  ! call initialization routine of the parent class
  call this%field%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_jay

subroutine init_field_djdxi( this, num_modes, dr, dxi, nd, nvp, part_shape )

  implicit none

  class( field_djdxi ), intent(inout) :: this
  integer, intent(in) :: num_modes, part_shape
  real, intent(in) :: dr, dxi
  integer, intent(in), dimension(2) :: nd, nvp

  integer, dimension(2,2) :: gc_num
  integer :: dim, i
  integer, dimension(2) :: ndp, noff
  character(len=20), save :: sname = "init_field_djdxi"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ndp = nd / nvp
  noff = (/0,0/)

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/0, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 2
  ! call initialization routine of the parent class
  call this%field%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_djdxi

end module field_src_class