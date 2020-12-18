module sim_beams_popas_class

use sim_beams_class
use beam3d_popas_class
use options_class
use input_class
use sysutil_module

implicit none

private

public :: sim_beams_popas

type, extends( sim_beams ) :: sim_beams_popas

  contains

  procedure :: alloc => alloc_sim_beams_popas

end type sim_beams_popas

character(len=32), save :: cls_name = 'sim_beams_popas'
integer, save :: cls_level = 2

contains

subroutine alloc_sim_beams_popas( this, input )

  implicit none

  class( sim_beams_popas ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  ! local data
  integer :: i
  character(len=32), save :: sname = 'alloc_sim_beams_popas'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.nbeams', this%num_beams )

  if ( .not. associated( this%beam ) ) then
    allocate( beam3d_popas :: this%beam( this%num_beams ) )
  endif

  call this%sim_beams%alloc( input )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_sim_beams_popas

end module sim_beams_popas_class