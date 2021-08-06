module sim_beams_tmplt_class

use sim_beams_class
use beam3d_tmplt_class
use options_class
use input_class
use sysutil_module

implicit none

private

public :: sim_beams_tmplt

type, extends( sim_beams ) :: sim_beams_tmplt

  ! private

  contains

  procedure :: alloc => alloc_sim_beams_tmplt

  ! overwrite the type-bound procedures if necessary
  ! procedure :: new   => init_sim_beams_tmplt
  ! procedure :: del   => end_sim_beams_tmplt

end type sim_beams_tmplt

character(len=32), save :: cls_name = 'sim_beams_tmplt'
integer, save :: cls_level = 2

contains

subroutine alloc_sim_beams_tmplt( this, input, opts )

  implicit none

  class( sim_beams_tmplt ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  ! local data
  integer :: i
  character(len=32), save :: sname = 'alloc_sim_beams_tmplt'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.nbeams', this%num_beams )

  if ( .not. associated( this%beam ) ) then
    allocate( beam3d_tmplt :: this%beam( this%num_beams ) )
  endif

  call this%sim_beams%alloc( input, opts )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_sim_beams_tmplt

end module sim_beams_tmplt_class