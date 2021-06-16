module simulation_tmplt_class

use simulation_class
use input_class
use sim_beams_tmplt_class
use sysutil_module
use options_class

implicit none

private

public :: simulation_tmplt

type, extends ( simulation ) :: simulation_tmplt

  ! add new members if necessary

  contains

  procedure :: alloc => alloc_simulation_tmplt

  ! overwrite the type-bound procedures if necessary
  ! procedure :: new   => init_simulation_tmplt
  ! procedure :: del   => end_simulation_tmplt
  ! procedure :: run   => run_simulation_tmplt

end type simulation_tmplt

character(len=18), save :: cls_name = 'simulation'
integer, save :: cls_level = 1

contains

subroutine alloc_simulation_tmplt( this, input, opts )

  implicit none

  class( simulation_tmplt ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  ! local data
  character(len=18), save :: sname = 'alloc_simulation_tmplt'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%beams ) ) then
    allocate( sim_beams_tmplt :: this%beams )
  endif

  call this%simulation%alloc( input, opts )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_simulation_tmplt

end module simulation_tmplt_class