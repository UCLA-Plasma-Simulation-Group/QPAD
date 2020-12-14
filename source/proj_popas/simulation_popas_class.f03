module simulation_popas_class

use simulation_class
use input_class
use sim_beams_popas_class
use diagnostics_popas_class
use sysutil_module

implicit none

private

public :: simulation_popas

type, extends ( simulation ) :: simulation_popas

  ! add new members if necessary

  contains

  procedure :: alloc => alloc_simulation_popas

  ! overwrite the type-bound procedures if necessary
  ! procedure :: new   => init_simulation_popas
  ! procedure :: del   => end_simulation_popas
  ! procedure :: run   => run_simulation_popas

end type simulation_popas

character(len=18), save :: cls_name = 'simulation'
integer, save :: cls_level = 1

contains

subroutine alloc_simulation_popas( this, input )

  implicit none

  class( simulation_popas ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  ! local data
  character(len=18), save :: sname = 'alloc_simulation_popas'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%beams ) ) then
    allocate( sim_beams_popas :: this%beams )
  endif

  if ( .not. associated( this%diag ) ) then
    allocate( sim_diag_popas :: this%diag )
  endif

  call this%simulation%alloc( input )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_simulation_popas

end module simulation_popas_class