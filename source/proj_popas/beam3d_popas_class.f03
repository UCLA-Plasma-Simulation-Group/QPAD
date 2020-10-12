module beam3d_popas_class

use sysutil_module
use beam3d_class
use part3d_popas_class
use part3d_class

implicit none

private

public :: beam3d_popas

type, extends( beam3d ) :: beam3d_popas

   private

   contains

   ! overwrite
   procedure :: alloc => alloc_beam3d_popas

   ! new procedures
   procedure :: get_emittance => get_emittance_beam3d_popas
   procedure :: get_ene_spread => get_ene_spread_beam3d_popas

end type

save

character(len=32) :: cls_name = 'beam3d_popas'
integer, parameter :: cls_level = 2

contains

subroutine alloc_beam3d_popas( this )

   implicit none

   class( beam3d_popas ), intent(inout) :: this
   ! local data
   character(len=32), save :: sname = 'alloc_beam3d_popas'

   call write_dbg( cls_name, sname, cls_level, 'starts' )
   if ( .not. associated( this%part ) ) then
      allocate( part3d_popas :: this%part )
   endif
   call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_beam3d_popas

subroutine get_emittance_beam3d_popas( this, fid, tstep, recv_tag, send_tag, id )

  implicit none

  class( beam3d_popas ), intent(inout) :: this
  integer, intent(in) :: fid, tstep, recv_tag, send_tag
  integer, intent(inout) :: id

  character(len=32), save :: sname = 'get_emittance_beam3d_popas'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select type ( obj => this%part )
  type is ( part3d_popas )
    call obj%get_emittance( fid, tstep, recv_tag, send_tag, id )
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_emittance_beam3d_popas

subroutine get_ene_spread_beam3d_popas( this, fid, tstep, recv_tag, send_tag, id )

  implicit none

  class( beam3d_popas ), intent(inout) :: this
  integer, intent(in) :: fid, tstep, recv_tag, send_tag
  integer, intent(inout) :: id

  character(len=32), save :: sname = 'get_ene_spread_beam3d_popas'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select type ( obj => this%part )
  type is ( part3d_popas )
    call obj%get_ene_spread( fid, tstep, recv_tag, send_tag, id )
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_ene_spread_beam3d_popas

end module beam3d_popas_class