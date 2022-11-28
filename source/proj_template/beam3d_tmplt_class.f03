module beam3d_tmplt_class

use parallel_module
use sysutil_module
use field_em_class
use beam3d_class
use part3d_comm
use part3d_class
use part3d_tmplt_class
use options_class
use input_class

implicit none

private

public :: beam3d_tmplt

type, extends( beam3d ) :: beam3d_tmplt

   private

   contains

   procedure :: alloc => alloc_beam3d_tmplt
   procedure :: push  => push_beam3d_tmplt

   ! overwrite the type-bound procedures if necessary
   ! procedure :: new   => init_beam3d_tmplt
   ! procedure :: del   => end_beam3d_tmplt
   ! procedure :: wr    => writehdf5_beam3d_tmplt
   ! procedure :: wrq   => writeq_beam3d_tmplt
   ! procedure :: wrst  => writerst_beam3d_tmplt
   ! procedure :: rrst  => readrst_beam3d_tmplt

   ! generic :: qdp  => qdeposit_beam3d_tmplt, qdeposit_beam3d_tmplt_finish
   ! procedure, private :: qdeposit_beam3d_tmplt, qdeposit_beam3d_tmplt_finish

end type

save

character(len=32) :: cls_name = 'beam3d_tmplt'
integer, parameter :: cls_level = 2

contains

subroutine alloc_beam3d_tmplt( this, input, opts, beam_id )

   implicit none

   class( beam3d_tmplt ), intent(inout) :: this
   type(input_json), intent(inout) :: input
   type(options), intent(in) :: opts
   integer, intent(in) :: beam_id
   ! local data
   character(len=32), save :: sname = 'alloc_beam3d_tmplt'

   call write_dbg( cls_name, sname, cls_level, 'starts' )
   if ( .not. associated( this%part ) ) then
      allocate( part3d_tmplt :: this%part )
   endif
   call this%beam3d%alloc( input, opts, beam_id )
   call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_beam3d_tmplt

subroutine push_beam3d_tmplt( this, ef, bf, tag, sid )

   implicit none

   class(beam3d_tmplt), intent(inout) :: this
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf
   integer, intent(in) :: tag
   integer, intent(inout) :: sid
   ! local data
   class(part3d), pointer :: part_ptr
   character(len=32), save :: sname = 'push_beam3d_tmplt'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   if (.not. this%evol) then
      call write_dbg(cls_name, sname, cls_level, 'ends')
      return
   endif

   select type ( part_ptr => this%part )
   type is ( part3d_tmplt )
      call part_ptr%push_tmplt( ef, bf )
   end select

   call this%part%update_bound()
   call move_part3d_comm( this%part, tag, sid )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_beam3d_tmplt

end module beam3d_tmplt_class