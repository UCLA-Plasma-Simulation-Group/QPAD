module species2d_subcyc_class

use sysutil_module
use species2d_class
use part2d_subcyc_class

implicit none

private

public :: species2d_subcyc

type, extends(species2d) :: species2d_subcyc

   private

   contains

   procedure :: alloc => alloc_species2d_subcyc
   procedure :: get_exp_fac_max => get_exp_fac_max_species2d_subcyc
   procedure :: clamp_exp_fac => clamp_exp_fac_species2d_subcyc
   
end type species2d_subcyc

save

character(len=10) :: cls_name = 'species2d_subcyc'
integer, parameter :: cls_level = 2

contains

subroutine alloc_species2d_subcyc(this)

   implicit none
   class(species2d_subcyc), intent(inout) :: this

   if (.not. associated(this%part)) then
      allocate(part2d_subcyc :: this%part)
   endif

end subroutine alloc_species2d_subcyc

subroutine get_exp_fac_max_species2d_subcyc(this, rst)

   implicit none
   class(species2d_subcyc), intent(inout) :: this
   real, intent(out) :: rst

   character(len=32), save :: sname = 'get_exp_fac_max_species2d_subcyc'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   select type (obj => this%part)
      type is (part2d_subcyc)
         rst = obj%get_exp_fac_max()
   end select
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine get_exp_fac_max_species2d_subcyc

subroutine clamp_exp_fac_species2d_subcyc(this, exp_fac_clamped)

   implicit none
   class(species2d_subcyc), intent(inout) :: this
   real, intent(in) :: exp_fac_clamped

   character(len=32), save :: sname = 'clamp_exp_fac_species2d_subcyc'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   select type (obj => this%part)
      type is (part2d_subcyc)
         call obj%clamp_exp_fac(exp_fac_clamped)
   end select
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine clamp_exp_fac_species2d_subcyc

end module species2d_subcyc_class