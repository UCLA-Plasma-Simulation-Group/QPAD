module neutral_subcyc_class

use sysutil_module
use neutral_class
use part2d_subcyc_class

implicit none

private

public :: neutral_subcyc

type, extends(neutral) :: neutral_subcyc
  private
  contains
  procedure :: alloc => alloc_neutral_subcyc
  procedure :: get_exp_fac_max => get_exp_fac_max_neutral_subcyc
end type neutral_subcyc

save

character(len=10) :: cls_name = 'neutral_subcyc'
integer, parameter :: cls_level = 2

contains

subroutine alloc_neutral_subcyc(this)

  implicit none
  class(neutral_subcyc), intent(inout) :: this

  if (.not. associated(this%part)) then
     allocate(part2d_subcyc :: this%part)
  endif

end subroutine alloc_neutral_subcyc

subroutine get_exp_fac_max_neutral_subcyc(this, rst)

  implicit none
  class(neutral_subcyc), intent(inout) :: this
  real, intent(out) :: rst

  real :: gam_pz_min_loc, gam_pz_min
  integer :: ierr, npp
  character(len=32), save :: sname = 'get_exp_fac_max_neutral_subcyc'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  select type (obj => this%part)
    type is (part2d_subcyc)
      rst = obj%get_exp_fac_max()
  end select
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine get_exp_fac_max_neutral_subcyc

end module neutral_subcyc_class