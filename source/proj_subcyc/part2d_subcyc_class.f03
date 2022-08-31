module part2d_subcyc_class

use parallel_module
use part2d_class
use mpi

implicit none

private

public :: part2d_subcyc

type, extends(part2d) :: part2d_subcyc

  contains
  procedure :: get_exp_fac_max => get_exp_fac_max_part2d_subcyc

end type part2d_subcyc

save

character(len=20), parameter :: cls_name = "part2d_subcyc"
integer, parameter :: cls_level = 2

contains

function get_exp_fac_max_part2d_subcyc(this) result(rst)

  implicit none
  class(part2d_subcyc), intent(inout) :: this
  real :: rst

  real :: exp_fac_max_loc, exp_fac_max
  integer :: ierr

  if (this%npp > 0) then
     exp_fac_max_loc = maxval(this%gamma(1:this%npp) / (this%gamma(1:this%npp) - this%p(3, 1:this%npp)))
     call mpi_allreduce(exp_fac_max_loc, exp_fac_max, 1, p_dtype_real, MPI_MAX, comm_loc(), ierr)
  else
     exp_fac_max = 1.0
  endif
  rst = exp_fac_max

end function get_exp_fac_max_part2d_subcyc

end module part2d_subcyc_class