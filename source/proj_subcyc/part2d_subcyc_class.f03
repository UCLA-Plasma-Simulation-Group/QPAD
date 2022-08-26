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

  real :: gam_pz_min_loc, gam_pz_min
  integer :: ierr

  if (this%npp > 0) then
     gam_pz_min_loc = minval(this%gamma(1:this%npp) - this%p(3, 1:this%npp))
     call mpi_allreduce(gam_pz_min_loc, gam_pz_min, 1, p_dtype_real, MPI_MIN, comm_loc(), ierr)
  else
     gam_pz_min = 1.0
  endif
  rst = 1.0 / gam_pz_min

end function get_exp_fac_max_part2d_subcyc

end module part2d_subcyc_class