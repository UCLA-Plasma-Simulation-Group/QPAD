module sim_plasma_subcyc_class

use species2d_subcyc_class
use neutral_subcyc_class
use input_class
use sim_plasma_class
use sysutil_module

implicit none

private

public :: sim_plasma_subcyc

type, extends(sim_plasma) :: sim_plasma_subcyc
  contains
  procedure :: alloc => alloc_sim_plasma_subcyc
  procedure :: get_exp_fac_max => get_exp_fac_max_sim_plasma_subcyc
  procedure :: clamp_exp_fac => clamp_exp_fac_sim_plasma_subcyc
end type sim_plasma_subcyc

character(len=18), save :: cls_name = 'sim_plasma_subcyc'
integer, save :: cls_level = 2

contains

subroutine alloc_sim_plasma_subcyc(this, input)

  implicit none
  class(sim_plasma_subcyc), intent(inout) :: this
  type(input_json), intent(inout) :: input

  integer :: i

  call input%get('simulation.nspecies', this%num_species)
  call input%get('simulation.nneutrals', this%num_neutrals)

  if (.not. associated(this%spe)) allocate(species2d_subcyc :: this%spe(this%num_species))
  if (.not. associated(this%neut)) allocate(neutral_subcyc :: this%neut(this%num_neutrals))

  do i = 1, this%num_species
    call this%spe(i)%alloc()
  enddo

  do i = 1, this%num_neutrals
    call this%neut(i)%alloc()
  enddo

end subroutine alloc_sim_plasma_subcyc

subroutine get_exp_fac_max_sim_plasma_subcyc(this, exp_fac_max)

  implicit none
  class(sim_plasma_subcyc), intent(inout) :: this
  real, intent(inout) :: exp_fac_max

  real :: fac_tmp
  integer :: ierr, i
  character(len=32), save :: sname = 'get_exp_fac_max_sim_plasma_subcyc'

  call write_dbg(cls_name, sname, cls_level, 'starts')

  exp_fac_max = 1.0
  do i = 1, this%num_species
    select type (obj => this%spe(i))
      type is (species2d_subcyc)
        call obj%get_exp_fac_max(fac_tmp)
    end select
    exp_fac_max = max(exp_fac_max, fac_tmp)
  enddo

  do i = 1, this%num_neutrals
    select type (obj => this%neut(i))
      type is (neutral_subcyc)
        call obj%get_exp_fac_max(fac_tmp)
    end select
    exp_fac_max = max(exp_fac_max, fac_tmp)
  enddo

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine get_exp_fac_max_sim_plasma_subcyc

subroutine clamp_exp_fac_sim_plasma_subcyc(this, exp_fac_clamped)

  implicit none
  class(sim_plasma_subcyc), intent(inout) :: this
  real, intent(in) :: exp_fac_clamped

  integer :: i
  character(len=32), save :: sname = 'clamp_exp_fac_sim_plasma_subcyc'

  call write_dbg(cls_name, sname, cls_level, 'starts')

  do i = 1, this%num_species
    select type (obj => this%spe(i))
      type is (species2d_subcyc)
        call obj%clamp_exp_fac(exp_fac_clamped)
    end select
  enddo

  do i = 1, this%num_neutrals
    select type (obj => this%neut(i))
      type is (neutral_subcyc)
        call obj%clamp_exp_fac(exp_fac_clamped)
    end select
  enddo

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine clamp_exp_fac_sim_plasma_subcyc

end module sim_plasma_subcyc_class