module sim_species_class

use species2d_class
use parallel_pipe_class
use grid_class

use system

implicit none

private

public :: sim_species

type sim_species

  ! private

  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()
! ===============================================================================
! THIS PART IS TO BE FINISHED
! ===============================================================================
  ! type( species ), dimension(:), allocatable :: spe
  ! type( fdist2d_wrap ), dimension(:), allocatable :: pf

  ! contains

  ! generic :: new => init_sim_species
  ! generic :: del => end_sim_species

  ! procedure, private :: init_sim_species, end_sim_species
! ===============================================================================

end type sim_species

character(len=18), save :: cls_name = 'sim_species'
integer, save :: cls_level = 2

end module sim_species_class