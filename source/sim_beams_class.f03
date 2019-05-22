module sim_beams_class

use beam3d_class
use parallel_pipe_class
use grid_class

use sys

implicit none

private

public :: sim_beams

type sim_beams

  ! private

  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()
! ===============================================================================
! THIS PART IS TO BE FINISHED
! ===============================================================================
  ! type( beam ), dimension(:), allocatable :: beam
  ! type( fdist3d_wrap ), dimension(:), allocatable :: pf

  ! contains

  ! generic :: new => init_sim_beams
  ! generic :: del => end_sim_beams

  ! procedure, private :: init_sim_beams, end_sim_beams
! ===============================================================================

end type sim_beams

character(len=18), save :: cls_name = 'sim_beams'
integer, save :: cls_level = 2


end module sim_beams_class