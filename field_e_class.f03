module field_e_class

use field_class
use field_solver_class
use ufield_class
use param

implicit none

private

type, extends( field ) :: field_e

  private

  class( field_solver ), dimension(:), pointer :: solver => null()

  contains

  generic :: new => init_field_e
  procedure :: del => end_field_e
  ! generic :: read_input => read_input_field_psi
  generic :: solve => solve_field_e

  procedure, private :: init_field_e
  procedure, private :: end_field_e
  procedure, private :: sort_src
  procedure, private :: sort_sol
  procedure, private :: solve_field_e

end type field_psi

end module field_e_class