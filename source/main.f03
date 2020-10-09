program quickpic_quasi3d

use sysutil_module
use parallel_module
use options_class
use input_class
use simulation_class

implicit none

type(options) :: opts
type(input_json) :: input_file

class( simulation ), pointer :: sim => null()

integer :: nstages

! initialize MPI and OpenMP environments
call init_parallel( n_threads=0 )

! set up standard output and error stream
call init_stdout( id_proc() )
call init_errors( eunit=2, idproc=id_proc(), monitor=5 )

! create input file
call input_file%new()

! initialize pipeline
call input_file%get( 'simulation.nodes(2)', nstages )
call init_pipeline( nstages )

! initialize simulation options
call opts%new( input_file )

sim => create_simulation( opts )

call sim%new( input_file, opts )
call sim%run()
call sim%del()

contains

function create_simulation(opts) result(sim)

  use param

  implicit none

  type(options), intent(in) :: opts
  class(simulation), pointer :: sim

  ! allocate simulation object according to algorithm
  select case ( opts%algorithm )
  case ( p_sim_standard )
    allocate( simulation :: sim )
  ! case ( p_sim_template )
  !   allocate( simulation_tmplt :: sim )
  case default
    allocate( simulation :: sim )
  end select

end function create_simulation

end program quickpic_quasi3d
