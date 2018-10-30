module simulation_class

use parallel_class
use parallel_pipe_class
use grid_class
use sim_fields_class
use sim_beams_class
use sim_species_class
use diagnostics_class

use input_class
use system
use param
use mpi

implicit none

private

public :: simulation

type simulation

  private

  type( input_json ), pointer :: input => null()
  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()

  type( sim_fields ) :: fields
  type( sim_species ) :: species
  type( sim_beams ) :: beams
  type( sim_diag ) :: diag
  real :: dr, dxi, dt
  integer :: iter, nstep3d, nstep2d, start3d, nbeams, nspecies, tstep
  integer :: num_modes, interp

  contains

  generic :: new => init_simulation
  generic :: del => end_simulation
  ! generic :: run => run_simulation

  procedure, private :: init_simulation, end_simulation
  ! procedure, private :: run_simulation

end type simulation

contains

subroutine init_simulation(this)

  implicit none

  class(simulation), intent(inout) :: this
  ! local data
  character(len=18), save :: sname = 'init_simulation:'
  character(len=18), save :: cls_name = 'simulation'
  integer, save :: cls_level = 1

  real :: min, max, n0, dr, dxi, dt, time
  integer :: nr, nz
  logical :: read_rst
  character(len=:), allocatable :: interp_str

  allocate( this%input )
  call this%input%new()
  this%pp => this%input%pp
  this%gp => this%input%gp

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%input%get( 'simulation.n0', n0 )

  call this%input%get( 'simulation.box.r(1)', min )
  call this%input%get( 'simulation.box.r(2)', max )
  nr = this%gp%get_nd(1)
  dr = ( max - min ) / real(nr)
  call this%input%get( 'simulation.box.z(1)', min )
  call this%input%get( 'simulation.box.z(2)', max )
  nz = this%gp%get_nd(2)
  dxi = ( max - min ) / real(nz)

  this%dr = dr
  this%dxi = dxi

  this%nstep2d = this%gp%get_ndp(2)

  call this%input%get( 'simulation.time', time )
  call this%input%get( 'simulation.dt', dt )
  this%nstep3d = time/dt
  this%dt = dt

  call this%input%get( 'simulation.read_restart', read_rst )
  if (read_rst) then
    call this%input%get( 'simulation.restart_timestep', this%start3d )
    this%start3d = this%start3d + 1
  else
    this%start3d = 1
  endif

  call this%input%get( 'simulation.iter', this%iter )
  call this%input%get( 'simulation.nbeams', this%nbeams )
  call this%input%get( 'simulation.nspecies', this%nspecies )
  call this%input%get( 'simulation.num_modes', this%num_modes )

  call this%input%get( 'simulation.interp', interp_str )
  select case ( trim(interp_str) )
  case ( 'linear' )
    this%interp = p_ps_linear
  case ( 'quadratic' )
    this%interp = p_ps_quadratic
  case default
    call write_err( 'Invalid interpolation type!' )
  end select

  call this%fields%new( this%pp, this%gp, this%dr, this%dxi, this%num_modes, this%interp )
! ===============================================================================
! THIS PART IS TO BE FINISHED
! ===============================================================================
  ! call this%beams%new( ... )
  ! call this%species%new( ... )
! ===============================================================================

  call this%diag%new( this%pp, this%input, this%fields, this%beams, this%species )
  ! call this%beams%new(this%in,this%fields)
  ! call this%species%new(this%in,this%fields,(this%start3d-1)*dt)

  ! call this%init_diag()                 

  ! allocate(this%tag_spe(this%nspecies),this%tag_beam(this%nbeams))
  ! allocate(this%id_spe(this%nspecies),this%id_beam(this%nbeams))
  ! allocate(this%id_bq(this%nbeams,3),this%tag_bq(this%nbeams,2))

  ! allocate(this%id(9+size(this%diag)))
  ! this%id(:) = MPI_REQUEST_NULL
  ! this%id_spe(:) = MPI_REQUEST_NULL
  ! this%id_beam(:) = MPI_REQUEST_NULL                 
  ! this%id_bq(:,:) = MPI_REQUEST_NULL                 

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_simulation

subroutine end_simulation(this)

  implicit none

  class( simulation ), intent(inout) :: this

  ! local data
  character(len=18), save :: sname = 'end_simulation'
  character(len=18), save :: cls_name = 'simulation'
  integer, save :: cls_level = 1
  integer :: ierr

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%fields%del()
! ===============================================================================
! THIS PART IS TO BE FINISHED
! ===============================================================================
  ! call this%beams%del()
  ! call this%species%del()
! ===============================================================================

  call this%diag%del()
  call this%gp%del()
  call this%pp%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_simulation

end module simulation_class