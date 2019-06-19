module simulation_class

use parallel_class
use parallel_pipe_class
use grid_class
use sim_fields_class
use sim_beams_class
use sim_species_class
use diagnostics_class
use field_class

use input_class
use sys
use param
use mpi

implicit none

private

public :: simulation

type simulation

  ! private

  type( input_json ), pointer :: input => null()
  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()

  type( sim_fields ) :: fields
  type( sim_species ) :: species
  type( sim_beams ) :: beams
  type( sim_diag ) :: diag
  real :: dr, dxi, dt
  integer :: iter, nstep3d, nstep2d, start3d, nbeams, nspecies, tstep
  integer :: ndump, max_mode

  ! pipeline parameters
  integer, dimension(:), allocatable :: tag_field, id_field
  integer, dimension(:), allocatable :: tag_spe, id_spe
  integer, dimension(:), allocatable :: tag_beam, id_beam
  integer, dimension(:), allocatable :: tag_diag, id_diag

  contains

  generic :: new => init_simulation
  generic :: del => end_simulation
  generic :: run => run_simulation

  procedure, private :: init_simulation, end_simulation
  procedure, private :: run_simulation

end type simulation

character(len=18), save :: cls_name = 'simulation'
integer, save :: cls_level = 1

contains

subroutine init_simulation(this)

  implicit none

  class(simulation), intent(inout) :: this
  ! local data
  character(len=18), save :: sname = 'init_simulation'

  real :: n0, dr, dxi, dt, time
  logical :: read_rst

  allocate( this%input )
  call this%input%new()
  this%pp => this%input%pp
  this%gp => this%input%gp

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%dr  = this%gp%get_dr()
  this%dxi = this%gp%get_dxi()
  this%nstep2d = this%gp%get_ndp(2)

  call this%input%get( 'simulation.n0', n0 )
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
  call this%input%get( 'simulation.max_mode', this%max_mode )

  call this%fields%new( this%input )
  call this%beams%new( this%input )
  call this%species%new( this%input, (this%start3d-1)*dt )

  call this%diag%new( this%pp, this%input, this%fields, this%beams, this%species )

  allocate( this%tag_field(10), this%id_field(10) )
  allocate( this%tag_beam(this%nbeams), this%id_beam(this%nbeams) )

  ! allocate(this%tag_spe(this%nspecies),this%tag_beam(this%nbeams))
  ! allocate(this%id_spe(this%nspecies),this%id_beam(this%nbeams))
  ! allocate(this%id_bq(this%nbeams,3),this%tag_bq(this%nbeams,2))

  ! allocate(this%id(9+size(this%diag)))
  ! this%id(:) = MPI_REQUEST_NULL
  ! this%id_spe(:) = MPI_REQUEST_NULL
  this%id_beam(:) = MPI_REQUEST_NULL
  ! this%id_bq(:,:) = MPI_REQUEST_NULL                 
  this%id_field = MPI_REQUEST_NULL

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_simulation

subroutine end_simulation(this)

  implicit none

  class( simulation ), intent(inout) :: this

  ! local data
  character(len=18), save :: sname = 'end_simulation'
  integer :: ierr

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%fields%del()
  call this%beams%del()
  call this%species%del()
  call this%diag%del()
  call this%gp%del()
  call this%pp%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

  ! call MPI_FINALIZE(ierr)

end subroutine end_simulation

subroutine run_simulation( this )

  implicit none

  class( simulation ), intent(inout) :: this

  integer :: i, j, k, l, ierr, id
  integer, dimension(MPI_STATUS_SIZE) :: istat
  character(len=32), save :: sname = 'run_simulation'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call start_tprof( 'total simulation time' )

  do i = this%start3d, this%nstep3d

    this%tstep = i
    call write_stdout( '3D step = '//num2str(i) )

    call this%fields%q_beam%as(0.0)
    call this%fields%q_spe%as(0.0)
    do k = 1, this%nbeams
      call this%beams%beam(k)%qdp( this%fields%q_beam )
    enddo


    this%tag_field(1) = ntag(); call this%fields%cu%pipe_recv( this%tag_field(1) )
    this%tag_field(2) = ntag(); call this%fields%b%pipe_recv( this%tag_field(2) )
    this%tag_field(3) = ntag(); call this%fields%e%pipe_recv( this%tag_field(3) )
    this%tag_field(4) = ntag(); call this%fields%psi%pipe_recv( this%tag_field(4), nslice=2 )

    call this%fields%cu%copy_slice( 1, p_copy_2to1 )
    this%fields%b     = 0.0
    this%fields%e     = 0.0
    this%fields%b_spe = 0.0
    this%fields%psi   = 0.0

    do j = 1, this%nstep2d

      call this%fields%q_beam%copy_slice( j+1, p_copy_2to1 )
      call this%fields%q_beam%smooth_f1()
      call this%fields%b_beam%solve( this%fields%q_beam )
      this%fields%q_spe = 0.0
      do k = 1, this%nspecies
        call this%species%spe(k)%qdp( this%fields%q_spe )
      enddo

      call this%fields%q_spe%get_q_ax2() ! this must be put before smooth
      call this%fields%q_spe%smooth_f1()

      call this%fields%q_spe%copy_slice( j+1, p_copy_1to2 )
      call this%fields%psi%solve( this%fields%q_spe )
      do k = 1, this%nspecies
        call this%species%spe(k)%extpsi( this%fields%psi )
      enddo
      call this%fields%e%solve( this%fields%psi, j+1 )
      call this%fields%b_spe%solve( this%fields%cu )

      do l = 1, this%iter

        call add_f1( this%fields%b_spe, this%fields%b_beam, this%fields%b )
        call this%fields%e%solve( this%fields%b, this%fields%psi )
        this%fields%cu = 0.0
        this%fields%acu = 0.0
        this%fields%amu = 0.0
        do k = 1, this%nspecies
          call this%species%spe(k)%amjdp( this%fields%e, this%fields%b, &
            this%fields%cu, this%fields%amu, this%fields%acu )
        enddo
        call this%fields%acu%smooth_f1()
        call this%fields%amu%smooth_f1()
        call this%fields%cu%smooth_f1()
        call this%fields%dcu%solve( this%fields%acu, this%fields%amu )
        call this%fields%b_spe%solve( this%fields%dcu, this%fields%cu )
        call this%fields%b_spe%solve( this%fields%cu )
        if ( l == this%iter ) then
          do k = 1, this%nspecies
            call this%species%spe(k)%cbq(j+1)
          enddo
          call this%fields%cu%copy_slice( j+1, p_copy_1to2 )
        endif

      enddo ! iteration

      call add_f1( this%fields%b_spe, this%fields%b_beam, this%fields%b )
      call this%fields%e%solve( this%fields%b, this%fields%psi )
      call dot_f1( this%dxi, this%fields%dcu )
      call add_f1( this%fields%dcu, this%fields%cu, (/1,2/), (/1,2/) )
      do k = 1, this%nspecies
        call this%species%spe(k)%push( this%fields%e, this%fields%b )
      enddo
      call this%fields%e%copy_slice( j+1, p_copy_1to2 )
      call this%fields%b%copy_slice( j+1, p_copy_1to2 )
      call this%fields%psi%copy_slice( j+1, p_copy_1to2 )

    enddo ! 2d loop

    call MPI_WAIT( this%id_field(1), istat, ierr )
    call this%fields%cu%pipe_send( this%tag_field(1), this%id_field(1) )
    call MPI_WAIT( this%id_field(2), istat, ierr )
    call this%fields%b%pipe_send( this%tag_field(2), this%id_field(2) )
    call MPI_WAIT( this%id_field(3), istat, ierr )
    call this%fields%e%pipe_send( this%tag_field(3), this%id_field(3) )
    call MPI_WAIT( this%id_field(4), istat, ierr )
    call this%fields%psi%pipe_send( this%tag_field(4), this%id_field(4), nslice=2 )

    do k = 1, this%nbeams
      this%tag_beam(k) = ntag()
      call MPI_WAIT( this%id_beam(k), istat, ierr )
      call this%beams%beam(k)%push( this%fields%e, this%fields%b, &
        this%tag_beam(k), this%tag_beam(k), this%id_beam(k) )
    enddo

    call this%diag%run( this%tstep, this%dt )

    do k = 1, this%nspecies
      call this%species%spe(k)%renew( i*this%dt )
    enddo

  enddo ! 3d loop

  call stop_tprof( 'total simulation time' )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine run_simulation

function ntag()

   implicit none
   integer, save :: tag = 0
   integer :: ntag
           
   ntag = tag
   tag = tag + 1

end function ntag

end module simulation_class