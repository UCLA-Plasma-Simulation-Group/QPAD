module simulation_class

use parallel_class
use parallel_pipe_class
use grid_class
use sim_fields_class
use sim_beams_class
use sim_species_class
use diagnostics_class
use field_class
use field_psi_class
use field_b_class
use field_e_class
use field_src_class
use beam3d_class
use species2d_class

use input_class
use sysutil
use param
use mpi

implicit none

private

public :: simulation

integer, parameter :: p_max_tag_num = 32

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
  integer, dimension(:), allocatable :: tag_bq, id_bq
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

  real :: n0, dt, time
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

  allocate( this%tag_field(p_max_tag_num), this%id_field(p_max_tag_num) )
  allocate( this%tag_beam(this%nbeams), this%id_beam(this%nbeams) )
  allocate( this%tag_spe(this%nspecies), this%id_spe(this%nspecies) )
  allocate( this%tag_bq(this%nbeams), this%id_bq(this%nbeams) )

  this%id_field(:) = MPI_REQUEST_NULL
  this%id_spe(:)   = MPI_REQUEST_NULL
  this%id_beam(:)  = MPI_REQUEST_NULL
  this%id_bq(:)    = MPI_REQUEST_NULL

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_simulation

subroutine end_simulation(this)

  implicit none

  class( simulation ), intent(inout) :: this

  ! local data
  character(len=18), save :: sname = 'end_simulation'

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

  integer :: i, j, k, l, ierr
  integer, dimension(MPI_STATUS_SIZE) :: istat
  character(len=32), save :: sname = 'run_simulation'

  type(field_psi), pointer :: psi
  type(field_e), pointer :: e_spe, e_beam, e
  type(field_b), pointer :: b_spe, b_beam, b
  type(field_jay), pointer :: cu, amu
  type(field_rho), pointer :: q_spe, q_beam
  type(field_djdxi), pointer :: dcu, acu
  type(beam3d), dimension(:), pointer :: beam
  type(species2d), dimension(:), pointer :: spe

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call start_tprof( 'total simulation time' )

  psi    => this%fields%psi
  e_spe  => this%fields%e_spe
  e_beam => this%fields%e_beam
  e      => this%fields%e
  b_spe  => this%fields%b_spe
  b_beam => this%fields%b_beam
  b      => this%fields%b
  cu     => this%fields%cu
  amu    => this%fields%amu
  q_spe  => this%fields%q_spe
  q_beam => this%fields%q_beam
  dcu    => this%fields%dcu
  acu    => this%fields%acu

  beam => this%beams%beam
  spe  => this%species%spe

  call this%diag%run( 0, this%dt )

  do i = this%start3d, this%nstep3d

    this%tstep = i
    call write_stdout( '3D step = '//num2str(i) )

    call q_beam%as(0.0)
    call q_spe%as(0.0)

    ! pipeline data transfer for beams
    do k = 1, this%nbeams
      this%tag_bq(k) = ntag()
      call beam(k)%qdp( q_beam, this%tag_bq(k), this%id_bq(k) )
    enddo

    ! pipeline data transfer for species
    do k = 1, this%nspecies
      this%tag_spe(k) = ntag()
      call spe(k)%precv( this%tag_spe(k) )
    enddo

    ! pipeline data transfer for fields
    this%tag_field(1) = ntag(); call cu%pipe_recv( this%tag_field(1) )
    this%tag_field(2) = ntag(); call b%pipe_recv( this%tag_field(2) )
    this%tag_field(3) = ntag(); call e%pipe_recv( this%tag_field(3) )
    this%tag_field(4) = ntag(); call psi%pipe_recv( this%tag_field(4), nslice=2 )

    call cu%copy_slice( 1, p_copy_2to1 )
    b     = 0.0
    e     = 0.0
    b_spe = 0.0
    psi   = 0.0

    do j = 1, this%nstep2d

      ! finish the beam qdeposit
      if ( j == this%nstep2d ) then
        do k = 1, this%nbeams
          call beam(k)%qdp( this%tag_bq(k) )
        enddo
      endif

      call q_beam%copy_slice( j+1, p_copy_2to1 )
      call q_beam%smooth_f1()
      call b_beam%solve( q_beam )
      q_spe = 0.0
      do k = 1, this%nspecies
        call spe(k)%qdp( q_spe )
      enddo

      ! call q_spe%get_q_ax2() ! this must be put before smooth
      call q_spe%smooth_f1()

      call q_spe%copy_slice( j+1, p_copy_1to2 )
      call psi%solve( q_spe )
      do k = 1, this%nspecies
        call spe(k)%extpsi( psi )
      enddo
      call e%solve( psi, j+1 )
      call b_spe%solve( cu )

      do l = 1, this%iter

        call add_f1( b_spe, b_beam, b )
        call e%solve( b, psi )
        cu = 0.0
        acu = 0.0
        amu = 0.0
        do k = 1, this%nspecies
          call spe(k)%amjdp( e, b, cu, amu, acu )
        enddo
        call acu%smooth_f1()
        call amu%smooth_f1()
        call cu%smooth_f1()
        call dcu%solve( acu, amu )
        call b_spe%solve( dcu, cu )
        call b_spe%solve( cu )
        if ( l == this%iter ) then
          do k = 1, this%nspecies
            call spe(k)%cbq(j+1)
          enddo
          call cu%copy_slice( j+1, p_copy_1to2 )
        endif

      enddo ! iteration

      call add_f1( b_spe, b_beam, b )
      call e%solve( b, psi )
      call dot_f1( this%dxi, dcu )
      call add_f1( dcu, cu, (/1,2/), (/1,2/) )
      do k = 1, this%nspecies
        call spe(k)%push( e, b )
      enddo
      call e%copy_slice( j+1, p_copy_1to2 )
      call b%copy_slice( j+1, p_copy_1to2 )
      call psi%copy_slice( j+1, p_copy_1to2 )

    enddo ! 2d loop

    ! pipeline for species
    do k = 1, this%nspecies
      call spe(k)%psend( this%tag_spe(k), this%id_spe(k) )
    enddo

    ! pipeline for fields
    call MPI_WAIT( this%id_field(1), istat, ierr )
    call cu%pipe_send( this%tag_field(1), this%id_field(1) )
    call MPI_WAIT( this%id_field(2), istat, ierr )
    call b%pipe_send( this%tag_field(2), this%id_field(2) )
    call MPI_WAIT( this%id_field(3), istat, ierr )
    call e%pipe_send( this%tag_field(3), this%id_field(3) )
    call MPI_WAIT( this%id_field(4), istat, ierr )
    call psi%pipe_send( this%tag_field(4), this%id_field(4), nslice=2 )

    ! pipeline for beams
    do k = 1, this%nbeams
      this%tag_beam(k) = ntag()
      call MPI_WAIT( this%id_beam(k), istat, ierr )
      call beam(k)%push( e, b, this%tag_beam(k), this%tag_beam(k), this%id_beam(k) )
    enddo

    call this%diag%run( this%tstep, this%dt )

    do k = 1, this%nspecies
      call MPI_WAIT( this%id_spe(k), istat, ierr )
      call spe(k)%renew( i*this%dt )
    enddo

  enddo ! 3d loop

  call stop_tprof( 'total simulation time' )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine run_simulation

end module simulation_class