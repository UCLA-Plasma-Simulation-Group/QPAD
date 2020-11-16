module simulation_class

use parallel_module
use options_class
use sim_fields_class
use sim_beams_class
use sim_species_class
use diagnostics_class
use field_class
use field_psi_class
use field_vpot_class
use field_b_class
use field_e_class
use field_src_class
use beam3d_class
use species2d_class

use input_class
use sysutil_module
use param
use mpi

use debug_tool

implicit none

private

public :: simulation

integer, parameter :: p_max_tag_num = 32

type simulation

  ! private

  class( sim_fields ),  pointer :: fields  => null()
  class( sim_species ), pointer :: species => null()
  class( sim_beams ),   pointer :: beams   => null()
  class( sim_diag ),    pointer :: diag    => null()

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

  procedure :: alloc => alloc_simulation
  procedure :: new   => init_simulation
  procedure :: del   => end_simulation
  procedure :: run   => run_simulation

end type simulation

character(len=18), save :: cls_name = 'simulation'
integer, save :: cls_level = 1

contains

subroutine alloc_simulation( this, input )

  implicit none

  class( simulation ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  ! local data
  character(len=18), save :: sname = 'alloc_simulation'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%fields ) )  allocate( sim_fields :: this%fields )
  if ( .not. associated( this%species ) ) allocate( sim_species :: this%species )
  if ( .not. associated( this%beams ) )   allocate( sim_beams :: this%beams )
  if ( .not. associated( this%diag ) )    allocate( sim_diag :: this%diag )

  call this%fields%alloc( input )
  call this%species%alloc( input )
  call this%beams%alloc( input )
  call this%diag%alloc( input )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_simulation

subroutine init_simulation(this, input, opts)

  implicit none

  class(simulation), intent(inout) :: this
  type(input_json), intent(inout) :: input
  type(options), intent(in) :: opts
  ! local data
  character(len=18), save :: sname = 'init_simulation'

  real :: n0, dt, time
  logical :: read_rst

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%dr  = opts%get_dr()
  this%dxi = opts%get_dxi()
  this%nstep2d = opts%get_ndp(2)

  call input%get( 'simulation.n0', n0 )
  call input%get( 'simulation.time', time )
  call input%get( 'simulation.dt', dt )
  this%nstep3d = time/dt
  this%dt = dt

  call input%get( 'simulation.read_restart', read_rst )
  if (read_rst) then
    call input%get( 'simulation.restart_timestep', this%start3d )
    this%start3d = this%start3d + 1
  else
    this%start3d = 1
  endif

  call input%get( 'simulation.iter', this%iter )
  call input%get( 'simulation.nbeams', this%nbeams )
  call input%get( 'simulation.nspecies', this%nspecies )
  call input%get( 'simulation.max_mode', this%max_mode )

  call this%fields%new( input, opts )
  call this%beams%new( input, opts )
  call this%species%new( input, opts, (this%start3d-1)*dt )

  call this%diag%new( input, opts, this%fields, this%beams, this%species )

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
  integer :: ierr
  character(len=18), save :: sname = 'end_simulation'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%fields%del()
  call this%beams%del()
  call this%species%del()
  call this%diag%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

  call mpi_finalize(ierr)

end subroutine end_simulation

subroutine run_simulation( this )

  implicit none

  class( simulation ), intent(inout) :: this

  integer :: i, j, k, l, ierr
  integer, dimension(MPI_STATUS_SIZE) :: istat
  character(len=32), save :: sname = 'run_simulation'

  type(field_psi), pointer :: psi
  type(field_vpot), pointer :: vpot
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
  vpot   => this%fields%vpot
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

  ! deposit beams and do diagnostics to see the initial distribution if it is
  ! a fresh run
  if ( this%start3d == 1 ) then

    call q_beam%as(0.0)
    call q_spe%as(0.0)
    ! pipeline data transfer for beams
    do k = 1, this%nbeams
      this%tag_bq(k) = ntag()
      call beam(k)%qdp( q_beam, this%tag_bq(k), this%id_bq(k) )
    enddo
    ! do k = 1, this%nbeams
    !   call beam(k)%qdp( this%tag_bq(k) )
    ! enddo
    call this%diag%run( 0, this%dt )

  endif

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
    ! this%tag_field(1) = ntag(); call cu%pipe_recv( this%tag_field(1) )
    ! this%tag_field(2) = ntag(); call b%pipe_recv( this%tag_field(2) )
    ! this%tag_field(3) = ntag(); call e%pipe_recv( this%tag_field(3) )
    ! this%tag_field(4) = ntag(); call psi%pipe_recv( this%tag_field(4), nslice=2 )
    this%tag_field(1) = ntag()
    call cu%pipe_recv( this%tag_field(1), p_mpi_forward, p_cell_inner )

    ! call cu%copy_slice( 1, p_copy_2to1 )
    b     = 0.0
    e     = 0.0
    b_spe = 0.0
    e_spe = 0.0
    psi   = 0.0

    do j = 1, this%nstep2d

      ! finish the beam qdeposit
      ! if ( j == this%nstep2d + 1 ) then
      !   do k = 1, this%nbeams
      !     call beam(k)%qdp( this%tag_bq(k) )
      !   enddo
      ! endif

      call q_beam%copy_slice( j, p_copy_2to1 )
      call q_beam%smooth_f1()
      call b_beam%solve( q_beam )
      q_spe = 0.0
      do k = 1, this%nspecies
        call spe(k)%qdp( q_spe )
      enddo

      call q_spe%copy_slice( j, p_copy_1to2 )
      call psi%solve( q_spe )
      ! do k = 1, this%nspecies
      !   call spe(k)%extpsi( psi )
      ! enddo
      ! call e%solve( psi, j )
      call b_spe%solve( cu )

      do l = 1, this%iter

        call add_f1( b_spe, b_beam, b )
        call e%solve( cu )
        call e%solve( b, psi )
        cu = 0.0
        acu = 0.0
        amu = 0.0
        do k = 1, this%nspecies
          call spe(k)%amjdp( e, b, cu, amu, acu )
        enddo

        call dcu%solve( acu, amu )
        call b_spe%solve( dcu, cu )
        call b_spe%solve( cu )
        if ( l == this%iter ) then
          do k = 1, this%nspecies
            call spe(k)%cbq(j)
          enddo
          call cu%copy_slice( j, p_copy_1to2 )
        endif

      enddo ! iteration

      call add_f1( b_spe, b_beam, b )
      call e_spe%solve( b_spe, psi )
      call e%solve( cu )
      call e%solve( b, psi )

      ! for vector potential diagnostics
      if ( this%diag%has_vpotz .or. this%diag%has_vpott ) then
        if ( this%diag%has_vpotz ) call vpot%solve_vpotz( cu )
        if ( this%diag%has_vpott ) call vpot%solve_vpott( cu )
        call vpot%copy_slice( j, p_copy_1to2 )
      endif

      call dot_f1( this%dxi, dcu )
      call add_f1( dcu, cu, (/1,2/), (/1,2/) )
      do k = 1, this%nspecies
        call spe(k)%push( e, b )
      enddo
      call e%copy_slice( j, p_copy_1to2 )
      call b%copy_slice( j, p_copy_1to2 )
      call psi%copy_slice( j, p_copy_1to2 )
      call b_spe%copy_slice( j, p_copy_1to2 )
      call e_spe%copy_slice( j, p_copy_1to2 )

      if ( j == 1 ) then
        ! call mpi_wait( this%id_field(1), istat, ierr )
        ! call cu%pipe_send( this%tag_field(1), this%id_field(1), p_mpi_backward )
        ! call mpi_wait( this%id_field(2), istat, ierr )
        ! call b%pipe_send( this%tag_field(2), this%id_field(2), p_mpi_backward )
        ! call mpi_wait( this%id_field(3), istat, ierr )
        ! call e%pipe_send( this%tag_field(3), this%id_field(3), p_mpi_backward )
        ! call mpi_wait( this%id_field(4), istat, ierr )
        ! call psi%pipe_send( this%tag_field(4), this%id_field(4), p_mpi_backward )

        call mpi_wait( this%id_field(2), istat, ierr )
        this%tag_field(2) = ntag()
        call b%pipe_send( this%tag_field(2), this%id_field(2), p_mpi_backward, p_cell_inner )
        call mpi_wait( this%id_field(3), istat, ierr )
        this%tag_field(3) = ntag()
        call e%pipe_send( this%tag_field(3), this%id_field(3), p_mpi_backward, p_cell_inner )
        call mpi_wait( this%id_field(4), istat, ierr )
        this%tag_field(4) = ntag()
        call psi%pipe_send( this%tag_field(4), this%id_field(4), p_mpi_backward, p_cell_inner )
      endif

      if ( j == this%nstep2d ) then
        call mpi_wait( this%id_field(1), istat, ierr )
        call cu%pipe_send( this%tag_field(1), this%id_field(1), p_mpi_forward, p_cell_inner )
      endif

    enddo ! 2d loop

    ! pipeline for species
    do k = 1, this%nspecies
      call spe(k)%psend( this%tag_spe(k), this%id_spe(k) )
    enddo

    ! pipeline for fields
    ! call mpi_wait( this%id_field(1), istat, ierr )
    ! call cu%pipe_send( this%tag_field(1), this%id_field(1) )
    ! call mpi_wait( this%id_field(2), istat, ierr )
    ! call b%pipe_send( this%tag_field(2), this%id_field(2) )
    ! call mpi_wait( this%id_field(3), istat, ierr )
    ! call e%pipe_send( this%tag_field(3), this%id_field(3) )
    ! call mpi_wait( this%id_field(4), istat, ierr )
    ! call psi%pipe_send( this%tag_field(4), this%id_field(4), nslice=2 )

    call b%pipe_recv( this%tag_field(2), p_mpi_backward, p_cell_guard )
    call e%pipe_recv( this%tag_field(3), p_mpi_backward, p_cell_guard )
    call psi%pipe_recv( this%tag_field(4), p_mpi_backward, p_cell_guard )

    ! pipeline for beams
    do k = 1, this%nbeams
      this%tag_beam(k) = ntag()
      call mpi_wait( this%id_beam(k), istat, ierr )
      call beam(k)%push( e, b, this%tag_beam(k), this%tag_beam(k), this%id_beam(k) )
    enddo

    call this%diag%run( this%tstep, this%dt )

    do k = 1, this%nspecies
      call mpi_wait( this%id_spe(k), istat, ierr )
      call spe(k)%renew( i*this%dt )
    enddo

  enddo ! 3d loop

  call stop_tprof( 'total simulation time' )

  call write_tprof()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine run_simulation

end module simulation_class