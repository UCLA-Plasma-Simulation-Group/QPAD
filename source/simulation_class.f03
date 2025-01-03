module simulation_class

use parallel_module
use options_class
use sim_fields_class
use sim_beams_class
use sim_plasma_class
use sim_lasers_class
use diagnostics_class
use field_class
use field_complex_class
use field_psi_class
use field_vpot_class
use field_b_class
use field_e_class
use field_src_class
use field_laser_class
use beam3d_class
use species2d_class
use neutral_class

use input_class
use sysutil_module
use param
use mpi

use debug_tool

implicit none

private

public :: simulation
public :: convergence_tester

integer, parameter :: p_max_tag_num = 32

type simulation

  ! private

  class( sim_fields ), pointer :: fields => null()
  class( sim_plasma ), pointer :: plasma => null()
  class( sim_beams ),  pointer :: beams  => null()
  class( sim_lasers ), pointer :: lasers => null()
  class( sim_diag ),   pointer :: diag   => null()

  real :: dr, dxi, dt
  real :: iter_reltol, iter_abstol
  integer :: iter_max, nstep3d, nstep2d, start2d, start3d, tstep
  integer :: nbeams, nspecies, nneutrals, nlasers
  integer :: ndump, max_mode

  ! pipeline parameters
  integer, dimension(:), allocatable :: tag_field, id_field
  integer, dimension(:), allocatable :: tag_spe, id_spe
  integer, dimension(:,:), allocatable :: tag_neut, id_neut
  integer, dimension(:), allocatable :: tag_beam, id_beam
  integer, dimension(:,:), allocatable :: tag_laser, id_laser
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

subroutine alloc_simulation( this, input, opts )

  implicit none

  class( simulation ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  ! local data
  character(len=18), save :: sname = 'alloc_simulation'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%fields ) ) allocate( sim_fields :: this%fields )
  if ( .not. associated( this%plasma ) ) allocate( sim_plasma :: this%plasma )
  if ( .not. associated( this%beams ) )  allocate( sim_beams :: this%beams )
  if ( .not. associated( this%lasers ) ) allocate( sim_lasers :: this%lasers )
  if ( .not. associated( this%diag ) )   allocate( sim_diag :: this%diag )

  call this%fields%alloc( input )
  call this%plasma%alloc( input )
  call this%beams%alloc( input, opts )
  call this%lasers%alloc( input, opts )
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
  integer :: rnd_seed, num_seeds
  integer, dimension(:), allocatable :: seed

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call write_stdout( 'Initializing simulation...' )

  ! initialize pseudo-random number sequence
  call input%get( 'simulation.random_seed', rnd_seed )
  if ( rnd_seed == 0 ) then
    ! OS generated seed
    call write_stdout( 'Using OS-generated seeds for pseudo-random numbers.' )
    call random_seed()
  else
    ! user specified seeds
    call write_stdout( 'Using user-specified seeds for pseudo-random numbers.' )
    call random_seed( size=num_seeds )
    allocate( seed(num_seeds) )
    seed = rnd_seed
    call random_seed( put=seed )
  endif

  this%dr  = opts%get_dr()
  this%dxi = opts%get_dxi()
  this%nstep2d = opts%get_ndp(2)
  this%start2d = opts%get_noff(2) + 1

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

  call input%get( 'simulation.iter_max', this%iter_max )
  call input%get( 'simulation.iter_reltol', this%iter_reltol )
  call input%get( 'simulation.iter_abstol', this%iter_abstol )
  call input%get( 'simulation.nbeams', this%nbeams )
  call input%get( 'simulation.nspecies', this%nspecies )
  call input%get( 'simulation.nneutrals', this%nneutrals )
  call input%get( 'simulation.nlasers', this%nlasers )
  call input%get( 'simulation.max_mode', this%max_mode )

  call write_stdout( 'Initializing fields...' )
  call this%fields%new( input, opts )

  call write_stdout( 'Initializing beams...' )
  call this%beams%new( input, opts )

  call write_stdout( 'Initializing lasers...' )
  call this%lasers%new( input, opts )

  call write_stdout( 'Initializing plasma...' )
  call this%plasma%new( input, opts, (this%start3d-1)*dt )

  call write_stdout( 'Initializing diagnostics...' )
  call this%diag%new( input, opts, this%fields, this%beams, this%plasma, this%lasers )

  call write_stdout( 'Initializing pipeline...' )
  allocate( this%tag_field(p_max_tag_num), this%id_field(p_max_tag_num) )
  allocate( this%tag_beam(this%nbeams), this%id_beam(this%nbeams) )
  allocate( this%tag_laser(2, this%nlasers), this%id_laser(2, this%nlasers) )
  allocate( this%tag_spe(this%nspecies), this%id_spe(this%nspecies) )
  allocate( this%tag_neut(4, this%nneutrals), this%id_neut(4, this%nneutrals) )
  allocate( this%tag_bq(this%nbeams), this%id_bq(this%nbeams) )

  this%id_field = MPI_REQUEST_NULL
  this%id_spe   = MPI_REQUEST_NULL
  this%id_neut  = MPI_REQUEST_NULL
  this%id_beam  = MPI_REQUEST_NULL
  this%id_laser = MPI_REQUEST_NULL
  this%id_bq    = MPI_REQUEST_NULL

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_simulation

subroutine end_simulation(this)

  implicit none

  class( simulation ), intent(inout) :: this

  ! local data
  integer :: ierr
  character(len=18), save :: sname = 'end_simulation'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call write_stdout( 'Terminating simulation...' )
  call this%fields%del()
  call this%beams%del()
  call this%lasers%del()
  call this%plasma%del()
  call this%diag%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

  call mpi_finalize(ierr)

end subroutine end_simulation

subroutine run_simulation( this )

  implicit none

  class( simulation ), intent(inout) :: this

  integer :: i, j, k, l, ierr
  real :: rel_res, abs_res
  integer, dimension(MPI_STATUS_SIZE) :: istat
  character(len=32), save :: sname = 'run_simulation'

  class(field_psi), pointer :: psi
  class(field_vpot), pointer :: vpot
  class(field_e), pointer :: e_spe, e_beam, e
  class(field_b), pointer :: b_spe, b_beam, b
  class(field), pointer :: chi
  class(field_jay), pointer :: cu, amu
  class(field_rho), pointer :: q_spe, q_beam
  class(field_djdxi), pointer :: dcu, acu
  class(field_laser), pointer :: laser_all
  class(beam3d), dimension(:), pointer :: beam
  class(field_laser), dimension(:), pointer :: laser
  class(species2d), dimension(:), pointer :: spe
  class(neutral), dimension(:), pointer :: neut

  call write_dbg( cls_name, sname, cls_level, 'starts' )

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

  beam      => this%beams%beam
  laser     => this%lasers%laser
  laser_all => this%lasers%laser_all
  chi       => this%lasers%chi
  spe       => this%plasma%spe
  neut      => this%plasma%neut

  ! deposit beams and do diagnostics to see the initial distribution if it is
  ! a fresh run
  call write_stdout( 'Starting simulation...' )
  if ( this%start3d == 1 ) then

    call q_beam%as(0.0)
    call q_spe%as(0.0)
    ! pipeline data transfer for beams
    do k = 1, this%nbeams
      this%tag_bq(k) = ntag()
      call beam(k)%qdp( q_beam, this%tag_bq(k), this%id_bq(k) )
    enddo

    call this%diag%run( 0, this%dt )

  endif

  call start_tprof( 'total simulation time' )

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

    ! pipeline data transfer for neutrals
    do k = 1, this%nneutrals
      ! tag 1 and 2 are for particle array and ion density transfer respectively
      this%tag_neut(1,k) = ntag()
      this%tag_neut(2,k) = ntag()
      this%tag_neut(3,k) = ntag()
      this%tag_neut(4,k) = ntag()
      call neut(k)%precv( this%tag_neut(1:4,k) )
    enddo

    b     = 0.0
    e     = 0.0
    b_spe = 0.0
    e_spe = 0.0
    psi   = 0.0
    cu    = 0.0
    acu   = 0.0
    amu   = 0.0
    do k = 1, this%nlasers
      call laser(k)%zero( only_f1=.true. )
    enddo

    ! pipeline data transfer for current and species B-field
    this%tag_field(1) = ntag()
    call cu%pipe_recv( this%tag_field(1), 'forward', 'replace' )
    this%tag_field(4) = ntag()
    call b_spe%pipe_recv( this%tag_field(4), 'forward', 'replace' )

    do j = 1, this%nstep2d

      call q_beam%copy_slice( j, p_copy_2to1 )
      call b_beam%solve( q_beam )
      q_spe = 0.0
      do k = 1, this%nspecies
        call spe(k)%qdp( q_spe )
      enddo

      do k = 1, this%nneutrals
        call neut(k)%qdp( q_spe )
        call neut(k)%ion_deposit( q_spe )
      enddo

      call psi%solve( q_spe )
      do k = 1, this%nspecies
        call spe(k)%interp_psi(psi)
      enddo
      call b_spe%solve( cu )

      call laser_all%zero( only_f1=.true. )
      do k = 1, this%nlasers
        call laser(k)%copy_slice( j, p_copy_2to1 )
        call laser(k)%set_grad( j )
        call laser_all%gather( laser(k) )
      enddo

      ! predictor-corrector iteration
      do l = 1, this%iter_max

        ! store the old Br
        call convergence_tester(b_spe, 2, 'record')

        call add_f1( b_spe, b_beam, b )
        call e%solve( cu )
        call e%solve( b, psi )
        cu = 0.0
        acu = 0.0
        amu = 0.0

        do k = 1, this%nspecies
          call spe(k)%amjdp( e, b, laser_all, cu, amu, acu, this%dxi )
        enddo

        do k = 1, this%nneutrals
          call neut(k)%amjdp( e, b, cu, amu, acu, this%dxi )
        enddo

        call dcu%solve( acu, amu )
        call b_spe%solve( dcu, cu )
        call b_spe%solve( cu )

        ! get the relative error between the old and new Br
        call convergence_tester(b_spe, 2, 'compare', rel_res=rel_res, abs_res=abs_res)
        if (rel_res < this%iter_reltol .or. abs_res < this%iter_abstol) exit

      enddo ! iteration

      ! deposit chi
      call this%lasers%deposit_chi( spe, j )

      do k = 1, this%nspecies
        call spe(k)%cbq(j)
      enddo
      do k = 1, this%nneutrals
        call neut(k)%cbq(j)
      enddo
      call cu%copy_slice( j, p_copy_1to2 )
      call add_f1( cu, q_spe, (/3/), (/1/) )
      call q_spe%copy_slice( j, p_copy_1to2 )

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

      ! send the last slice of current and species B-field to the next stage
      if ( j == this%nstep2d ) then
        call mpi_wait( this%id_field(1), istat, ierr )
        call cu%pipe_send( this%tag_field(1), this%id_field(1), 'forward' )
        call mpi_wait( this%id_field(4), istat, ierr )
        call b_spe%pipe_send( this%tag_field(4), this%id_field(4), 'forward' )
      endif

      ! advance species particles
      do k = 1, this%nspecies
        call spe(k)%push_u( e, b, laser_all, this%dxi )
        call spe(k)%push_x( this%dxi )
        ! call spe(k)%sort( this%start2d + j - 1 )
      enddo

      ! ionize and advance particles of neutrals
      do k = 1, this%nneutrals
        call neut(k)%update( e, psi, i*this%dt )
        call neut(k)%push_u( e, b, this%dxi )
        call neut(k)%push_x( this%dxi )
        ! call neut(k)%push( e, b, laser_all )
        ! TODO: add sorting
      enddo

      call e%copy_slice( j, p_copy_1to2 )
      call b%copy_slice( j, p_copy_1to2 )
      call psi%copy_slice( j, p_copy_1to2 )
      call b_spe%copy_slice( j, p_copy_1to2 )
      call e_spe%copy_slice( j, p_copy_1to2 )

      ! send the first slice of E and B field back to the last stage for 3D 
      ! particle push
      if ( j == 1 ) then
        call mpi_wait( this%id_field(2), istat, ierr )
        this%tag_field(2) = ntag()
        call b%pipe_send( this%tag_field(2), this%id_field(2), 'backward', 'inner' )
        call mpi_wait( this%id_field(3), istat, ierr )
        this%tag_field(3) = ntag()
        call e%pipe_send( this%tag_field(3), this%id_field(3), 'backward', 'inner' )
      endif

    enddo ! 2d loop

    ! pipeline for species
    do k = 1, this%nspecies
      call spe(k)%psend( this%tag_spe(k), this%id_spe(k) )
    enddo

    ! pipeline for neutrals
    do k = 1, this%nneutrals
      call neut(k)%psend( this%tag_neut(1:4,k), this%id_neut(1:4,k) )
    enddo

    ! pipeline for E and B fields
    call b%pipe_recv( this%tag_field(2), 'backward', 'guard', 'replace' )
    call e%pipe_recv( this%tag_field(3), 'backward', 'guard', 'replace' )

    ! advance laser fields
    call this%lasers%advance()

    ! pipeline for beams
    do k = 1, this%nbeams
      this%tag_beam(k) = ntag()
      call mpi_wait( this%id_beam(k), istat, ierr )
      call beam(k)%push( e, b, this%tag_beam(k), this%id_beam(k) )
    enddo

    call this%diag%run( this%tstep, this%dt )

    ! renew species for next 3D step
    do k = 1, this%nspecies
      call mpi_wait( this%id_spe(k), istat, ierr )
      call spe(k)%renew( i*this%dt )
    enddo

    ! renew neutrals for next 3D step
    do k = 1, this%nneutrals
      call mpi_wait( this%id_neut(1,k), istat, ierr )
      call mpi_wait( this%id_neut(2,k), istat, ierr )
      call mpi_wait( this%id_neut(3,k), istat, ierr )
      call mpi_wait( this%id_neut(4,k), istat, ierr )
      call neut(k)%renew( i*this%dt )
    enddo

  enddo ! 3d loop

  call stop_tprof( 'total simulation time' )

  call write_tprof()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine run_simulation

subroutine convergence_tester(fld, dim, operation, rel_res, abs_res)

  use ufield_class
  use field_class

  implicit none
  class(field), intent(in) :: fld
  integer, intent(in) :: dim
  character(len=*), intent(in) :: operation
  real, intent(out), optional :: rel_res, abs_res

  integer :: max_mode, mode, nrp, ierr
  real, dimension(:), allocatable, save :: fld_re_old, fld_im_old
  real :: fld_old_norm, norm_tmp, res
  type(ufield), dimension(:), pointer :: fld_re_ptr => null(), fld_im_ptr => null()

  max_mode = fld%get_max_mode()
  fld_re_ptr => fld%get_rf_re()
  if (max_mode > 0) fld_im_ptr => fld%get_rf_im()
  nrp = fld_re_ptr(0)%get_ndp(1)

  if (.not. allocated(fld_re_old)) then
    allocate(fld_re_old(nrp), fld_im_old(nrp))
  endif

  select case (trim(operation))

    case ("record")

      ! In the "record" mode, the subroutine stores the sum of field values of the given dimension.
      fld_re_old = 0.0
      fld_im_old = 0.0
      do mode = 0, max_mode
        fld_re_old = fld_re_old + abs(fld_re_ptr(mode)%f1(dim, 1:nrp))
        if (mode == 0) cycle
        fld_im_old = fld_im_old + abs(fld_im_ptr(mode)%f1(dim, 1:nrp))
      enddo

    case ("compare")

      ! In the "compare" mode, the subroutine compares the relative error between the old and new field values.
      if (.not. present(rel_res) .or. .not. present(abs_res)) then
        call write_err("Parameter 'rel_res' and 'abs_res' must be given for 'compare' operation.")
      endif

      ! 2-norm
      ! norm_tmp = norm2(fld_re_old)**2 + norm2(fld_im_old)**2
      ! call mpi_allreduce(norm_tmp, fld_old_norm, 1, p_dtype_real, MPI_SUM, comm_loc(), ierr)

      ! infinity-norm
      norm_tmp = sqrt(maxval(fld_re_old**2 + fld_im_old**2))
      call mpi_allreduce(norm_tmp, fld_old_norm, 1, p_dtype_real, MPI_MAX, comm_loc(), ierr)
      
      ! DEBUG
      ! if (id_proc_loc() == 0) then
      !   print *, "fld_old_norm = ", fld_old_norm
      ! endif

      do mode = 0, max_mode
        fld_re_old = fld_re_old - abs(fld_re_ptr(mode)%f1(dim, 1:nrp))
        if (mode == 0) cycle
        fld_im_old = fld_im_old - abs(fld_im_ptr(mode)%f1(dim, 1:nrp))
      enddo

      ! 2-norm
      ! norm_tmp = norm2(fld_re_old)**2 + norm2(fld_im_old)**2
      ! call mpi_allreduce(norm_tmp, res, 1, p_dtype_real, MPI_SUM, comm_loc(), ierr)
      ! residue = sqrt(res / fld_old_norm)

      ! infinity-norm
      norm_tmp = sqrt(maxval(fld_re_old**2 + fld_im_old**2))
      call mpi_allreduce(norm_tmp, abs_res, 1, p_dtype_real, MPI_MAX, comm_loc(), ierr)
      ! deal with divided-by-zero error
      if (fld_old_norm > epsilon(1.0d0)) then
        rel_res = abs_res / fld_old_norm
      else
        rel_res = huge(1.0)
      endif

    case default
      call write_err("Invalid operation mode!")

  end select

end subroutine convergence_tester

end module simulation_class