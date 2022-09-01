module simulation_subcyc_class

use simulation_class
use parallel_module
use options_class
use sim_fields_class
use sim_beams_class
use sim_plasma_subcyc_class
use sim_lasers_class
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
use species2d_subcyc_class
use neutral_class
use neutral_subcyc_class

use input_class
use sysutil_module
use param
use mpi

use debug_tool

implicit none

private

public :: simulation_subcyc

integer, parameter :: p_max_tag_num = 32

type, extends(simulation) :: simulation_subcyc

  ! private

  ! maximum expansion factor of time step
  real :: exp_fac_max
  ! clamped expansion factor
  real :: exp_fac_clamped
  ! minimum 2D time step
  real :: dt_2d_min

  contains

  procedure :: alloc => alloc_simulation_subcyc
  procedure :: new   => init_simulation_subcyc
  procedure :: run   => run_simulation_subcyc

end type simulation_subcyc

character(len=18), save :: cls_name = 'simulation_subcyc'
integer, save :: cls_level = 1

contains

subroutine alloc_simulation_subcyc(this, input, opts)

  implicit none

  class(simulation_subcyc), intent(inout) :: this
  type(input_json), intent(inout) :: input
  type(options), intent(in) :: opts
  ! local data
  character(len=18), save :: sname = 'alloc_simulation_subcyc'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  if (.not. associated(this%plasma)) allocate(sim_plasma_subcyc :: this%plasma)
  call this%simulation%alloc(input, opts)
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine alloc_simulation_subcyc

subroutine init_simulation_subcyc(this, input, opts)

  implicit none
  class(simulation_subcyc), intent(inout) :: this
  type(input_json), intent(inout) :: input
  type(options), intent(in) :: opts
  ! local data
  character(len=18), save :: sname = 'init_simulation_subcyc'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call input%get('simulation.dt_2d_min', this%dt_2d_min)
  call input%get('simulation.expansion_fac_max', this%exp_fac_max)
  call input%get('simulation.expansion_fac_clamped', this%exp_fac_clamped)
  call this%simulation%new(input, opts)
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_simulation_subcyc

subroutine run_simulation_subcyc(this)

  implicit none
  class(simulation_subcyc), intent(inout) :: this

  integer :: i, j, k, l, ierr, n_subcyc, i_subcyc
  real :: rel_res, abs_res, exp_fac_max, dxi_subcyc, exp_fac_max_3dloop
  integer, dimension(MPI_STATUS_SIZE) :: istat
  character(len=32), save :: sname = 'run_simulation_subcyc'

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

  call write_dbg(cls_name, sname, cls_level, 'starts')

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
  call write_stdout('Starting simulation...')
  if (this%start3d == 1) then

    call q_beam%as(0.0)
    call q_spe%as(0.0)
    ! pipeline data transfer for beams
    do k = 1, this%nbeams
      this%tag_bq(k) = ntag()
      call beam(k)%qdp(q_beam, this%tag_bq(k), this%id_bq(k))
    enddo

    call this%diag%run(0, this%dt)

  endif

  call start_tprof('total simulation time')

  do i = this%start3d, this%nstep3d

    this%tstep = i
    call write_stdout('3D step = '//num2str(i))

    call q_beam%as(0.0)
    call q_spe%as(0.0)

    ! pipeline data transfer for beams
    do k = 1, this%nbeams
      this%tag_bq(k) = ntag()
      call beam(k)%qdp(q_beam, this%tag_bq(k), this%id_bq(k))
    enddo

    ! pipeline data transfer for species
    do k = 1, this%nspecies
      this%tag_spe(k) = ntag()
      call spe(k)%precv(this%tag_spe(k))
    enddo

    ! pipeline data transfer for neutrals
    do k = 1, this%nneutrals
      ! tag 1 and 2 are for particle array and ion density transfer respectively
      this%tag_neut(1,k) = ntag()
      this%tag_neut(2,k) = ntag()
      this%tag_neut(3,k) = ntag()
      this%tag_neut(4,k) = ntag()
      call neut(k)%precv(this%tag_neut(1:4,k))
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
      call laser(k)%zero(only_f1=.true.)
    enddo

    ! pipeline data transfer for current and species B-field
    this%tag_field(1) = ntag()
    call cu%pipe_recv(this%tag_field(1), 'forward', 'replace')
    this%tag_field(4) = ntag()
    call b_spe%pipe_recv(this%tag_field(4), 'forward', 'replace')

    ! DEBUG
    exp_fac_max_3dloop = 0.0

    do j = 1, this%nstep2d

      call q_beam%copy_slice(j, p_copy_2to1)
      call b_beam%solve(q_beam)

      call laser_all%zero(only_f1=.true.)
      do k = 1, this%nlasers
        call laser(k)%copy_slice(j, p_copy_2to1)
        call laser(k)%set_grad(j)
        call laser_all%gather(laser(k))
      enddo

      ! get the maximum 2D time step expansion factor at last 2D step
      select type (obj => this%plasma)
        type is (sim_plasma_subcyc)
          call obj%get_exp_fac_max(exp_fac_max)
          ! DEBUG
          ! exp_fac_max_3dloop = max(exp_fac_max_3dloop, exp_fac_max)
      end select
      call get_subcyc_step(exp_fac_max, this%exp_fac_max, this%dxi, &
        this%dt_2d_min, dxi_subcyc, n_subcyc)

      ! begin subcycling
      do i_subcyc = 1, n_subcyc

        q_spe = 0.0
        do k = 1, this%nspecies
          call spe(k)%qdp(q_spe)
        enddo

        do k = 1, this%nneutrals
          call neut(k)%qdp(q_spe)
          call neut(k)%ion_deposit(q_spe)
        enddo

        call psi%solve(q_spe)
        do k = 1, this%nspecies
          call spe(k)%interp_psi(psi)
        enddo
        call b_spe%solve(cu)

        ! predictor-corrector iteration
        do l = 1, this%iter_max

          ! store the old Br
          call convergence_tester(b_spe, 2, 'record')

          call add_f1(b_spe, b_beam, b)
          call e%solve(cu)
          call e%solve(b, psi)
          cu = 0.0
          acu = 0.0
          amu = 0.0

          do k = 1, this%nspecies
            call spe(k)%amjdp(e, b, laser_all, cu, amu, acu, dxi_subcyc)
          enddo

          do k = 1, this%nneutrals
            call neut(k)%amjdp(e, b, cu, amu, acu, dxi_subcyc)
          enddo

          ! select type (obj => this%plasma)
          !   type is (sim_plasma_subcyc)
          !     call obj%clamp_exp_fac(this%exp_fac_clamped)
          ! end select

          call dcu%solve(acu, amu)
          call b_spe%solve(dcu, cu)
          call b_spe%solve(cu)

          ! get the relative error between the old and new Br
          call convergence_tester(b_spe, 2, 'compare', rel_res=rel_res, abs_res=abs_res)
          if (rel_res < this%iter_reltol .or. abs_res < this%iter_abstol) exit

        enddo ! iteration

        call add_f1(b_spe, b_beam, b)
        call e_spe%solve(b_spe, psi)
        call e%solve(cu)
        call e%solve(b, psi)

        ! advance species particles
        do k = 1, this%nspecies
          call spe(k)%push_u(e, b, laser_all, dxi_subcyc)
          select type (obj => spe(k))
            type is (species2d_subcyc)
              call obj%clamp_exp_fac(this%exp_fac_clamped)
              ! DEBUG
              ! call obj%get_exp_fac_max(exp_fac_max)
              ! exp_fac_max_3dloop = max(exp_fac_max_3dloop, exp_fac_max)
          end select
          call spe(k)%push_x(dxi_subcyc)
          ! call spe(k)%sort(this%start2d + j - 1)
        enddo

        ! ionize and advance particles of neutrals
        do k = 1, this%nneutrals
          call neut(k)%update(e, psi, i*this%dt)
          call neut(k)%push_u(e, b, dxi_subcyc)
          select type (obj => neut(k))
            type is (neutral_subcyc)
              call obj%clamp_exp_fac(this%exp_fac_clamped)
          end select
          call neut(k)%push_x(dxi_subcyc)
          ! call neut(k)%push(e, b, laser_all)
          ! TODO: add sorting
        enddo

      enddo ! subcycling

      ! deposit chi
      call this%lasers%deposit_chi(spe, j)

      do k = 1, this%nspecies
        call spe(k)%cbq(j)
      enddo
      do k = 1, this%nneutrals
        call neut(k)%cbq(j)
      enddo
      call cu%copy_slice(j, p_copy_1to2)
      call add_f1(cu, q_spe, (/3/), (/1/))
      call q_spe%copy_slice(j, p_copy_1to2)

      ! for vector potential diagnostics
      if (this%diag%has_vpotz .or. this%diag%has_vpott) then
        if (this%diag%has_vpotz) call vpot%solve_vpotz(cu)
        if (this%diag%has_vpott) call vpot%solve_vpott(cu)
        call vpot%copy_slice(j, p_copy_1to2)
      endif

      call dot_f1(this%dxi, dcu)
      call add_f1(dcu, cu, (/1,2/), (/1,2/))

      ! send the last slice of current and species B-field to the next stage
      if (j == this%nstep2d) then
        call mpi_wait(this%id_field(1), istat, ierr)
        call cu%pipe_send(this%tag_field(1), this%id_field(1), 'forward')
        call mpi_wait(this%id_field(4), istat, ierr)
        call b_spe%pipe_send(this%tag_field(4), this%id_field(4), 'forward')
      endif

      call e%copy_slice(j, p_copy_1to2)
      call b%copy_slice(j, p_copy_1to2)
      call psi%copy_slice(j, p_copy_1to2)
      call b_spe%copy_slice(j, p_copy_1to2)
      call e_spe%copy_slice(j, p_copy_1to2)

      ! send the first slice of E and B field back to the last stage for 3D 
      ! particle push
      if (j == 1) then
        call mpi_wait(this%id_field(2), istat, ierr)
        this%tag_field(2) = ntag()
        call b%pipe_send(this%tag_field(2), this%id_field(2), 'backward', 'inner')
        call mpi_wait(this%id_field(3), istat, ierr)
        this%tag_field(3) = ntag()
        call e%pipe_send(this%tag_field(3), this%id_field(3), 'backward', 'inner')
      endif

    enddo ! 2d loop

    ! DEBUG
    ! call write_stdout("fac = "//num2str(exp_fac_max_3dloop), only_root=.false.)

    ! pipeline for species
    do k = 1, this%nspecies
      call spe(k)%psend(this%tag_spe(k), this%id_spe(k))
    enddo

    ! pipeline for neutrals
    do k = 1, this%nneutrals
      call neut(k)%psend(this%tag_neut(1:4,k), this%id_neut(1:4,k))
    enddo

    ! pipeline for E and B fields
    call b%pipe_recv(this%tag_field(2), 'backward', 'guard', 'replace')
    call e%pipe_recv(this%tag_field(3), 'backward', 'guard', 'replace')

    ! advance laser fields
    call this%lasers%advance()

    ! pipeline for beams
    do k = 1, this%nbeams
      this%tag_beam(k) = ntag()
      call mpi_wait(this%id_beam(k), istat, ierr)
      call beam(k)%push(e, b, this%tag_beam(k), this%id_beam(k))
    enddo

    call this%diag%run(this%tstep, this%dt)

    ! renew species for next 3D step
    do k = 1, this%nspecies
      call mpi_wait(this%id_spe(k), istat, ierr)
      call spe(k)%renew(i*this%dt)
    enddo

    ! renew neutrals for next 3D step
    do k = 1, this%nneutrals
      call mpi_wait(this%id_neut(1,k), istat, ierr)
      call mpi_wait(this%id_neut(2,k), istat, ierr)
      call mpi_wait(this%id_neut(3,k), istat, ierr)
      call mpi_wait(this%id_neut(4,k), istat, ierr)
      call neut(k)%renew(i*this%dt)
    enddo

  enddo ! 3d loop

  call stop_tprof('total simulation time')

  call write_tprof()

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine run_simulation_subcyc

subroutine get_subcyc_step(exp_fac, exp_fac_max, dt, dt_min, dt_subcyc, n_subcyc)

  implicit none
  real, intent(in) :: exp_fac, exp_fac_max, dt, dt_min
  real, intent(out) :: dt_subcyc
  integer, intent(out) :: n_subcyc

  if (exp_fac > exp_fac_max) then
    n_subcyc = ceiling(exp_fac / exp_fac_max)
    dt_subcyc = dt / n_subcyc
    if (dt_subcyc < dt_min) then
      n_subcyc = floor(dt / dt_min)
      dt_subcyc = dt / n_subcyc
    endif
  else
    n_subcyc = 1
    dt_subcyc = dt
  endif

end subroutine get_subcyc_step

end module simulation_subcyc_class