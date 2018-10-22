module simulation_class

use parallel_class
use parallel_pipe_class
use grid_class
use field_psi_class
use field_e_class
use field_b_class
use field_src_class

use input_class
use system
use param
use mpi

implicit none

private

public :: simulation

type sim_fields

  private

  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()

  type( field_psi ), allocatable :: psi
  type( field_b ), allocatable :: b_spe, b_beam
  type( field_e ), allocatable :: e_spe, e_beam
  type( field_jay ), allocatable :: jay
  type( field_rho ), allocatable :: q_spe, q_beam
  type( field_djdxi ), allocatable :: djdxi

  contains

  generic :: new => init_sim_fields
  generic :: del => end_sim_fields

  procedure, private :: init_sim_fields, end_sim_fields

end type sim_fields

type simulation

  private

  type( input_json ), pointer :: input => null()
  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()

  type( sim_fields ) :: fields
  real :: dr, dxi, dt
  integer :: iter, nstep3d, nstep2d, start3d, nbeams, nspecies, tstep
  integer :: num_modes, interp

  contains

  generic :: new => init_simulation
  generic :: del => end_simulation
  ! generic :: go => go_simulation

  procedure, private :: init_simulation, end_simulation
  ! procedure, private :: init_diag, diag_simulation
  ! procedure, private :: go_simulation

end type simulation

contains

subroutine init_sim_fields( this, input, dr, dxi, num_modes, part_shape )

  implicit none

  class( sim_fields ), intent(inout) :: this
  type( input_json ), pointer, intent(inout) :: input
  real, intent(in) :: dr, dxi
  integer, intent(in) :: num_modes, part_shape

  ! local data
  character(len=18), save :: sname = 'init_sim_fields'
  character(len=18), save :: cls_name = 'sim_fields'
  integer, save :: cls_level = 1
  character(len=20) :: s1, s2, s3
  character(len=:), allocatable :: ff
  integer :: i,n,ndump,j,k,l,m
  integer :: entity

  this%gp => input%gp
  this%pp => input%pp

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  allocate( this%psi, this%e_spe, this%b_spe, this%e_beam, this%b_beam, &
    this%jay, this%q_spe, this%q_beam, this%djdxi )

  call this%psi%new( this%pp, this%gp, dr, dxi, num_modes, part_shape )
  call this%jay%new( this%pp, this%gp, dr, dxi, num_modes, part_shape )
  call this%q_spe%new( this%pp, this%gp, dr, dxi, num_modes, part_shape )
  call this%q_beam%new( this%pp, this%gp, dr, dxi, num_modes, part_shape )
  call this%djdxi%new( this%pp, this%gp, dr, dxi, num_modes, part_shape )
  entity = p_entity_plasma
  call this%e_spe%new( this%pp, this%gp, dr, dxi, num_modes, part_shape, entity )
  call this%b_spe%new( this%pp, this%gp, dr, dxi, num_modes, part_shape, entity )
  entity = p_entity_beam
  call this%e_beam%new( this%pp, this%gp, dr, dxi, num_modes, part_shape, entity )
  call this%b_beam%new( this%pp, this%gp, dr, dxi, num_modes, part_shape, entity )

  ! call input%get('simulation.nspecies',n)

  ! loop1: do i = 1, n
  !   write (s1, '(I4.4)') i
  !   call input%info('species('//trim(s1)//').diag',n_children=m)
  !   do j = 1, m
  !      write (s2, '(I4.4)') j
  !      call input%get('species('//trim(s1)//').diag'//'('//trim(s2)//').ndump',ndump)
  !      if (ndump>0) then
  !         call input%info('species('//trim(s1)//').diag'//'('//trim(s2)//').name',n_children=l)
  !         do k = 1, l
  !            write (s3, '(I4.4)') k
  !            if(allocated(ff)) deallocate(ff)
  !            call input%get('species('//trim(s1)//').diag'//'('//trim(s2)//').name'&
  !            &//'('//trim(s3)//')',ff)
  !            if (ff == 'jx' .or. ff == 'jy' .or. ff == 'jz') then
  !               allocate(this%cu3d)
  !               call this%cu3d%new(this%p,this%err,this%sp3,dim=1)
  !               exit loop1
  !            end if
  !         end do
  !      end if
  !   end do
  ! end do loop1

  ! call input%info('field.diag',n_children=n)

  ! loop2: do i = 1, n
  !   write (s1,'(I4.4)') i
  !   call input%get('field.diag('//trim(s1)//').ndump',ndump)
  !   if (ndump > 0) then
  !      call input%info('field.diag('//trim(s1)//').name',n_children=m)
  !      do j = 1, m
  !         write (s2,'(I4.4)') j
  !         if(allocated(ff)) deallocate(ff)
  !         call input%get('field.diag('//trim(s1)//').name('//trim(s2)//')',ff)
  !         if (ff == 'psi') then
  !            allocate(this%psi3d)
  !            call this%psi3d%new(this%p,this%err,this%sp3,dim=1)
  !            exit loop2
  !         end if
  !      end do
  !   end if
  ! end do loop2

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_fields

subroutine end_sim_fields( this )

  implicit none

  class(sim_fields), intent(inout) :: this
  ! local data
  character(len=18), save :: sname = 'end_sim_fields'
  character(len=18), save :: cls_name = 'sim_fields'
  integer, save :: cls_level = 1
  integer :: i, n

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%psi%del()
  call this%e_spe%del()
  call this%b_spe%del()
  call this%e_beam%del()
  call this%b_beam%del()
  call this%jay%del()
  call this%q_spe%del()
  call this%q_beam%del()
  call this%djdxi%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_fields

subroutine init_simulation(this)

  implicit none

  class(simulation), intent(inout) :: this
  ! local data
  character(len=18), save :: sname = 'init_simulation:'
  character(len=18), save :: cls_name = 'simulation'
  integer, save :: cls_level = 0

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

  call this%fields%new( this%input, this%dr, this%dxi, this%num_modes, this%interp )
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
  integer, save :: cls_level = 0
  integer :: ierr

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%fields%del()
  ! call this%beams%del()
  ! call this%species%del()
  call this%gp%del()
  call this%pp%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_simulation

end module simulation_class