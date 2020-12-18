module sim_beams_class

use beam3d_class
use parallel_module
use options_class
use fdist3d_class
use input_class
use field_class
use field_src_class
use hdf5io_class
use param
use sysutil_module
use part3d_comm

implicit none

private

public :: sim_beams

type sim_beams

  ! private

  class( beam3d ), dimension(:), pointer :: beam => null()
  type( fdist3d_wrap ), dimension(:), pointer :: pf => null()

  integer :: num_beams

  contains

  procedure :: alloc => alloc_sim_beams
  procedure :: new   => init_sim_beams
  procedure :: del   => end_sim_beams

end type sim_beams

character(len=18), save :: cls_name = 'sim_beams'
integer, save :: cls_level = 2

contains

subroutine alloc_sim_beams( this, input )

  implicit none

  class( sim_beams ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  ! local data
  integer :: i
  character(len=18), save :: sname = 'alloc_sim_beams'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.nbeams', this%num_beams )

  if ( .not. associated( this%beam ) ) then
    allocate( beam3d :: this%beam( this%num_beams ) )
  endif

  do i = 1, this%num_beams
    call this%beam(i)%alloc()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_sim_beams

subroutine init_sim_beams( this, input, opts )

  implicit none

  class( sim_beams ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts

  ! local data
  character(len=18), save :: sname = 'init_sim_beams'
  integer :: i, part_dim
  real :: qm, qbm, dt, amm = 0.0
  logical :: read_rst, has_spin
  integer :: rst_timestep, ps, sm_type, sm_ord, ierr, max_mode, npf, push_type
  type(hdf5file) :: file_rst
  character(len=:), allocatable :: str

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.dt', dt )
  call input%get( 'simulation.read_restart', read_rst )
  ! call input%get( 'simulation.nbeams', n )
  call input%get( 'simulation.max_mode', max_mode )

  ! read interpolation type
  call input%get( 'simulation.interpolation', str )
  select case ( trim(str) )
  case ( 'linear' )
    ps = p_ps_linear
  case default
    call write_err( 'Invalid interpolation type! Only "linear" are supported currently.' )
  end select

  ! read smooth parameters
  call input%get( 'simulation.smooth_type', str )
  select case ( trim(str) )
  case ( 'none' )
    sm_type = p_smooth_none
  case ( 'binomial' )
    sm_type = p_smooth_binomial
  case ( 'compensated' )
    sm_type = p_smooth_compensated
  case default
    call write_err( 'Invalid smooth type! Only "binomial" and "compensated" are supported currently.' )
  end select
  call input%get( 'simulation.smooth_order', sm_ord )

  ! allocate( this%beam(n), this%pf(n) )
  allocate( this%pf( this%num_beams ) )

  ! initialize beam profile for each beam
  do i = 1, this%num_beams

    call input%get( 'beam('//num2str(i)//').profile', npf )
    select case ( npf )
    case (0)
       allocate( fdist3d_000 :: this%pf(i)%p )
       call this%pf(i)%p%new( input, i )
    case (1)
       allocate( fdist3d_001 :: this%pf(i)%p )
       call this%pf(i)%p%new( input, i )
    case (2)
       allocate( fdist3d_002 :: this%pf(i)%p )
       call this%pf(i)%p%new( input, i )
    case (100)
       allocate( fdist3d_100 :: this%pf(i)%p )
       call this%pf(i)%p%new( input, i )
  ! Add new distributions right above this line
    case default
      call write_err( 'Invalid beam profile!' )
    end select

  enddo

  ! initialize beam particle manager
  do i = 1, this%num_beams
    part_dim = p_x_dim + p_p_dim + 1
    call input%get( 'beam('//num2str(i)//').has_spin', has_spin )
    if ( has_spin ) part_dim = part_dim + p_s_dim
    call set_part3d_comm( part_dim, this%pf(i)%p%getnpmax() )
  enddo
  call init_part3d_comm( opts )

  ! initialize beams
  do i = 1, this%num_beams

    call input%get( 'beam('//num2str(i)//').q', qm )
    call input%get( 'beam('//num2str(i)//').m', qbm )
    qbm = qm/qbm

    push_type = p_push_reduced
    call input%get( 'beam('//num2str(i)//').push_type', str )
    select case ( trim(str) )
    case ( 'reduced' )
      push_type = p_push_reduced
    case ( 'boris' )
      push_type = p_push_boris
    case default
      call write_err( 'Invalid pusher type! Only "reduced" and "boris" are supported currently.' )
    end select

    call input%get( 'beam('//num2str(i)//').has_spin', has_spin )
    if (has_spin) then
      call input%get( 'beam('//num2str(i)//').anom_mag_moment', amm )
    endif

    call this%beam(i)%new( opts, max_mode, ps, this%pf(i)%p, &
      qbm, dt, push_type, sm_type, sm_ord, has_spin, amm )

    if ( read_rst ) then
      call input%get( 'simulation.restart_timestep', rst_timestep )
      call file_rst%new(&
        filename = './RST/Beam'//num2str(i,2)//'/',&
        dataname = 'RST-beam'//num2str(i,2)//'-'//num2str(id_proc(),6),&
        n = rst_timestep)
      call this%beam(i)%rrst(file_rst)
    endif

  enddo

  call MPI_BARRIER( comm_world(), ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_beams

subroutine end_sim_beams( this )

  implicit none

  class( sim_beams ), intent(inout) :: this

  integer :: i, n
  character(len=18), save :: sname = 'end_sim_beams'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  n = size( this%beam )

  do i = 1, n
    call this%beam(i)%del()
  enddo

  call end_part3d_comm()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_beams

end module sim_beams_class