module sim_beams_class

use beam3d_class
use parallel_pipe_class
use grid_class
use fdist3d_class
use input_class
use field_class
use field_src_class
use hdf5io_class
use param
use sysutil

implicit none

private

public :: sim_beams

type sim_beams

  ! private

  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()
  type( beam3d ), dimension(:), pointer :: beam => null()
  type( fdist3d_wrap ), dimension(:), pointer :: pf => null()

  contains

  generic :: new => init_sim_beams
  generic :: del => end_sim_beams

  procedure, private :: init_sim_beams, end_sim_beams

end type sim_beams

character(len=18), save :: cls_name = 'sim_beams'
integer, save :: cls_level = 2

contains

subroutine init_sim_beams( this, input )

  implicit none

  class( sim_beams ), intent(inout) :: this
  type( input_json ), pointer, intent(inout) :: input

  ! local data
  character(len=18), save :: sname = 'init_sim_beams'
  integer :: i, n
  ! real, dimension(3,100) :: arg
  ! logical :: quiet
  real :: qm, qbm, dt
  logical :: read_rst
  integer :: rst_timestep, ps, sm_type, sm_ord, ierr, max_mode, npf
  type(hdf5file) :: file_rst
  character(len=:), allocatable :: str

  this%gp => input%gp
  this%pp => input%pp

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.dt', dt )
  call input%get( 'simulation.read_restart', read_rst )
  call input%get( 'simulation.nbeams', n )
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

  allocate( this%beam(n), this%pf(n) )

  do i = 1, n

    call input%get( 'beam('//num2str(i)//').profile', npf )
    select case ( npf )
    case (0)
       allocate( fdist3d_000 :: this%pf(i)%p )
       call this%pf(i)%p%new( input, i )
    case (1)
       allocate( fdist3d_001 :: this%pf(i)%p )
       call this%pf(i)%p%new( input, i )
  ! Add new distributions right above this line
    case default
      call write_err( 'Invalid beam profile!' )
    end select

    call input%get( 'beam('//num2str(i)//').q', qm )
    call input%get( 'beam('//num2str(i)//').m', qbm )
    qbm = qm/qbm

    call this%beam(i)%new( this%pp, this%gp, max_mode, ps, this%pf(i)%p, &
      qbm, dt, 7, sm_type, sm_ord )

    if ( read_rst ) then
      call input%get( 'simulation.restart_timestep', rst_timestep )
      call file_rst%new(&
        filename = './RST/Beam'//num2str(i,2)//'/',&
        dataname = 'RST-beam'//num2str(i,2)//'-'//num2str(this%pp%getidproc(),6),&
        n = rst_timestep)
      call this%beam(i)%rrst(file_rst)
    endif

  enddo

  call MPI_BARRIER( this%pp%getlworld(), ierr )

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

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_beams

end module sim_beams_class