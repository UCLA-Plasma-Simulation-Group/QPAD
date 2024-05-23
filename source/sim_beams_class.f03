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

  integer :: num_beams

  contains

  procedure :: alloc => alloc_sim_beams
  procedure :: new   => init_sim_beams
  procedure :: del   => end_sim_beams

end type sim_beams

character(len=18), save :: cls_name = 'sim_beams'
integer, save :: cls_level = 2

contains

subroutine alloc_sim_beams( this, input, opts )

  implicit none

  class( sim_beams ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  ! local data
  integer :: i
  character(len=18), save :: sname = 'alloc_sim_beams'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.nbeams', this%num_beams )

  if ( .not. associated( this%beam ) ) then
    allocate( beam3d :: this%beam( this%num_beams ) )
  endif

  do i = 1, this%num_beams
    call this%beam(i)%alloc( input, opts, i )
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
  real :: dt
  logical :: read_rst
  integer :: rst_timestep, ps, sm_ord, ierr, max_mode
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
  case ( 'quadratic' )
    ps = p_ps_quadratic
    call write_stdout('Attention! Quadratic interpolation for beams not implemented!')
  case default
    call write_err( 'Invalid interpolation type! Only "linear" and "quadratic" are supported &
      &currently.' )
  end select

  ! initialize beam particle manager
  do i = 1, this%num_beams
    part_dim = p_x_dim + p_p_dim + 1
    if ( this%beam(i)%pf%has_spin ) part_dim = part_dim + p_s_dim
    call set_part3d_comm( part_dim, this%beam(i)%pf%npmax )
  enddo
  call init_part3d_comm( opts )

  ! initialize beams
  do i = 1, this%num_beams

    sm_ord = 0
    if(input%found('beams('//num2str(i)//').smooth_order')) then
      call input%get('beams('//num2str(i)//').smooth_order', sm_ord)
    endif

    call this%beam(i)%new( input, opts, max_mode, ps, dt, sm_ord, i )

    if ( read_rst ) then
      call input%get( 'simulation.restart_timestep', rst_timestep )
      call file_rst%new(&
        filename = './RST/Beam'//num2str(i,2)//'/', &
        dataname = 'RST-beam'//num2str(i,2)//'-'//num2str(id_proc(),6), &
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
