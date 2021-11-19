module sim_lasers_class

use parallel_module
use options_class
use sysutil_module
use param
use input_class
use field_laser_class
use kwargs_class

implicit none

private

public :: sim_lasers

type, public :: sim_lasers
  
  class( field_laser ), dimension(:), pointer :: laser => null()
  integer :: num_lasers

  contains

  procedure :: alloc => alloc_sim_lasers
  procedure :: new   => init_sim_lasers
  procedure :: del   => end_sim_lasers

end type sim_lasers

character(len=32), save :: cls_name = 'sim_lasers'
integer, save :: cls_level = 2

contains

subroutine alloc_sim_lasers( this, input )

  implicit none
  class( sim_lasers ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  ! local data
  integer :: max_mode, i
  character(len=32), save :: sname = 'alloc_sim_lasers'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.nlasers', this%num_lasers )
  call input%get( 'simulation.max_mode', max_mode )

  if ( .not. associated( this%laser ) ) then
    allocate( field_laser :: this%laser( this%num_lasers ) )
  endif

  do i = 1, this%num_lasers
    call this%laser(i)%alloc( max_mode )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_sim_lasers

subroutine init_sim_lasers( this, input, opts )

  implicit none

  class( sim_lasers ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts

  ! local data
  character(len=32), save :: sname = 'init_sim_lasers'
  character(len=:), allocatable :: str
  integer :: max_mode, iter, i
  integer, dimension(2,2) :: gc_num
  real :: k0 
  type( kw_list ) :: kwargs

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.max_mode', max_mode )
  call kwargs%append( 'iter', 0 )
  call kwargs%append( 'k0', 10.0 )

  ! WARNING: the following number of cell only support iter=1
  gc_num(:,1) = (/1, 1/)
  gc_num(:,2) = (/1, 1/)

  do i = 1, this%num_lasers
    call input%get( 'laser' // num2str(i) // '.iteration', iter )
    call input%get( 'laser' // num2str(i) // '.k0', k0 )
    call kwargs%set( 'iter', iter )
    call kwargs%set( 'k0', k0 )
    call this%laser(i)%new( opts, 1, max_mode, gc_num, only_f1=.false., kwargs=kwargs )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_lasers

subroutine end_sim_lasers( this )

  implicit none
  class(sim_lasers), intent(inout) :: this
  ! local data
  character(len=32), save :: sname = 'end_sim_lasers'
  integer :: i

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  do i = 1, this%num_lasers
    call this%laser(i)%del()
  enddo
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_lasers

end module sim_lasers_class