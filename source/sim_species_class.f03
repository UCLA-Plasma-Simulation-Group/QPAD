module sim_species_class

use species2d_class
use parallel_module
use options_class
use fdist2d_class
use input_class
use param
use sysutil
use part2d_comm

implicit none

private

public :: sim_species

type sim_species

  ! private

  type( species2d ), dimension(:), pointer :: spe => null()
  type( fdist2d_wrap ), dimension(:), pointer :: pf => null()

  contains

  generic :: new => init_sim_species
  generic :: del => end_sim_species

  procedure, private :: init_sim_species, end_sim_species

end type sim_species

character(len=18), save :: cls_name = 'sim_species'
integer, save :: cls_level = 2

contains

subroutine init_sim_species( this, input, opts, s )

  implicit none

  class( sim_species ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  real, intent(in) :: s

  ! local data
  character(len=18), save :: sname = 'init_sim_species'
  integer :: i, n
  ! real, dimension(3,100) :: arg
  ! logical :: quiet
  real :: qm, qbm
  integer :: ps, sm_type, sm_ord, max_mode, npf, part_dim
  ! type(hdf5file) :: file_rst
  character(len=:), allocatable :: str

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.nspecies', n )
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

  allocate( this%spe(n), this%pf(n) )

  do i = 1, n

    call input%get( 'species('//num2str(i)//').profile', npf )
    select case ( npf )
    case (0)
       allocate( fdist2d_000 :: this%pf(i)%p )
       call this%pf(i)%p%new( input, opts, i )
    case (12)
       allocate( fdist2d_012 :: this%pf(i)%p )
       call this%pf(i)%p%new( input, opts, i )
    ! Add new distributions right above this line
    case default
       call write_err( 'Invalid beam profile!' )
    end select

  enddo

  ! initialize species particle manager
  part_dim = 8
  do i = 1, n
    call set_part2d_comm( part_dim, npmax = this%pf(i)%p%getnpmax() )
  enddo
  call init_part2d_comm( opts )

  do i = 1, n

    call input%get('species('//num2str(i)//').q',qm)
    call input%get('species('//num2str(i)//').m',qbm)
    qbm = qm/qbm
    call this%spe(i)%new( opts, this%pf(i)%p, ps, max_mode,&
      qbm, part_dim, s, sm_type, sm_ord )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_species

subroutine end_sim_species( this )

  implicit none

  class( sim_species ), intent(inout) :: this

  integer :: i, n
  character(len=18), save :: sname = 'end_sim_species'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  n = size( this%spe )

  do i = 1, n
    call this%spe(i)%del()
  enddo

  call end_part2d_comm()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_species

end module sim_species_class