module sim_species_class

use species2d_class
use parallel_pipe_class
use grid_class
use fdist2d_class
use input_class
use param
use sys

implicit none

private

public :: sim_species

type sim_species

  ! private

  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()
  type( species2d ), dimension(:), allocatable :: spe
  type( fdist2d_wrap ), dimension(:), allocatable :: pf

  contains

  generic :: new => init_sim_species
  generic :: del => end_sim_species

  procedure, private :: init_sim_species, end_sim_species

end type sim_species

character(len=18), save :: cls_name = 'sim_species'
integer, save :: cls_level = 2

contains

subroutine init_sim_species( this, input, s )

  implicit none

  class( sim_species ), intent(inout) :: this
  type( input_json ), pointer, intent(inout) :: input
  real, intent(in) :: s

  ! local data
  character(len=18), save :: sname = 'init_sim_species'
  integer :: i, n
  ! real, dimension(3,100) :: arg
  ! logical :: quiet
  real :: qm, qbm
  integer :: ps, sm_type, sm_ord, max_mode, npf
  ! type(hdf5file) :: file_rst
  character(len=:), allocatable :: str

  this%gp => input%gp
  this%pp => input%pp
  ! pqb => qb

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
       call this%pf(i)%p%new( input, i )
    case (12)
       allocate( fdist2d_012 :: this%pf(i)%p )
       call this%pf(i)%p%new( input, i )
    ! Add new distributions right above this line
    case default
       call write_err( 'Invalid beam profile!' )
    end select

    call input%get('species('//num2str(i)//').q',qm)
    call input%get('species('//num2str(i)//').m',qbm)
    qbm = qm/qbm
    call this%spe(i)%new( this%pp, this%gp, this%pf(i)%p, ps, max_mode,&
      qbm, 8, s, sm_type, sm_ord )

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

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_species

end module sim_species_class