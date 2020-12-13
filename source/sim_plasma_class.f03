module sim_plasma_class

use species2d_class
use neutral_class
use parallel_module
use options_class
use fdist2d_class
use input_class
use param
use sysutil_module
use part2d_comm

implicit none

private

public :: sim_plasma

type sim_plasma

  class( species2d ), dimension(:), pointer :: spe => null()
  class( neutral ), dimension(:), pointer :: neut => null()
  type( fdist2d_wrap ), dimension(:), pointer :: pf_spe => null(), pf_neut => null()

  integer :: num_species, num_neutrals

  contains

  procedure :: alloc => alloc_sim_plasma
  procedure :: new   => init_sim_plasma
  procedure :: del   => end_sim_plasma

end type sim_plasma

character(len=18), save :: cls_name = 'sim_plasma'
integer, save :: cls_level = 2

contains

subroutine alloc_sim_plasma( this, input )

  implicit none

  class( sim_plasma ), intent(inout) :: this
  type( input_json ), intent(inout) :: input

  integer :: i

  call input%get( 'simulation.nspecies', this%num_species )
  call input%get( 'simulation.nneutrals', this%num_neutrals )

  if ( .not. associated( this%spe ) ) allocate( species2d :: this%spe( this%num_species ) )
  if ( .not. associated( this%neut ) ) allocate( neutral :: this%neut( this%num_neutrals ) )

  do i = 1, this%num_species
    call this%spe(i)%alloc()
  enddo

  do i = 1, this%num_neutrals
    call this%neut(i)%alloc()
  enddo

end subroutine alloc_sim_plasma

subroutine init_sim_plasma( this, input, opts, s )

  implicit none

  class( sim_plasma ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  real, intent(in) :: s

  ! local data
  character(len=18), save :: sname = 'init_sim_plasma'
  real :: qm, qbm, omega_p, np
  integer :: i, ps, sm_type, sm_ord, max_mode, npf, part_dim, elem, ion_max
  character(len=:), allocatable :: str

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.max_mode', max_mode )
  call input%get( 'simulation.n0', np )
  omega_p = sqrt(np) * 5.641460231180626d4

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

  allocate( this%pf_spe( this%num_species ), this%pf_neut( this%num_neutrals ) )

  do i = 1, this%num_species

    call input%get( 'species('//num2str(i)//').profile', npf )
    select case ( npf )
    case (0)
       allocate( fdist2d_000 :: this%pf_spe(i)%p )
       call this%pf_spe(i)%p%new( input, opts, i, 'species' )
    case (12)
       allocate( fdist2d_012 :: this%pf_spe(i)%p )
       call this%pf_spe(i)%p%new( input, opts, i, 'species' )
    ! Add new distributions right above this line
    case default
       call write_err( 'Invalid species profile!' )
    end select

  enddo

  do i = 1, this%num_neutrals

    call input%get( 'neutrals('//num2str(i)//').profile', npf )
    select case ( npf )
    case (0)
       allocate( fdist2d_000 :: this%pf_neut(i)%p )
       call this%pf_neut(i)%p%new( input, opts, i, 'neutrals' )
    case (12)
       allocate( fdist2d_012 :: this%pf_neut(i)%p )
       call this%pf_neut(i)%p%new( input, opts, i, 'neutrals' )
    ! Add new distributions right above this line
    case default
       call write_err( 'Invalid neutral profile!' )
    end select

  enddo

  ! initialize 2D particle manager
  part_dim = 8
  ! loop over all the 2D particle profile to get the buffer size
  do i = 1, this%num_species
    call set_part2d_comm( part_dim, npmax = this%pf_spe(i)%p%getnpmax() )
  enddo
  do i = 1, this%num_neutrals
    call set_part2d_comm( part_dim, npmax = this%pf_neut(i)%p%getnpmax() )
  enddo
  call init_part2d_comm( opts )

  do i = 1, this%num_species

    call input%get( 'species('//num2str(i)//').q', qm )
    call input%get( 'species('//num2str(i)//').m', qbm )
    qbm = qm / qbm
    call this%spe(i)%new( opts, this%pf_spe(i)%p, ps, max_mode, &
      qbm, s, sm_type, sm_ord )

  enddo

  do i = 1, this%num_neutrals

    call input%get( 'neutrals('//num2str(i)//').q', qm )
    call input%get( 'neutrals('//num2str(i)//').m', qbm )
    call input%get( 'neutrals('//num2str(i)//').element', elem )
    call input%get( 'neutrals('//num2str(i)//').ion_max', ion_max )
    qbm = qm / qbm
    call this%neut(i)%new( opts, this%pf_neut(i)%p, max_mode, elem, ion_max, &
      qbm, omega_p, s, sm_type, sm_ord )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_plasma

subroutine end_sim_plasma( this )

  implicit none

  class( sim_plasma ), intent(inout) :: this

  integer :: i
  character(len=18), save :: sname = 'end_sim_plasma'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 1, this%num_species
    call this%spe(i)%del()
  enddo

  do i = 1, this%num_neutrals
    call this%neut(i)%del()
  enddo

  call end_part2d_comm()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_plasma

end module sim_plasma_class