module sim_plasma_class

use species2d_class
use neutral_class
use neutral2_class
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
  class( neutral2 ), dimension(:), pointer :: neut2 => null()
  type( fdist2d ), dimension(:), pointer :: pf_spe => null(), pf_neut => null(),pf_neut2 => null()

  integer :: num_species, num_neutrals, num_neutral2s

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
  call input%get( 'simulation.nneutral2s', this%num_neutral2s )

  if ( .not. associated( this%spe ) ) allocate( species2d :: this%spe( this%num_species ) )
  if ( .not. associated( this%neut ) ) allocate( neutral :: this%neut( this%num_neutrals ) )
  if ( .not. associated( this%neut2 ) ) allocate( neutral2 :: this%neut2( this%num_neutral2s ) )

  do i = 1, this%num_species
    call this%spe(i)%alloc()
  enddo

  do i = 1, this%num_neutrals
    call this%neut(i)%alloc()
  enddo

!   do i = 1, this%num_neutral2s
!     call this%neut2(i)%alloc()
!   enddo

end subroutine alloc_sim_plasma

subroutine init_sim_plasma( this, input, opts, s )

  implicit none

  class( sim_plasma ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  real, intent(in) :: s

  ! local data
  character(len=18), save :: sname = 'init_sim_plasma'
  real :: qm, qbme, qbm, qbmi, omega_p, np
  integer :: i, ps, sm_type, sm_ord, max_mode, npf, part_dim, elem, ion_max, push_type, v, sec
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

  ! initialize profiles of species and neutrals
  allocate( this%pf_spe( this%num_species ), this%pf_neut( this%num_neutrals ), this%pf_neut2( this%num_neutral2s) )

  do i = 1, this%num_species
    call this%pf_spe(i)%new( input, opts, 'species', i )
  enddo

  do i = 1, this%num_neutrals
    call this%pf_neut(i)%new( input, opts, 'neutrals', i )
  enddo

  do i = 1, this%num_neutral2s
    call this%pf_neut2(i)%new( input, opts, 'neutral2s', i )
  enddo

  ! initialize 2D particle manager
  part_dim = 15
  ! loop over all the 2D particle profile to get the buffer size
  do i = 1, this%num_species
    call set_part2d_comm( part_dim, npmax = this%pf_spe(i)%npmax )
  enddo
  do i = 1, this%num_neutrals
    call set_part2d_comm( part_dim, npmax = this%pf_neut(i)%npmax )
  enddo
  do i = 1, this%num_neutral2s
    call set_part2d_comm( part_dim, npmax = this%pf_neut2(i)%npmax )
  enddo
  call init_part2d_comm( opts )

  do i = 1, this%num_species

    call input%get( 'species('//num2str(i)//').q', qm )
    call input%get( 'species('//num2str(i)//').m', qbm )
    qbm = qm / qbm

    push_type = p_push2_robust
    call input%get( 'species('//num2str(i)//').push_type', str )
    select case ( trim(str) )
    case ( 'robust' )
      push_type = p_push2_robust
    case ( 'clamp' )
      push_type = p_push2_clamp
    case ( 'robust-subcycling' )
      push_type = p_push2_robust_subcyc
    case ( 'explicit' )
      push_type = p_push2_explicit
    case default
      call write_err( 'Invalid pusher type! Only "robust", "clamp" and "robust-subcycling" &
        &are supported currently.' )
    end select

    call this%spe(i)%new( opts, this%pf_spe(i), ps, max_mode, &
      qbm, s, push_type, sm_type, sm_ord )

  enddo

  do i = 1, this%num_neutrals

    call input%get( 'neutrals('//num2str(i)//').q', qm )
    call input%get( 'neutrals('//num2str(i)//').m', qbm )
    call input%get( 'neutrals('//num2str(i)//').element', elem )
    call input%get( 'neutrals('//num2str(i)//').ion_max', ion_max )
    qbm = qm / qbm

    push_type = p_push2_robust
    call input%get( 'neutrals('//num2str(i)//').push_type', str )
    select case ( trim(str) )
    case ( 'robust' )
      push_type = p_push2_robust
    case ( 'clamp' )
      push_type = p_push2_clamp
    case ( 'robust-subcycling' )
      push_type = p_push2_robust_subcyc
    case default
      call write_err( 'Invalid pusher type! Only "robust", "clamp" and "robust-subcycling" &
        &are supported currently.' )
    end select

    call this%neut(i)%new( opts, this%pf_neut(i), max_mode, elem, ion_max, &
      qbm, omega_p, s, push_type, sm_type, sm_ord )

  enddo

  do i = 1, this%num_neutral2s

    call input%get( 'neutral2s('//num2str(i)//').q', qm )
    call input%get( 'neutral2s('//num2str(i)//').m', qbm )
    call input%get( 'neutral2s('//num2str(i)//').mi', qbmi)
    call input%get( 'neutral2s('//num2str(i)//').element', elem )
    call input%get( 'neutral2s('//num2str(i)//').ion_max', ion_max )
    call input%get( 'neutral2s('//num2str(i)//').v', v )
    call input%get( 'neutral2s('//num2str(i)//').sec', sec )
    qbme = qm / qbm
    qbm = -qm/qbmi

    push_type = p_push2_robust
    call input%get( 'neutral2s('//num2str(i)//').push_type', str )
    select case ( trim(str) )
    case ( 'robust' )
      push_type = p_push2_robust
    case ( 'clamp' )
      push_type = p_push2_clamp
    case ( 'robust-subcycling' )
      push_type = p_push2_robust_subcyc
    case default
      call write_err( 'Invalid pusher type! Only "robust", "clamp" and "robust-subcycling" &
        &are supported currently.' )
    end select

    call this%neut2(i)%new( opts, this%pf_neut2(i), max_mode, elem, ion_max, v, sec, &
      qbm, qbme, omega_p, s, push_type, sm_type, sm_ord )

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

  do i = 1, this%num_neutral2s
    call this%neut2(i)%del()
  enddo

  call end_part2d_comm()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_plasma

end module sim_plasma_class