module sim_fields_class

use parallel_pipe_class
use grid_class
use field_psi_class
use field_e_class
use field_b_class
use field_vpot_class
use field_src_class
use sysutil
use param
use input_class

implicit none

private

public :: sim_fields

type sim_fields

  ! private

  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()

  type( field_psi ), pointer :: psi => null()
  type( field_vpot ), pointer :: vpot => null() ! just for diagnostic
  type( field_b ), pointer :: b_spe => null(), b_beam => null(), b => null()
  type( field_e ), pointer :: e_spe => null(), e_beam => null(), e => null()
  type( field_jay ), pointer :: cu => null(), amu => null()
  type( field_rho ), pointer :: q_spe => null(), q_beam => null()
  type( field_djdxi ), pointer :: dcu => null(), acu => null()

  contains

  generic :: new => init_sim_fields
  generic :: del => end_sim_fields

  procedure, private :: init_sim_fields, end_sim_fields

end type sim_fields

character(len=18), save :: cls_name = 'sim_fields'
integer, save :: cls_level = 2

contains

subroutine init_sim_fields( this, input )

  implicit none

  class( sim_fields ), intent(inout) :: this
  type( input_json ), pointer, intent(inout) :: input

  ! local data
  character(len=18), save :: sname = 'init_sim_fields'
  character(len=:), allocatable :: str
  integer :: entity, max_mode, ps, bnd, sm_type, sm_ord

  this%gp => input%gp
  this%pp => input%pp

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  allocate( this%psi, &
    this%e_spe, this%e_beam, this%e, &
    this%b_spe, this%b_beam, this%b, &
    this%cu, this%acu, this%amu, this%dcu, &
    this%q_spe, this%q_beam )

  call input%get( 'simulation.max_mode', max_mode )

  ! read field boundary type
  call input%get( 'simulation.field_boundary', str )
  select case ( trim(str) )
  case ( 'open' )
    bnd = p_bnd_open
  case ( 'zero' )
    bnd = p_bnd_zero
  case default
    call write_err( 'Invalid field boundary type! Only "open" and "zero" are supported currently.' )
  end select

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

  call this%psi%new(    this%pp, this%gp, max_mode, ps, bnd )
  call this%q_spe%new(  this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%q_beam%new( this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%cu%new(     this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%dcu%new(    this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%acu%new(    this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%amu%new(    this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  entity = p_entity_plasma
  call this%e_spe%new(  this%pp, this%gp, max_mode, ps, bnd, entity )
  call this%b_spe%new(  this%pp, this%gp, max_mode, ps, bnd, entity )
  call this%e%new(      this%pp, this%gp, max_mode, ps, bnd, entity )
  call this%b%new(      this%pp, this%gp, max_mode, ps, bnd, entity )
  entity = p_entity_beam
  call this%e_beam%new( this%pp, this%gp, max_mode, ps, bnd, entity )
  call this%b_beam%new( this%pp, this%gp, max_mode, ps, bnd, entity )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_fields

subroutine end_sim_fields( this )

  implicit none

  class(sim_fields), intent(inout) :: this
  ! local data
  character(len=18), save :: sname = 'end_sim_fields'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%psi%del()
  call this%e_spe%del()
  call this%b_spe%del()
  call this%e_beam%del()
  call this%b_beam%del()
  call this%cu%del()
  call this%dcu%del()
  call this%acu%del()
  call this%amu%del()
  call this%q_spe%del()
  call this%q_beam%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_fields

end module sim_fields_class