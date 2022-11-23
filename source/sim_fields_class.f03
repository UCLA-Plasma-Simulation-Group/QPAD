module sim_fields_class

use parallel_module
use options_class
use field_psi_class
use field_e_class
use field_b_class
use field_vpot_class
use field_src_class
use sysutil_module
use param
use input_class

implicit none

private

public :: sim_fields

type sim_fields

  ! private

  class( field_psi ), pointer :: psi => null()
  class( field_vpot ), pointer :: vpot => null() ! just for diagnostic
  class( field_b ), pointer :: b_spe => null(), b_beam => null(), b => null()
  class( field_e ), pointer :: e_spe => null(), e_beam => null(), e => null()
  class( field_jay ), pointer :: cu => null(), amu => null(), gamma => null()
  class( field_rho ), pointer :: q_spe => null(), q_beam => null()
  !class( field_djdxi ), pointer :: dcu => null(), acu => null()

  contains

  procedure :: alloc => alloc_sim_fields
  procedure :: new   => init_sim_fields
  procedure :: del   => end_sim_fields

end type sim_fields

character(len=18), save :: cls_name = 'sim_fields'
integer, save :: cls_level = 2

contains

subroutine alloc_sim_fields( this, input )

  implicit none

  class( sim_fields ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  ! local data
  integer :: max_mode
  character(len=32), save :: sname = 'alloc_sim_fields'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%psi ) )    allocate( field_psi :: this%psi )
  if ( .not. associated( this%b_spe ) )  allocate( field_b :: this%b_spe )
  if ( .not. associated( this%b_beam ) ) allocate( field_b :: this%b_beam )
  if ( .not. associated( this%b ) )      allocate( field_b :: this%b )
  if ( .not. associated( this%e_spe ) )  allocate( field_e :: this%e_spe )
  if ( .not. associated( this%e_beam ) ) allocate( field_e :: this%e_beam )
  if ( .not. associated( this%e ) )      allocate( field_e :: this%e )
  if ( .not. associated( this%cu ) )     allocate( field_jay :: this%cu )
  if ( .not. associated( this%amu ) )    allocate( field_jay :: this%amu )
  if ( .not. associated( this%gamma ) )  allocate( field_rho :: this%gamma)
  if ( .not. associated( this%q_spe ) )  allocate( field_rho :: this%q_spe )
  if ( .not. associated( this%q_beam ) ) allocate( field_rho :: this%q_beam )
!   if ( .not. associated( this%dcu ) )    allocate( field_djdxi :: this%dcu )
!   if ( .not. associated( this%acu ) )    allocate( field_djdxi :: this%acu )

  call input%get( 'simulation.max_mode', max_mode )

  call this%psi%alloc( max_mode )
  call this%b_spe%alloc( max_mode )
  call this%b_beam%alloc( max_mode )
  call this%b%alloc( max_mode )
  call this%e_spe%alloc( max_mode )
  call this%e_beam%alloc( max_mode )
  call this%e%alloc( max_mode )
  call this%cu%alloc( max_mode )
  call this%amu%alloc( max_mode )
  call this%gamma%alloc( max_mode )
  call this%q_spe%alloc( max_mode )
  call this%q_beam%alloc( max_mode )
!   call this%dcu%alloc( max_mode )
!   call this%acu%alloc( max_mode )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_sim_fields

subroutine init_sim_fields( this, input, opts )

  implicit none

  class( sim_fields ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts

  ! local data
  character(len=18), save :: sname = 'init_sim_fields'
  character(len=:), allocatable :: str
  integer :: entity, max_mode, ps, bnd, sm_type, sm_ord

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ! allocate( this%psi, &
  !   this%e_spe, this%e_beam, this%e, &
  !   this%b_spe, this%b_beam, this%b, &
  !   this%cu, this%acu, this%amu, this%dcu, &
  !   this%q_spe, this%q_beam )

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

  call this%psi%new(    opts, max_mode, ps, bnd )
  call this%q_spe%new(  opts, max_mode, ps, sm_type, sm_ord )
  call this%q_beam%new( opts, max_mode, ps, sm_type, sm_ord )
  call this%cu%new(     opts, max_mode, ps, sm_type, sm_ord )
!   call this%dcu%new(    opts, max_mode, ps, sm_type, sm_ord )
!   call this%acu%new(    opts, max_mode, ps, sm_type, sm_ord )
  call this%amu%new(    opts, max_mode, ps, sm_type, sm_ord )
  call this%gamma%new(    opts, max_mode, ps, sm_type, sm_ord )
  entity = p_entity_plasma
  call this%e_spe%new(  opts, max_mode, ps, bnd, entity )
  call this%b_spe%new(  opts, max_mode, ps, bnd, entity )
  call this%e%new(      opts, max_mode, ps, bnd, entity )
  call this%b%new(      opts, max_mode, ps, bnd, entity )
  entity = p_entity_beam
  call this%e_beam%new( opts, max_mode, ps, bnd, entity )
  call this%b_beam%new( opts, max_mode, ps, bnd, entity )

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
!   call this%dcu%del()
!   call this%acu%del()
  call this%amu%del()
  call this%gamma%del()
  call this%q_spe%del()
  call this%q_beam%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_fields

end module sim_fields_class