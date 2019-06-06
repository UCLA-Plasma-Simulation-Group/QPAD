module sim_fields_class

use parallel_pipe_class
use grid_class
use field_psi_class
use field_e_class
use field_b_class
use field_src_class
use sys
use param
use input_class

implicit none

private

public :: sim_fields

type sim_fields

  ! private

  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()

  type( field_psi ), allocatable :: psi
  type( field_b ), allocatable :: b_spe, b_beam, b
  type( field_e ), allocatable :: e_spe, e_beam, e
  type( field_jay ), allocatable :: cu, amu
  type( field_rho ), allocatable :: q_spe, q_beam
  type( field_djdxi ), allocatable :: dcu, acu

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
  character(len=20) :: s1, s2, s3
  character(len=:), allocatable :: str
  integer :: i,n,ndump,j,k,l,m
  integer :: entity, max_mode, ps, bnd, sm_type, sm_ord
  real :: dr, dxi, tol, min, max

  this%gp => input%gp
  this%pp => input%pp

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  allocate( this%psi, &
    this%e_spe, this%e_beam, this%e, &
    this%b_spe, this%b_beam, this%b, &
    this%cu, this%acu, this%amu, this%dcu, &
    this%q_spe, this%q_beam )

  call input%get( 'simulation.max_mode', max_mode )
  call input%get( 'simulation.solver_tol', tol )
  
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

  ! call input%get( 'simulation.box.r(1)', min )
  ! call input%get( 'simulation.box.r(2)', max )
  ! dr = ( max - min ) / this%gp%get_nd(1)
  ! call input%get( 'simulation.box.z(1)', min )
  ! call input%get( 'simulation.box.z(2)', max )
  ! dxi = ( max - min ) / this%gp%get_nd(2)

  call this%psi%new(    this%pp, this%gp, max_mode, ps, bnd, tol )
  call this%q_spe%new(  this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%q_beam%new( this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%cu%new(     this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%dcu%new(    this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%acu%new(    this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  call this%amu%new(    this%pp, this%gp, max_mode, ps, sm_type, sm_ord )
  entity = p_entity_plasma
  call this%e_spe%new(  this%pp, this%gp, max_mode, ps, bnd, entity, tol )
  call this%b_spe%new(  this%pp, this%gp, max_mode, ps, bnd, entity, tol )
  call this%e%new(      this%pp, this%gp, max_mode, ps, bnd, entity, tol )
  call this%b%new(      this%pp, this%gp, max_mode, ps, bnd, entity, tol )
  entity = p_entity_beam
  call this%e_beam%new( this%pp, this%gp, max_mode, ps, bnd, entity, tol )
  call this%b_beam%new( this%pp, this%gp, max_mode, ps, bnd, entity, tol )

  ! call input%get('simulation.nspecies',n)

  ! loop1: do i = 1, n
  !   write (s1, '(I4.4)') i
  !   call input%info('species('//trim(s1)//').diag',n_children=m)
  !   do j = 1, m
  !      write (s2, '(I4.4)') j
  !      call input%get('species('//trim(s1)//').diag'//'('//trim(s2)//').ndump',ndump)
  !      if (ndump>0) then
  !         call input%info('species('//trim(s1)//').diag'//'('//trim(s2)//').name',n_children=l)
  !         do k = 1, l
  !            write (s3, '(I4.4)') k
  !            if(allocated(ff)) deallocate(ff)
  !            call input%get('species('//trim(s1)//').diag'//'('//trim(s2)//').name'&
  !            &//'('//trim(s3)//')',ff)
  !            if (ff == 'jx' .or. ff == 'jy' .or. ff == 'jz') then
  !               allocate(this%cu3d)
  !               call this%cu3d%new(this%p,this%err,this%sp3,dim=1)
  !               exit loop1
  !            end if
  !         end do
  !      end if
  !   end do
  ! end do loop1

  ! call input%info('field.diag',n_children=n)

  ! loop2: do i = 1, n
  !   write (s1,'(I4.4)') i
  !   call input%get('field.diag('//trim(s1)//').ndump',ndump)
  !   if (ndump > 0) then
  !      call input%info('field.diag('//trim(s1)//').name',n_children=m)
  !      do j = 1, m
  !         write (s2,'(I4.4)') j
  !         if(allocated(ff)) deallocate(ff)
  !         call input%get('field.diag('//trim(s1)//').name('//trim(s2)//')',ff)
  !         if (ff == 'psi') then
  !            allocate(this%psi3d)
  !            call this%psi3d%new(this%p,this%err,this%sp3,dim=1)
  !            exit loop2
  !         end if
  !      end do
  !   end if
  ! end do loop2

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_fields

subroutine end_sim_fields( this )

  implicit none

  class(sim_fields), intent(inout) :: this
  ! local data
  character(len=18), save :: sname = 'end_sim_fields'
  integer :: i, n

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