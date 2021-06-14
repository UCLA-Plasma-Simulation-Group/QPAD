module fdist3d_std_lib

use input_class
use sysutil_module

implicit none

private

public :: set_prof_uniform, get_den_uniform
public :: set_prof_gaussian, get_den_gaussian

contains

! ------------------------------------------------------------------------------
! UNIFORM PROFILES
! ------------------------------------------------------------------------------

subroutine set_prof_uniform( input, sect_name, dim, prof_pars )
  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  integer, intent(in) :: dim
  real, intent(inout), dimension(:), pointer :: prof_pars
  ! placeholder for uniform profile, do nothing
  return
end subroutine set_prof_uniform

subroutine get_den_uniform( x, prof_pars, den_value )
  implicit none
  real, intent(in) :: x
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: den_value
  den_value = 1.0
end subroutine get_den_uniform

! ------------------------------------------------------------------------------
! GAUSSIAN PROFILES
! ------------------------------------------------------------------------------

subroutine set_prof_gaussian( input, sect_name, dim, prof_pars )
  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  integer, intent(in) :: dim
  real, intent(inout), dimension(:), pointer :: prof_pars
  ! local
  real :: z0

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  allocate( prof_pars(2) )

  call input%get( trim(sect_name) // '.gauss_center(' // num2str(dim) // ')', prof_pars(1) )
  call input%get( trim(sect_name) // '.gauss_sigma(' // num2str(dim) // ')', prof_pars(2) )
  if ( dim == 3 ) then
    call input%get( 'simulation.box.z(1)', z0 )
    prof_pars(1) = prof_pars(1) - z0
  endif

end subroutine set_prof_gaussian

subroutine get_den_gaussian( x, prof_pars, den_value )
  implicit none
  real, intent(in) :: x
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: den_value

  real :: sigma, mu

  mu    = prof_pars(1)  
  sigma = prof_pars(2)

  den_value = exp( -0.5 * ( (x-mu)/sigma )**2 )

end subroutine get_den_gaussian

end module fdist3d_std_lib