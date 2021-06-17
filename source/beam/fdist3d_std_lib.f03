module fdist3d_std_lib

use input_class
use sysutil_module

implicit none

private

public :: set_prof_uniform, get_den_uniform
public :: set_prof_gaussian, get_den_gaussian
public :: set_prof_pw_linear, get_den_pw_linear

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

! ------------------------------------------------------------------------------
! PIECEWISE LINEAR PROFILES
! ------------------------------------------------------------------------------

subroutine set_prof_pw_linear( input, sect_name, dim, prof_pars )
  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  integer, intent(in) :: dim
  real, intent(inout), dimension(:), pointer :: prof_pars
  ! local
  integer :: len_x, i
  real, dimension(:), allocatable :: fx, x
  real :: z0

  call input%get( trim(sect_name) // '.piecewise_fx' // num2str(dim), fx )
  call input%get( trim(sect_name) // '.piecewise_x' // num2str(dim), x )

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  len_x = size(x)
  allocate( prof_pars( 2 * len_x ) )

  do i = 2, len_x
    if ( x(i) <= x(i-1) ) then
      call write_err( 'Position array must be monotonically increasing!' )
    endif
  enddo

  do i = 1, len_x
    prof_pars(i) = x(i)
    prof_pars(len_x+i) = fx(i)
  enddo
  if ( dim == 3 ) then
    call input%get( 'simulation.box.z(1)', z0 )
    prof_pars(1:len_x) = prof_pars(1:len_x) - z0
  endif

  deallocate( fx, x )

end subroutine set_prof_pw_linear

subroutine get_den_pw_linear( x, prof_pars, den_value )
  implicit none
  real, intent(in) :: x
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: den_value

  integer :: len_x, i
  real, dimension(:), pointer :: x_array => null(), fx_array => null()

  len_x = size( prof_pars ) / 2

  ! use pointer to associate with the segment of data to avoid hard-copy
  x_array  => prof_pars( 1 : len_x )
  fx_array => prof_pars( len_x+1 : 2*len_x )

  ! check boundaries
  if ( x < x_array(1) .or. x > x_array(len_x) ) then
    call write_err( 'Density profile at current position is undefined!' )
  endif

  do i = 2, len_x
    if ( x <= x_array(i) ) then
      den_value = fx_array(i-1) + ( fx_array(i) - fx_array(i-1) ) / &
        ( x_array(i) - x_array(i-1) ) * ( x - x_array(i-1) )
      exit
    endif
  enddo

end subroutine get_den_pw_linear

end module fdist3d_std_lib