module fdist3d_std_lib

use input_class
use sysutil_module

implicit none

private

public :: set_prof_uniform, get_den_uniform
public :: set_prof_gaussian, get_den_gaussian
public :: set_prof_parabolic, get_den_parabolic
public :: set_prof_rational, get_den_rational
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
! PARABOLIC PROFILES
! ------------------------------------------------------------------------------

subroutine set_prof_parabolic( input, sect_name, dim, prof_pars )
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

  call input%get( trim(sect_name) // '.parabolic_center(' // num2str(dim) // ')', prof_pars(1) )
  call input%get( trim(sect_name) // '.parabolic_radius(' // num2str(dim) // ')', prof_pars(2) )
  if ( dim == 3 ) then
    call input%get( 'simulation.box.z(1)', z0 )
    prof_pars(1) = prof_pars(1) - z0
  endif

end subroutine set_prof_parabolic

subroutine get_den_parabolic( x, prof_pars, den_value )
  implicit none
  real, intent(in) :: x
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: den_value

  real :: x0, r0, tmp

  x0 = prof_pars(1)  
  r0 = prof_pars(2)

  den_value = 0.0
  tmp = (x-x0)**2 / (r0*r0)
  if ( tmp < 1.0 ) den_value = 1.0 - tmp

end subroutine get_den_parabolic

! ------------------------------------------------------------------------------
! RATIONAL PROFILES
! ------------------------------------------------------------------------------

subroutine set_prof_rational( input, sect_name, dim, prof_pars )
  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  integer, intent(in) :: dim
  real, intent(inout), dimension(:), pointer :: prof_pars
  ! local
  integer :: len_pn, len_pd, i
  real, dimension(:), allocatable :: pn, pd
  real :: z0, p_x0

  call input%get( trim(sect_name) // '.poly_numerator' // num2str(dim), pn )
  call input%get( trim(sect_name) // '.poly_denominator' // num2str(dim), pd )
  call input%get( trim(sect_name) // '.poly_x0(' // num2str(dim) // ')', p_x0 )

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  len_pn = size( pn )
  len_pd = size( pd )
  allocate( prof_pars( len_pn + len_pd + 3 ) )

  prof_pars(1) = p_x0
  prof_pars(2) = real(len_pn)
  prof_pars(3) = real(len_pd)
  prof_pars(4:len_pn+3) = pn
  prof_pars(len_pn+4:len_pn+len_pd+3) = pd
  if ( dim == 3 ) then
    call input%get( 'simulation.box.z(1)', z0 )
    prof_pars(1) = prof_pars(1) - z0
  endif

  deallocate( pn, pd )

end subroutine set_prof_rational

subroutine get_den_rational( x, prof_pars, den_value )
  implicit none
  real, intent(in) :: x
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: den_value

  integer :: len_pn, len_pd, i
  real, dimension(:), pointer :: pn => null(), pd => null()
  real :: p_x0, pn_val, pd_val

  p_x0 = prof_pars(1)
  len_pn = int( prof_pars(2) )
  len_pd = int( prof_pars(3) )
  pn => prof_pars(4:len_pn+3)
  pd => prof_pars(len_pn+4:len_pn+len_pd+3)

  pn_val = 0.0
  do i = 1, len_pn
    pn_val = pn_val + pn(i) * ( x - p_x0 )**(i-1)
  enddo
  pd_val = 0.0
  do i = 1, len_pd
    pd_val = pd_val + pd(i) * ( x - p_x0 )**(i-1)
  enddo

  den_value = pn_val / pd_val

end subroutine get_den_rational

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
    call write_wrn( 'Density profile at current position is undefined! Set zero &
      &charge for this particle.' )
    den_value = 0.0
    return
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