module fdist3d_rnd_lib

use input_class
use sysutil_module
use math_module

implicit none

private

public :: set_prof_uniform, get_rndpos_uniform
public :: set_prof_gaussian, get_rndpos_gaussian
public :: set_prof_parabolic, get_rndpos_parabolic
public :: set_prof_pw_linear, get_rndpos_pw_linear

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
  ! local
  real :: z0

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  allocate( prof_pars(2) )

  call input%get( trim(sect_name) // '.range' // num2str(dim) // '(1)', prof_pars(1) )
  call input%get( trim(sect_name) // '.range' // num2str(dim) // '(2)', prof_pars(2) )
  if ( dim == 3 ) then
    call input%get( 'simulation.box.z(1)', z0 )
    prof_pars(1) = prof_pars(1) - z0
    prof_pars(2) = prof_pars(2) - z0
  endif

end subroutine set_prof_uniform

subroutine get_rndpos_uniform( prof_pars, pos )
  implicit none
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: pos
  
  real :: tmp, range_l, range_u

  range_l = prof_pars(1)
  range_u = prof_pars(2)
  call random_number(tmp)
  pos = tmp * ( range_u - range_l ) + range_l

end subroutine get_rndpos_uniform

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

subroutine get_rndpos_gaussian( prof_pars, pos )
  implicit none
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: pos

  real :: sigma, mu

  mu    = prof_pars(1)  
  sigma = prof_pars(2)
  pos   = ranorm() * sigma + mu

end subroutine get_rndpos_gaussian

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

subroutine get_rndpos_parabolic( prof_pars, pos )
  implicit none
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: pos

  real :: r0, x0, x, f, dice

  x0 = prof_pars(1)  
  r0 = prof_pars(2)
  call random_number( dice )

  x = 2.0 * cos( (pi + acos(2.0*dice-1.0)) / 3.0 )
  pos = x0 + x * r0

end subroutine get_rndpos_parabolic

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
  real, dimension(:), allocatable :: pdf, cdf, x
  real :: z0

  call input%get( trim(sect_name) // '.piecewise_fx' // num2str(dim), pdf )
  call input%get( trim(sect_name) // '.piecewise_x' // num2str(dim), x )

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  len_x = size(x)
  allocate( prof_pars( 3 * len_x ), cdf(len_x) )

  ! check data validity
  do i = 2, len_x
    if ( x(i) <= x(i-1) ) then
      call write_err( 'Position array must be monotonically increasing!' )
    endif
  enddo

  ! integrate the PDF to obtain the CDF
  cdf(1) = 0.0
  do i = 2, len_x
    cdf(i) = cdf(i-1) + 0.5 * ( pdf(i) + pdf(i-1) ) * ( x(i) - x(i-1) )
  enddo

  ! normalize
  pdf = pdf / cdf(len_x)
  cdf = cdf / cdf(len_x)

  do i = 1, len_x
    prof_pars(i)         = x(i)
    prof_pars(len_x+i)   = pdf(i)
    prof_pars(2*len_x+i) = cdf(i)
  enddo

  if ( dim == 3 ) then
    call input%get( 'simulation.box.z(1)', z0 )
    prof_pars(1:len_x) = prof_pars(1:len_x) - z0
  endif

  deallocate( pdf, cdf, x )

end subroutine set_prof_pw_linear

subroutine get_rndpos_pw_linear( prof_pars, pos )

  implicit none
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: pos

  integer :: len_x, i
  real, dimension(:), pointer :: x => null(), pdf => null(), cdf => null()
  real :: dice, a, b, c

  len_x = size( prof_pars ) / 3

  ! use pointer to associate with the segment of data to avoid hard-copy
  x   => prof_pars( 1 : len_x )
  pdf => prof_pars( len_x+1 : 2*len_x )
  cdf => prof_pars( 2*len_x+1 : 3*len_x )

  call random_number( dice )

  do i = 2, len_x
    if ( cdf(i) > dice ) then
      a = 0.5 * ( pdf(i) - pdf(i-1) ) / ( x(i) - x(i-1) )
      b = pdf(i-1)
      c = cdf(i-1) - dice
      pos = 0.5 * ( sqrt( b*b - 4.0*a*c ) - b ) / a + x(i-1)
      exit
    endif
  enddo

end subroutine get_rndpos_pw_linear

end module fdist3d_rnd_lib