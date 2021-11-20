module profile_laser_lib

use input_class
use sysutil_module
use kwargs_class
! use param

implicit none
private

public :: set_prof_perp_gaussian, get_prof_perp_gaussian
public :: set_prof_lon_sin2, get_prof_lon_sin2

contains

! ------------------------------------------------------------------------------
! GAUSSIAN TRANSVERSE PROFILES
! ------------------------------------------------------------------------------
subroutine set_prof_perp_gaussian( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  type(kw_list), intent(inout) :: prof_pars

  real :: val

  ! call input%get( trim(sect_name) // '.center', val )
  ! call prof_pars%append( 'center', val )
  call input%get( trim(sect_name) // '.w0', val )
  call prof_pars%append( 'w0', val )
  call input%get( trim(sect_name) // '.focal_distance', val )
  call prof_pars%append( 'f_dist', val )
  call input%get( trim(sect_name) // '.lon_center', val )
  call prof_pars%append( 'lon_center', val )

end subroutine set_prof_perp_gaussian

subroutine get_prof_perp_gaussian( r, z, k0, prof_pars, mode, ar_re, ar_im, ai_re, ai_im )

  implicit none
  real, intent(in) :: r, z, k0
  type(kw_list), intent(in) :: prof_pars
  integer, intent(in) :: mode
  real, intent(out) :: ar_re, ar_im, ai_re, ai_im

  real :: w0, zr, curv, f_dist, lon_center, gouy_shift, z_shift, z2, zr2, w, phase, r2, amp

  call prof_pars%get( 'w0', w0 )
  call prof_pars%get( 'f_dist', f_dist )
  call prof_pars%get( 'lon_center', lon_center )

  if ( mode == 0 ) then

    z_shift = lon_center - z - f_dist
    r2 = r * r
    z2 = z_shift * z_shift
    zr = 0.5 * k0 * w0 * w0
    zr2 = zr * zr
    curv = z_shift / ( z2 + zr2 )
    w = w0 * sqrt( 1.0 + z2 / zr2 )
    gouy_shift = atan2( z_shift, zr )
    phase = 0.5 * k0 * r2 * curv - gouy_shift
    amp = w0 / w * exp(-r2 / (w*w))

    ar_re = amp * cos(phase)
    ar_im = 0.0
    ai_re = -amp * sin(phase)
    ai_im = 0.0

  else

    ar_re = 0.0
    ar_im = 0.0
    ai_re = 0.0
    ai_im = 0.0

  endif

end subroutine get_prof_perp_gaussian

! ------------------------------------------------------------------------------
! SIN2 LONGITUDINAL PROFILES
! ------------------------------------------------------------------------------
subroutine set_prof_lon_sin2( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  type(kw_list), intent(inout) :: prof_pars

  real :: val

  call input%get( trim(sect_name) // '.lon_center', val )
  call prof_pars%append( 'lon_center', val )
  call input%get( trim(sect_name) // '.t_rise', val )
  call prof_pars%append( 't_rise', val )
  call input%get( trim(sect_name) // '.t_flat', val )
  call prof_pars%append( 't_flat', val )
  call input%get( trim(sect_name) // '.t_fall', val )
  call prof_pars%append( 't_fall', val )

end subroutine set_prof_lon_sin2

subroutine get_prof_lon_sin2( z, prof_pars, env )

  implicit none
  real, intent(in) :: z
  type(kw_list), intent(in) :: prof_pars
  real, intent(out) :: env

  real :: center, t_rise, t_fall, t_flat, flat_start, flat_end
  real, parameter :: pih = 1.570796326794897

  call prof_pars%get( 'lon_center', center )
  call prof_pars%get( 't_rise', t_rise )
  call prof_pars%get( 't_flat', t_flat )
  call prof_pars%get( 't_fall', t_fall )
  flat_start = center - 0.5 * t_flat
  flat_end   = center + 0.5 * t_flat

  if ( z < flat_start - t_rise ) then
    env = 0.0
  elseif ( z < flat_start ) then
    env = cos( (z - flat_start) / t_rise * pih )
    env = env * env
  elseif ( z < flat_end ) then
    env = 1.0
  elseif ( z < flat_end + t_fall ) then
    env = cos( (z - flat_end) / t_fall * pih )
    env = env * env
  else
    env = 0.0
  endif

end subroutine get_prof_lon_sin2

end module profile_laser_lib