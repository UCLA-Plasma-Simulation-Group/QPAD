module profile_laser_lib

use input_class
use sysutil_module
use kwargs_class
use math_module

implicit none
private

public :: set_prof_perp_gaussian, get_prof_perp_gaussian
public :: set_prof_perp_laguerre, get_prof_perp_laguerre
public :: set_prof_lon_sin2, get_prof_lon_sin2
public :: set_prof_lon_poly, get_prof_lon_poly

integer, parameter :: laguerre_p_max = 5
integer, parameter :: laguerre_l_max = 5
real, dimension(0:laguerre_p_max, 0:laguerre_l_max), parameter :: laguerre_norm_fac = &
  reshape( (/ &
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
    .428881942195642, .587209029015734, .715480974477903, .824715297721105, .921317171436064, 1.00881288923537, &
    .367879441171442, .606530658706018, .849813805803178, 1.09328889163892, 1.33673290295862, 1.58013471092357, &
    .409916278920296, .771854776322994, 1.20512863331838, 1.69577853394018, 2.23759926979643, 2.82618005138511, &
    .541341132211037, 1.13063145394112, 1.92864111604763, 2.92737142922967, 4.12615415717722, 5.52484537726148, &
    .811173616512526, 1.84462325501195, 3.39116963304902, 5.49617869099731, 8.21093909438522, 11.5834999136547 &
  /), (/6, 6/) )

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
! LAGUERRE-GAUSSIAN TRANSVERSE PROFILES
! ------------------------------------------------------------------------------
subroutine set_prof_perp_laguerre( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  type(kw_list), intent(inout) :: prof_pars

  real :: rval
  integer :: ival

  call input%get( trim(sect_name) // '.w0', rval )
  call prof_pars%append( 'w0', rval )
  call input%get( trim(sect_name) // '.focal_distance', rval )
  call prof_pars%append( 'f_dist', rval )
  call input%get( trim(sect_name) // '.lon_center', rval )
  call prof_pars%append( 'lon_center', rval )
  call input%get( trim(sect_name) // '.radial_index', ival )
  call prof_pars%append( 'radial_index', ival )
  call input%get( trim(sect_name) // '.phi_index', ival )
  call prof_pars%append( 'phi_index', ival )

end subroutine set_prof_perp_laguerre

subroutine get_prof_perp_laguerre( r, z, k0, prof_pars, mode, ar_re, ar_im, ai_re, ai_im )

  implicit none
  real, intent(in) :: r, z, k0
  type(kw_list), intent(in) :: prof_pars
  integer, intent(in) :: mode
  real, intent(out) :: ar_re, ar_im, ai_re, ai_im

  real :: w0, zr, curv, f_dist, lon_center, gouy_shift, z_shift, z2, zr2, phase, w2, r2, r2_iw2, amp
  integer :: p, l
  real, parameter :: sqrt2 = 1.414213562373095

  call prof_pars%get( 'w0', w0 )
  call prof_pars%get( 'f_dist', f_dist )
  call prof_pars%get( 'lon_center', lon_center )
  call prof_pars%get( 'radial_index', p )
  call prof_pars%get( 'phi_index', l )

  if ( mode == l ) then

    z_shift = lon_center - z - f_dist
    r2 = r * r
    z2 = z_shift * z_shift
    zr = 0.5 * k0 * w0 * w0
    zr2 = zr * zr
    curv = z_shift / ( z2 + zr2 )
    w2 = w0 * w0 * ( 1.0 + z2 / zr2 )
    r2_iw2 = r2 / w2
    gouy_shift = real( 1 + 2 * p + abs(l) ) * atan2( z_shift, zr )
    phase = 0.5 * k0 * r2 * curv - gouy_shift

    ! This is the definition from Wikipedia
    ! amp = w0 / sqrt(w2) * exp(-r2_iw2) * ( sqrt2 * sqrt(r2_iw2) ) ** abs(l) &
    !       * sqrt( 2.0 * factorial(p) / ( pi * factorial(p+abs(l)) ) ) &
    !       * laguerre( p, abs(l), 2.0 * r2_iw2 )

    ! This is the definition from OSIRIS
    amp = w0 / sqrt(w2) * exp(-r2_iw2) * sqrt(r2_iw2) ** abs(l) &
          * laguerre( p, abs(l), 2.0 * r2_iw2 )
          ! * laguerre( p, abs(l), 2.0 * r2_iw2 )
    if ( p <= laguerre_p_max .and. l <= laguerre_l_max ) then
      amp = amp * 1.0
      ! amp = amp / laguerre_norm_fac( p, abs(l) )
    else
      call write_stdout( '[WARNING] The laser peak intensity is not normalized! &
        &The maximum mode supporting normalization is (p,l)=(5,5).' )
    endif

    if ( mode == 0 ) then
      ar_re = amp * cos(phase)
      ar_im = amp * sin(phase)
      ai_re = 0.0
      ai_im = 0.0
    else
      ar_re =  0.5 * amp * cos(phase)
      ar_im =  0.5 * amp * sin(phase)
      ai_re = -ar_im
      ai_im =  ar_re
    endif

  else

    ar_re = 0.0
    ar_im = 0.0
    ai_re = 0.0
    ai_im = 0.0

  endif

end subroutine get_prof_perp_laguerre

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

! ------------------------------------------------------------------------------
! POLYNOMIAL LONGITUDINAL PROFILES
! ------------------------------------------------------------------------------
subroutine set_prof_lon_poly( input, sect_name, prof_pars )

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

end subroutine set_prof_lon_poly

subroutine get_prof_lon_poly( z, prof_pars, env )

  implicit none
  real, intent(in) :: z
  type(kw_list), intent(in) :: prof_pars
  real, intent(out) :: env

  real :: center, t_rise, t_fall, t_flat, flat_start, flat_end, t

  call prof_pars%get( 'lon_center', center )
  call prof_pars%get( 't_rise', t_rise )
  call prof_pars%get( 't_flat', t_flat )
  call prof_pars%get( 't_fall', t_fall )
  flat_start = center - 0.5 * t_flat
  flat_end   = center + 0.5 * t_flat

  if ( z < flat_start - t_rise ) then
    env = 0.0
  elseif ( z < flat_start ) then
    t = ( z + t_rise - flat_start ) / t_rise
    env = 10.0 * t**3 - 15.0 * t**4 + 6.0 * t**5
  elseif ( z < flat_end ) then
    env = 1.0
  elseif ( z < flat_end + t_fall ) then
    t = ( flat_end + t_fall - z ) / t_fall
    env = 10.0 * t**3 - 15.0 * t**4 + 6.0 * t**5
  else
    env = 0.0
  endif

end subroutine get_prof_lon_poly

end module profile_laser_lib