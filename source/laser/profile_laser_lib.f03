module profile_laser_lib

use input_class
use sysutil_module
use kwargs_class
use math_module
use m_fparser
implicit none
private

public :: set_prof_perp_gaussian, get_prof_perp_gaussian
public :: set_prof_perp_laguerre, get_prof_perp_laguerre
public :: set_prof_lon_sin2, get_prof_lon_sin2
public :: set_prof_lon_poly, get_prof_lon_poly
public :: set_prof_lon_pw_linear, get_prof_lon_pw_linear
public :: set_prof_lon_cubic_spline, get_prof_lon_cubic_spline
public :: set_prof_perp_astrl_analytic, get_prof_perp_astrl_analytic
public :: set_prof_perp_astrl_discrete, get_prof_perp_astrl_discrete, normalize_a0_astrl_discrete
public :: set_prof_const, get_prof_const

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

end subroutine set_prof_perp_gaussian

subroutine get_prof_perp_gaussian( r, z, t, k, k0, prof_pars, mode, ar_re, ar_im, ai_re, ai_im )

  implicit none
  real, intent(in) :: r, z, t, k, k0
  type(kw_list), intent(in) :: prof_pars
  integer, intent(in) :: mode
  real, intent(out) :: ar_re, ar_im, ai_re, ai_im

  real :: w0, zr, curv, f_dist, gouy_shift, z_shift, z2, zr2, w, phase, r2, amp

  call prof_pars%get( 'w0', w0 )
  call prof_pars%get( 'f_dist', f_dist )

  if ( mode == 0 ) then

    z_shift = -1.0 * (z + f_dist)
    r2 = r * r
    z2 = z_shift * z_shift
    zr = 0.5 * k * w0 * w0
    zr2 = zr * zr
    curv = z_shift / ( z2 + zr2 )
    w = w0 * sqrt( 1.0 + z2 / zr2 )
    gouy_shift = atan2( z_shift, zr )
    phase = 0.5 * k * r2 * curv - gouy_shift - (k - k0) * t
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
  call input%get( trim(sect_name) // '.radial_index', ival )
  call prof_pars%append( 'radial_index', ival )
  call input%get( trim(sect_name) // '.phi_index', ival )
  call prof_pars%append( 'phi_index', ival )

end subroutine set_prof_perp_laguerre

subroutine get_prof_perp_laguerre( r, z, t, k, k0, prof_pars, mode, ar_re, ar_im, ai_re, ai_im )

  implicit none
  real, intent(in) :: r, z, t, k, k0
  type(kw_list), intent(in) :: prof_pars
  integer, intent(in) :: mode
  real, intent(out) :: ar_re, ar_im, ai_re, ai_im

  real :: w0, zr, curv, f_dist, gouy_shift, z_shift, z2, zr2, phase, w2, r2, r2_iw2, amp
  integer :: p, l
  real, parameter :: sqrt2 = 1.414213562373095

  call prof_pars%get( 'w0', w0 )
  call prof_pars%get( 'f_dist', f_dist )
  call prof_pars%get( 'radial_index', p )
  call prof_pars%get( 'phi_index', l )

  if ( mode == l ) then

    z_shift = -1.0 * ( z + f_dist )
    r2 = r * r
    z2 = z_shift * z_shift
    zr = 0.5 * k * w0 * w0
    zr2 = zr * zr
    curv = z_shift / ( z2 + zr2 )
    w2 = w0 * w0 * ( 1.0 + z2 / zr2 )
    r2_iw2 = r2 / w2
    gouy_shift = real( 1 + 2 * p + abs(l) ) * atan2( z_shift, zr )
    phase = 0.5 * k * r2 * curv - gouy_shift - (k - k0) * t

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
      ar_im = 0.0
      ai_re = -amp * sin(phase)
      ai_im = 0.0
    else
      ! Notice the different order here for convenience
      ar_re =  0.5 * amp * cos(phase)
      ai_re = -0.5 * amp * sin(phase)
      ar_im = -ai_re
      ai_im = ar_re
    endif

  else

    ar_re = 0.0
    ar_im = 0.0
    ai_re = 0.0
    ai_im = 0.0

  endif

end subroutine get_prof_perp_laguerre

! ------------------------------------------------------------------------------
! ASTRL ANALYTIC TRANSVERSE PROFILES
! ------------------------------------------------------------------------------
subroutine set_prof_perp_astrl_analytic( input, sect_name, prof_pars, math_funcs )

  implicit none

  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  type(kw_list), intent(inout) :: prof_pars
  character(len=:), allocatable :: read_str
  type(t_fparser), dimension(:), pointer, intent(inout) :: math_funcs
  real :: val
  real(p_k_fparse) , dimension(1) :: fparser_arr
  integer :: ierr
  real :: rval
  integer :: ival

  if ( .not. associated( math_funcs ) ) then
    allocate( t_fparser :: math_funcs( 3 ) )
  endif

  call input%get( trim(sect_name) // '.w0_math_func', read_str )
  call setup(math_funcs(1), trim(read_str), (/'xi'/), ierr)

  call input%get( trim(sect_name) // '.a0_math_func', read_str )
  call setup(math_funcs(2), trim(read_str), (/'xi'/), ierr)

  call input%get( trim(sect_name) // '.s0_math_func', read_str )
  call setup(math_funcs(3), trim(read_str), (/'xi'/), ierr)

end subroutine set_prof_perp_astrl_analytic

subroutine get_prof_perp_astrl_analytic( r, z, t, k, k0, prof_pars, math_funcs, mode, ar_re, ar_im, ai_re, ai_im )

  implicit none
  real, intent(in) :: r, z, t, k, k0
  type(kw_list), intent(in) :: prof_pars
  integer, intent(in) :: mode
  real, intent(out) :: ar_re, ar_im, ai_re, ai_im
  character(len=:), allocatable :: read_str
  real :: a0,w0, zr, curv, f_dist, gouy_shift, z_shift, z2, zr2, w, phase, r2, amp
  real(p_k_fparse), dimension(1) :: fparser_arr 
  type(t_fparser), dimension(:), pointer, intent(inout) :: math_funcs

  fparser_arr(1) = z
  w0 = eval(math_funcs(1), fparser_arr)
  a0 = eval(math_funcs(2), fparser_arr)
  f_dist = eval(math_funcs(3), fparser_arr)

  if ( mode == 0 ) then

    z_shift = -1.0 * (z + f_dist)
    r2 = r * r
    z2 = z_shift * z_shift
    zr = 0.5 * k * w0 * w0
    zr2 = zr * zr
    curv = z_shift / ( z2 + zr2 )
    w = w0 * sqrt( 1.0 + z2 / zr2 )
    gouy_shift = atan2( z_shift, zr )
    phase = 0.5 * k * r2 * curv - gouy_shift - (k - k0) * z
    amp = w0 / w * exp(-r2 / (w*w))

    ar_re = a0 * amp * cos(phase)
    ar_im = 0.0
    ai_re = -a0 * amp * sin(phase)
    ai_im = 0.0

  else

    ar_re = 0.0
    ar_im = 0.0
    ai_re = 0.0
    ai_im = 0.0

  endif

end subroutine get_prof_perp_astrl_analytic


! ------------------------------------------------------------------------------
! ASTRL DISCRETE TRANSVERSE PROFILES
! ------------------------------------------------------------------------------
subroutine set_prof_perp_astrl_discrete( input, sect_name, prof_pars, math_funcs )

  implicit none

  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  type(kw_list), intent(inout) :: prof_pars
  character(len=:), allocatable :: read_str
  type(t_fparser), dimension(:), pointer, intent(inout) :: math_funcs
  real :: val
  real(p_k_fparse) , dimension(1) :: fparser_arr
  integer :: ierr
  integer :: ival

  if ( .not. associated( math_funcs ) ) then
    allocate( t_fparser :: math_funcs( 4 ) )
  endif

  call input%get( trim(sect_name) // '.w0_math_func', read_str )
  call setup(math_funcs(1), trim(read_str), (/'xi'/), ierr)

  call input%get( trim(sect_name) // '.a0_math_func', read_str )
  call setup(math_funcs(2), trim(read_str), (/'xi'/), ierr)

  call input%get( trim(sect_name) // '.s0_math_func', read_str )
  call setup(math_funcs(3), trim(read_str), (/'xi'/), ierr)

  call input%get( trim(sect_name) // '.pulselet_math_func', read_str )
  call setup(math_funcs(4), trim(read_str), (/'xi'/), ierr)

  call input%get( trim(sect_name) // '.pulselet_delay', val )
  call prof_pars%append( 'pulselet_delay', val )

  call input%get( trim(sect_name) // '.pulselet_range', val )
  call prof_pars%append( 'pulselet_range', val )

  call input%get( trim(sect_name) // '.pulselet_offset', val )
  call prof_pars%append( 'pulselet_offset', val )

end subroutine set_prof_perp_astrl_discrete

! compute normalization so that the a0 of the pulse will equal the a0 of the input deck
! at the specified xi_norm, z_norm
! the requested value of a0 is used as input, then the variable a0 stores the required
! normalization to achieve the originally requested a0.
subroutine normalize_a0_astrl_discrete( r_norm, xi_norm, z_norm, &
     k0, chirp_coefs, prof_pars, math_funcs, mode_norm, a0 ) 

  implicit none 

  real, intent(in) :: r_norm, xi_norm, z_norm, k0
  real, dimension(:), intent(in) :: chirp_coefs
  type(kw_list), intent(inout) :: prof_pars
  type(t_fparser), dimension(:), pointer, intent(inout) :: math_funcs
  integer, intent(in) :: mode_norm
  real, intent(inout) :: a0 

  real :: k
  integer :: l 
  logical :: if_norm
  real :: ar_re, ar_im, ai_re, ai_im, a0_unnormalized

  ! longitudinal frequency chirp
  k = k0
  do l = 1, size(chirp_coefs)
     k = k + chirp_coefs(l) * xi_norm ** (l - 1)
  enddo
  
  call get_prof_perp_astrl_discrete( r_norm, xi_norm, z_norm, k, k0, prof_pars, math_funcs, &
       0, ar_re, ar_im, ai_re, ai_im )
     
  a0_unnormalized = sqrt( ar_re**2 + ai_re**2 ) 
  
  a0 = a0 / a0_unnormalized
  
end subroutine normalize_a0_astrl_discrete

subroutine get_prof_perp_astrl_discrete( r, z, t, k, k0, prof_pars, math_funcs, mode, ar_re, ar_im, ai_re, ai_im )

  implicit none
  real, intent(in) :: r, z, t, k, k0
  type(kw_list), intent(in) :: prof_pars
  integer, intent(in) :: mode
  real, intent(out) :: ar_re, ar_im, ai_re, ai_im
  character(len=:), allocatable :: read_str
  real :: a0,w0, zr, curv, f_dist, gouy_shift, z_shift, z2, zr2, w, phase, r2, amp
  real(p_k_fparse), dimension(1) :: fparser_arr 
  type(t_fparser), dimension(:), pointer, intent(inout) :: math_funcs
  real :: pulselet_delay, pulselet_range, pulselet_offset, xi_pulselet_min, xi_pulselet, pulselet_weight
  integer :: j, npulselets

  ar_re = 0.0
  ar_im = 0.0
  ai_re = 0.0
  ai_im = 0.0
  
  call prof_pars%get( 'pulselet_range', pulselet_range )
  call prof_pars%get( 'pulselet_delay', pulselet_delay )
  call prof_pars%get( 'pulselet_offset', pulselet_offset )

  ! snap to grid of pulselets 
  npulselets = ceiling( 2 * pulselet_range / pulselet_delay ) + 1
  xi_pulselet_min = floor( z / pulselet_delay ) * pulselet_delay &
       - (npulselets/2)*pulselet_delay + mod( pulselet_offset, pulselet_delay ) 

  if ( mode == 0 ) then

     do j = 1, npulselets
        
        xi_pulselet = xi_pulselet_min + (j-1) * pulselet_delay

        fparser_arr(1) = xi_pulselet
        w0 = eval(math_funcs(1), fparser_arr)
        a0 = eval(math_funcs(2), fparser_arr)
        f_dist = eval(math_funcs(3), fparser_arr)

        fparser_arr(1) = z - xi_pulselet
        pulselet_weight = eval(math_funcs(4), fparser_arr)

        ! move z_shift by t. t=0 is always used for initialization. t != 0 is used for normalization.
        z_shift = -1.0 * (z + f_dist - t)
        r2 = r * r
        z2 = z_shift * z_shift
        zr = 0.5 * k * w0 * w0
        zr2 = zr * zr
        curv = z_shift / ( z2 + zr2 )
        w = w0 * sqrt( 1.0 + z2 / zr2 )
        gouy_shift = atan2( z_shift, zr )
        phase = 0.5 * k * r2 * curv - gouy_shift - (k - k0) * z
        amp = w0 / w * exp(-r2 / (w*w))

        ! accumulate contribution from this pulse 
        ar_re = ar_re + pulselet_weight * a0 * amp * cos(phase)
        ai_re = ai_re - pulselet_weight * a0 * amp * sin(phase)
        
     enddo

  else

    ar_re = 0.0
    ar_im = 0.0
    ai_re = 0.0
    ai_im = 0.0

  endif

end subroutine get_prof_perp_astrl_discrete


! ------------------------------------------------------------------------------
! CONST LONGITUDINAL PROFILES (used with astrl pulses)
! ------------------------------------------------------------------------------
subroutine set_prof_const( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  type(kw_list), intent(inout) :: prof_pars

  real :: val

  ! do nothing

end subroutine set_prof_const

subroutine get_prof_const( z, prof_pars, env )

  implicit none
  real, intent(in) :: z
  type(kw_list), intent(in) :: prof_pars
  real, intent(out) :: env

  env = 1.0

end subroutine get_prof_const

! ------------------------------------------------------------------------------
! SIN2 LONGITUDINAL PROFILES
! ------------------------------------------------------------------------------
subroutine set_prof_lon_sin2( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  type(kw_list), intent(inout) :: prof_pars

  real :: val

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

  real :: t_rise, t_fall, t_flat, flat_start, flat_end
  real, parameter :: pih = 1.570796326794897

  call prof_pars%get( 't_rise', t_rise )
  call prof_pars%get( 't_flat', t_flat )
  call prof_pars%get( 't_fall', t_fall )
  flat_start = - 0.5 * t_flat
  flat_end   =   0.5 * t_flat

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

  real :: t_rise, t_fall, t_flat, flat_start, flat_end, t     

  call prof_pars%get( 't_rise', t_rise )
  call prof_pars%get( 't_flat', t_flat )
  call prof_pars%get( 't_fall', t_fall )
  flat_start = - 0.5 * t_flat
  flat_end   =   0.5 * t_flat

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


! ------------------------------------------------------------------------------
! ARBITRARY LONGITUDINAL PROFILES PIECEWISE LINEAR
! ------------------------------------------------------------------------------
subroutine set_prof_lon_pw_linear( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars

  integer :: len_t, i
  real, dimension(:), allocatable :: ft, t

  call input%get( trim(sect_name)//'.piecewise_ft', ft )
  call input%get( trim(sect_name)//'.piecewise_t', t )

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  len_t = size(t)
  do i = 2, len_t
    if ( t(i) <= t(i-1) ) then
      call write_err( 'Array of longitudinal position must be monotonically increasing!' )
    endif
  enddo

  allocate( prof_pars( 2 * len_t ) )

  do i = 1, len_t
    prof_pars(i) = t(i)
    prof_pars(len_t+i) = ft(i)
  enddo

  deallocate( ft, t )

end subroutine set_prof_lon_pw_linear

subroutine get_prof_lon_pw_linear( z, prof_pars, env )

  implicit none
  real, intent(in) :: z
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: env
  
  integer :: len_t, i
  real, dimension(:), pointer :: t_array => null(), ft_array => null()

  len_t = size( prof_pars ) / 2

  ! use pointer to associate with the segment of data to avoid hard-copy
  t_array  => prof_pars( 1 : len_t )
  ft_array => prof_pars( len_t+1 : 2*len_t )

  if ( z < t_array(1) .or. z > t_array(len_t) ) then
    env = 0.0
  else 
    env = 0.0
    do i = 2, len_t
      if ( z <= t_array(i) ) then
        env = ft_array(i-1) + ( ft_array(i) - ft_array(i-1) ) / &
          ( t_array(i) - t_array(i-1) ) * ( z - t_array(i-1) )
        exit
      endif
    enddo
  endif

end subroutine get_prof_lon_pw_linear

! ------------------------------------------------------------------------------
! ARBITRARY LONGITUDINAL PROFILES CUBIC SPLINE
! ------------------------------------------------------------------------------
subroutine set_prof_lon_cubic_spline( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars

  integer :: i, n
  real, dimension(:), allocatable :: t, ft, diag, off_diag_lower, off_diag_upper, rhs, d
  character(len=:), allocatable :: kind

  call input%get( trim(sect_name)//'.piecewise_ft', ft )
  call input%get( trim(sect_name)//'.piecewise_t', t )
  call input%get( trim(sect_name)//'.endpoints_type', kind )

  n = size(ft)
  do i = 2, n
    if ( t(i) <= t(i-1) ) then
      call write_err( 'Array of longitudinal position must be monotonically increasing!' )
    endif
  enddo
  ! note that the first element of off_diag_lower and the last element of off_diag_upper are useless.
  allocate(d(n), diag(n), off_diag_lower(n), off_diag_upper(n), rhs(n))
  diag = 4.0
  off_diag_lower = 1.0
  off_diag_upper = 1.0
  rhs(2:n-1) = 3.0 * (ft(3:n) - ft(1:n-2))

  ! solve the derivatives of cubic spline function at each points
  select case (kind)
    case ('natural')
      ! The second-order derivatives of the endpoints are zero
      diag(1) = 2.0; diag(n) = 2.0
      rhs(1) = 3.0 * (ft(2) - ft(1))
      rhs(n) = 3.0 * (ft(n) - ft(n-1))
      call tdma(n, off_diag_lower, diag, off_diag_upper, rhs, d)
    case ('clamped')
      ! The first-order derivatives of the endpoints are zero
      call tdma(n-2, off_diag_lower(2:n-1), diag(2:n-1), off_diag_upper(2:n-1), rhs(2:n-1), d(2:n-1))
      d(1) = 0.0; d(n) = 0.0
    case ('not-a-knot')
      ! The third-order derivatives of the endpoints are equal to those of the adjacent points.
      diag(1) = 2.0; diag(n) = 2.0
      off_diag_lower(n) = 4.0
      off_diag_upper(1) = 4.0
      rhs(1) = -5.0 * ft(1) + 4.0 * ft(2) + ft(3)
      rhs(n) = -5.0 * ft(n-2) + 4.0 * ft(n-1) + ft(n)
      call tdma(n, off_diag_lower, diag, off_diag_upper, rhs, d)
    case default
      call write_err("Invaild endpoints type of cubic spline profile.")
  end select

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )
  allocate( prof_pars(3 * n) )
  prof_pars(    1:  n) = t
  prof_pars(  n+1:2*n) = ft
  prof_pars(2*n+1:3*n) = d

end subroutine set_prof_lon_cubic_spline

subroutine get_prof_lon_cubic_spline( z, prof_pars, env )

  implicit none
  real, intent(in) :: z
  real, intent(in), dimension(:), pointer :: prof_pars
  real, intent(out) :: env
  
  integer :: i, n
  real :: t_norm, a, b, c, d
  real, dimension(:), pointer :: t_ptr => null(), ft_ptr => null(), dft_ptr => null()

  n = size(prof_pars) / 3
  ! use pointer to associate with the segment of data to avoid hard-copy
  t_ptr   => prof_pars(    1:  n)
  ft_ptr  => prof_pars(  n+1:2*n)
  dft_ptr => prof_pars(2*n+1:3*n)

  env = 0.0
  if ( z <= t_ptr(1) .or. z > t_ptr(n) ) then
    env = 0.0
  else
    env = 0.0
    do i = 1, n-1
      if ( z <= t_ptr(i+1) ) then
        t_norm = (z - t_ptr(i)) / (t_ptr(i+1) - t_ptr(i))
        env = ft_ptr(i) + dft_ptr(i) * t_norm &
          + ( 3.0 * (ft_ptr(i+1) - ft_ptr(i)) - 2.0 * dft_ptr(i) - dft_ptr(i+1) ) * t_norm**2 &
          + ( 2.0 * (ft_ptr(i) - ft_ptr(i+1)) + dft_ptr(i) + dft_ptr(i+1) ) * t_norm**3
        exit
      endif
    enddo
  endif

end subroutine get_prof_lon_cubic_spline

end module profile_laser_lib
