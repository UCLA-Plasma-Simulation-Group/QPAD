module fdist2d_lib

use input_class
use sysutil_module
use param
use math_module
use m_fparser
implicit none

private

public :: set_prof_perp_uniform, get_den_perp_uniform
public :: set_prof_perp_para_chl, get_den_perp_para_chl
public :: set_prof_perp_hllw_chl, get_den_perp_hllw_chl
public :: set_prof_lon_uniform, get_den_lon_uniform
public :: set_prof_lon_pw_linear, get_den_lon_pw_linear
public :: set_prof_lon_sine, get_den_lon_sine
public :: set_prof_lon_cubic_spline, get_den_lon_cubic_spline
public :: set_prof_analytic
public :: get_den_perp_analytic, get_den_lon_analytic

contains

! ------------------------------------------------------------------------------
! SECTION OF PERPENDICULAR PROFILES
! ------------------------------------------------------------------------------

subroutine set_prof_perp_uniform( input, sect_name, prof_pars )
  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars
  ! placeholder for uniform profile, do nothing
  return
end subroutine set_prof_perp_uniform

subroutine get_den_perp_uniform( x, s, prof_pars_perp, prof_pars_lon, den_value , math_func )
  implicit none
  real, intent(in), dimension(2) :: x
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_perp, prof_pars_lon
  real, intent(out) :: den_value
  type(t_fparser), pointer, intent(inout) :: math_func
  den_value = 1.0
end subroutine get_den_perp_uniform

subroutine set_prof_perp_para_chl( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  allocate( prof_pars(5) )

  call input%get( trim(sect_name) // '.channel_n0', prof_pars(1) )
  call input%get( trim(sect_name) // '.channel_depth', prof_pars(2) )
  call input%get( trim(sect_name) // '.channel_r0', prof_pars(3) )
  call input%get( trim(sect_name) // '.channel_width', prof_pars(4) )
  call input%get( trim(sect_name) // '.channel_rmax', prof_pars(5) )
end subroutine set_prof_perp_para_chl

subroutine get_den_perp_para_chl( x, s, prof_pars_perp, prof_pars_lon, den_value , math_func )

  implicit none
  real, intent(in), dimension(2) :: x
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_perp, prof_pars_lon
  real, intent(out) :: den_value
  type(t_fparser), pointer, intent(inout) :: math_func

  real :: r2, n0, n_depth, r02, r2_max

  real :: rlim
  n0      = prof_pars_perp(1)
  n_depth = prof_pars_perp(2)
  r02     = prof_pars_perp(3)**2
  r2_max  = prof_pars_perp(4)**2
  r2      = x(1)**2 + x(2)**2

  rlim = prof_pars_perp(5)

  if(r2 > rlim **2) then 
    den_value = 0
  else if ( r2 < r2_max ) then
    den_value = n0 + n_depth * r2 / r02
  else
    den_value = n0 + n_depth * r2_max / r02
  endif

end subroutine get_den_perp_para_chl

subroutine set_prof_perp_hllw_chl( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  allocate( prof_pars(3) )

  call input%get( trim(sect_name)//'.channel_rmin', prof_pars(1) )
  call input%get( trim(sect_name)//'.channel_rmax', prof_pars(2) )
  call input%get( trim(sect_name)//'.channel_depth', prof_pars(3) )

end subroutine set_prof_perp_hllw_chl

subroutine get_den_perp_hllw_chl( x, s, prof_pars_perp, prof_pars_lon, den_value , math_func )

  implicit none
  real, intent(in), dimension(2) :: x
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_perp, prof_pars_lon
  real, intent(out) :: den_value
  type(t_fparser), pointer, intent(inout) :: math_func

  real :: rmin, rmax, n_depth, r

  rmin    = prof_pars_perp(1)
  rmax    = prof_pars_perp(2)
  n_depth = prof_pars_perp(3)
  r       = sqrt( x(1)**2 + x(2)**2 )

  if ( r < rmax .and. r > rmin ) then
    den_value = n_depth
  else
    den_value = 0.0
  endif

end subroutine get_den_perp_hllw_chl

subroutine set_prof_analytic( input, sect_name, prof_pars, math_func )
  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars
  type(t_fparser), pointer, intent(inout) :: math_func
  character(len=:), allocatable :: read_str
  integer :: ierr

  if ( .not. associated( math_func ) ) then
    allocate( t_fparser :: math_func )
  endif
  call input%get( trim(sect_name) // '.math_func', read_str )
  call setup(math_func, trim(read_str), (/'x','y','z'/), ierr)
  

end subroutine set_prof_analytic

subroutine get_den_perp_analytic( x, s, prof_pars_perp, prof_pars_lon, den_value, math_func )
  implicit none
  real, intent(in), dimension(2) :: x
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_perp, prof_pars_lon
  real, intent(out) :: den_value

  real(p_k_fparse), dimension(3) :: fparser_arr 
  type(t_fparser), pointer, intent(inout) :: math_func

  fparser_arr(1) = x(1)
  fparser_arr(2) = x(2)
  fparser_arr(3) = s
  den_value = eval(math_func, fparser_arr)

end subroutine get_den_perp_analytic

! ------------------------------------------------------------------------------
! SECTION OF LONGITUDINAL PROFILES
! ------------------------------------------------------------------------------

subroutine set_prof_lon_uniform( input, sect_name, prof_pars )
  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars
  ! placeholder for uniform profile, do nothing
  return
end subroutine set_prof_lon_uniform

subroutine get_den_lon_uniform( s, prof_pars_lon, den_value )
  implicit none
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_lon
  real, intent(out) :: den_value
  den_value = 1.0
end subroutine get_den_lon_uniform

subroutine get_den_lon_analytic( s, prof_pars_lon, den_value )
  implicit none
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_lon
  real, intent(out) :: den_value
  den_value = 1.0
end subroutine get_den_lon_analytic

subroutine set_prof_lon_pw_linear( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars

  integer :: len_s, i
  real, dimension(:), allocatable :: fs, s

  call input%get( trim(sect_name)//'.piecewise_fs', fs )
  call input%get( trim(sect_name)//'.piecewise_s', s )

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  len_s = size(s)
  do i = 2, len_s
    if ( s(i) <= s(i-1) ) then
      call write_err( 'Array of longitudinal position must be monotonically increasing!' )
    endif
  enddo

  allocate( prof_pars( 2 * len_s ) )

  do i = 1, len_s
    prof_pars(i) = s(i)
    prof_pars(len_s+i) = fs(i)
  enddo

  deallocate( fs, s )

end subroutine set_prof_lon_pw_linear

subroutine get_den_lon_pw_linear( s, prof_pars_lon, den_value )

  implicit none
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_lon
  real, intent(out) :: den_value
  
  integer :: len_s, i
  real, dimension(:), pointer :: s_array => null(), fs_array => null()

  len_s = size( prof_pars_lon ) / 2

  ! use pointer to associate with the segment of data to avoid hard-copy
  s_array  => prof_pars_lon( 1 : len_s )
  fs_array => prof_pars_lon( len_s+1 : 2*len_s )

  if ( s <= s_array(1) .or. s > s_array(len_s) ) then
    den_value = 0.0
  else
    den_value = 0.0
    do i = 2, len_s
      if ( s <= s_array(i) ) then
        den_value = fs_array(i-1) + ( fs_array(i) - fs_array(i-1) ) / &
          ( s_array(i) - s_array(i-1) ) * ( s - s_array(i-1) )
        exit
      endif
    enddo
  endif

end subroutine get_den_lon_pw_linear

subroutine set_prof_lon_sine( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars

  integer :: i
  real :: start, length
  real, dimension(:), allocatable :: coefs

  call input%get( trim(sect_name)//'.start', start )
  call input%get( trim(sect_name)//'.length', length )
  call input%get( trim(sect_name)//'.sine_coefs', coefs )

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  allocate( prof_pars( 2 * size(coefs) + 2 ) )

  prof_pars(1) = start
  prof_pars(2) = length
  prof_pars(3:) = coefs

end subroutine set_prof_lon_sine

subroutine get_den_lon_sine( s, prof_pars_lon, den_value )

  implicit none
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_lon
  real, intent(out) :: den_value
  
  integer :: i
  real, dimension(:), pointer :: coefs => null()
  real :: start, length

  start = prof_pars_lon(1)
  length = prof_pars_lon(2)
  ! use pointer to associate with the segment of data to avoid hard-copy
  coefs => prof_pars_lon(3:)

  den_value = 0.0
  if (s > start .and. s <= start + length) then
    do i = 1, size(coefs)
      den_value = den_value + coefs(i) * sin(pi * real(i) * (s - start) / length)
    enddo
  endif
  den_value = abs(den_value)

end subroutine get_den_lon_sine

subroutine set_prof_lon_cubic_spline( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars

  integer :: i, n
  real, dimension(:), allocatable :: s, fs, diag, off_diag_lower, off_diag_upper, rhs, d
  character(len=:), allocatable :: kind

  call input%get( trim(sect_name)//'.piecewise_fs', fs )
  call input%get( trim(sect_name)//'.piecewise_s', s )
  call input%get( trim(sect_name)//'.endpoints_type', kind )

  n = size(fs)
  do i = 2, n
    if ( s(i) <= s(i-1) ) then
      call write_err( 'Array of longitudinal position must be monotonically increasing!' )
    endif
  enddo
  ! note that the first element of off_diag_lower and the last element of off_diag_upper are useless.
  allocate(d(n), diag(n), off_diag_lower(n), off_diag_upper(n), rhs(n))
  diag = 4.0
  off_diag_lower = 1.0
  off_diag_upper = 1.0
  rhs(2:n-1) = 3.0 * (fs(3:n) - fs(1:n-2))

  ! solve the derivatives of cubic spline function at each points
  ! according to https://mathworld.wolfram.com/CubicSpline.html
  select case (kind)
    case ('natural')
      ! The second-order derivatives of the endpoints are zero
      diag(1) = 2.0; diag(n) = 2.0
      rhs(1) = 3.0 * (fs(2) - fs(1))
      rhs(n) = 3.0 * (fs(n) - fs(n-1))
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
      rhs(1) = -5.0 * fs(1) + 4.0 * fs(2) + fs(3)
      rhs(n) = -5.0 * fs(n-2) + 4.0 * fs(n-1) + fs(n)
      call tdma(n, off_diag_lower, diag, off_diag_upper, rhs, d)
    case default
      call write_err("Invaild endpoints type of cubic spline profile.")
  end select

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )
  allocate( prof_pars(3 * n) )

  prof_pars(    1:  n) = s
  prof_pars(  n+1:2*n) = fs
  prof_pars(2*n+1:3*n) = d

end subroutine set_prof_lon_cubic_spline

subroutine get_den_lon_cubic_spline( s, prof_pars_lon, den_value )

  implicit none
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_lon
  real, intent(out) :: den_value
  
  integer :: i, n
  real :: s_norm, a, b, c, d
  real, dimension(:), pointer :: s_ptr => null(), fs_ptr => null(), dfs_ptr => null()

  n = size(prof_pars_lon) / 3
  ! use pointer to associate with the segment of data to avoid hard-copy
  s_ptr   => prof_pars_lon(    1:  n)
  fs_ptr  => prof_pars_lon(  n+1:2*n)
  dfs_ptr => prof_pars_lon(2*n+1:3*n)

  den_value = 0.0
  if ( s <= s_ptr(1) .or. s > s_ptr(n) ) then
    den_value = 0.0
  else
    den_value = 0.0
    do i = 1, n-1
      if ( s <= s_ptr(i+1) ) then
        s_norm = (s - s_ptr(i)) / (s_ptr(i+1) - s_ptr(i))
        den_value = fs_ptr(i) + dfs_ptr(i) * s_norm &
          + ( 3.0 * (fs_ptr(i+1) - fs_ptr(i)) - 2.0 * dfs_ptr(i) - dfs_ptr(i+1) ) * s_norm**2 &
          + ( 2.0 * (fs_ptr(i) - fs_ptr(i+1)) + dfs_ptr(i) + dfs_ptr(i+1) ) * s_norm**3
        exit
      endif
    enddo
  endif

end subroutine get_den_lon_cubic_spline

end module fdist2d_lib