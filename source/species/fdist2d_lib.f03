module fdist2d_lib

use input_class
use sysutil_module
use param

implicit none

private

public :: set_prof_perp_uniform, get_den_perp_uniform
public :: set_prof_perp_para_chl, get_den_perp_para_chl
public :: set_prof_perp_hllw_chl, get_den_perp_hllw_chl
public :: set_prof_lon_uniform, get_den_lon_uniform
public :: set_prof_lon_pw_linear, get_den_lon_pw_linear
public :: set_prof_lon_sine, get_den_lon_sine

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

subroutine get_den_perp_uniform( x, s, prof_pars_perp, prof_pars_lon, den_value )
  implicit none
  real, intent(in), dimension(2) :: x
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_perp, prof_pars_lon
  real, intent(out) :: den_value
  den_value = 1.0
end subroutine get_den_perp_uniform

subroutine set_prof_perp_para_chl( input, sect_name, prof_pars )

  implicit none
  type( input_json ), intent(inout) :: input
  character(len=*), intent(in) :: sect_name
  real, intent(inout), dimension(:), pointer :: prof_pars

  ! this should never be called
  if ( associated(prof_pars) ) deallocate( prof_pars )

  allocate( prof_pars(4) )

  call input%get( trim(sect_name) // '.channel_n0', prof_pars(1) )
  call input%get( trim(sect_name) // '.channel_depth', prof_pars(2) )
  call input%get( trim(sect_name) // '.channel_r0', prof_pars(3) )
  call input%get( trim(sect_name) // '.channel_width', prof_pars(4) )

end subroutine set_prof_perp_para_chl

subroutine get_den_perp_para_chl( x, s, prof_pars_perp, prof_pars_lon, den_value )

  implicit none
  real, intent(in), dimension(2) :: x
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_perp, prof_pars_lon
  real, intent(out) :: den_value

  real :: r2, n0, n_depth, r02, r2_max

  n0      = prof_pars_perp(1)
  n_depth = prof_pars_perp(2)
  r02     = prof_pars_perp(3)**2
  r2_max  = prof_pars_perp(4)**2
  r2      = x(1)**2 + x(2)**2

  if ( r2 < r2_max ) then
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

subroutine get_den_perp_hllw_chl( x, s, prof_pars_perp, prof_pars_lon, den_value )

  implicit none
  real, intent(in), dimension(2) :: x
  real, intent(in) :: s
  real, intent(in), dimension(:), pointer :: prof_pars_perp, prof_pars_lon
  real, intent(out) :: den_value

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

  ! check boundaries
  if ( s < s_array(1) .or. s > s_array(len_s) ) then
    call write_err( 'Longitudinal density profile at current position is undefined!' )
  endif

  do i = 2, len_s
    if ( s_array(i) <= s_array(i-1) ) then
      call write_err( 'Array of longitudinal position must be monotonically increasing!' )
    endif

    if ( s <= s_array(i) ) then
      den_value = fs_array(i-1) + ( fs_array(i) - fs_array(i-1) ) / &
        ( s_array(i) - s_array(i-1) ) * ( s - s_array(i-1) )
      exit
    endif
  enddo

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

end subroutine get_den_lon_sine

end module fdist2d_lib