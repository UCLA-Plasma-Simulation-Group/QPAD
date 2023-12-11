module profile_laser_class

! use parallel_module
use sysutil_module
use ufield_class
use kwargs_class
use input_class
use options_class
use profile_laser_lib
use m_fparser

implicit none
private

type, public :: profile_laser

  private

  ! Intensity profiles in perpendicular and longitudinal directions
  integer, dimension(2) :: prof_type

  ! Initialization geometry
  ! integer :: geom

  ! Grid offset and number of cells
  integer :: noff_r, noff_z, nr, nz, nrp, nzp

  ! Cell size
  real :: dr, dz

  ! Position of the lower boundary in the longitudinal direction
  real :: z0

  ! Longitudinal central position
  real :: lon_center

  ! Global intensity
  real :: a0

  ! central frequency
  real :: k0

  ! astrl math functions
  type(t_fparser), dimension(:), pointer :: math_funcs => null()
  ! functions by order are w0_math_func, a0_math_func, and s0_math_func

  ! frequency chirp
  real, dimension(:), allocatable :: chirp_coefs

  ! Parameter list of the intensity profiles in 1, 2, 3 directions
  type( kw_list ) :: prof_perp_pars, prof_lon_pars

  ! Parameter list of longitudinal intensity profile in piecewise linear/cubic spline
  real, dimension(:), pointer :: prof_lon_pars_lin => null()
  real, dimension(:), pointer :: prof_lon_pars_cub => null()

  ! Procedure of setting the density profiles
  procedure( set_prof_intf ), nopass, pointer :: set_prof_perp => null()
  procedure( set_prof_intf ), nopass, pointer :: set_prof_lon => null()
  procedure( set_prof_perp_astrl_analytic_intf ), nopass, pointer :: set_prof_perp_astrl_analytic => null()
  procedure( set_prof_array_intf ), nopass, pointer :: set_prof_lon_lin => null()
  procedure( set_prof_array_intf ), nopass, pointer :: set_prof_lon_cub => null()

  ! Procedure of getting the beam density
  procedure( get_prof_perp_intf ), nopass, pointer :: get_prof_perp => null()
  procedure( get_prof_perp_astrl_analytic_intf ), nopass, pointer :: get_prof_perp_astrl_analytic => null()
  procedure( get_prof_lon_intf ), nopass, pointer :: get_prof_lon => null()
  procedure( get_prof_array_intf ), nopass, pointer :: get_prof_lon_lin => null()
  procedure( get_prof_array_intf ), nopass, pointer :: get_prof_lon_cub => null()
  contains

  procedure :: new => init_profile_laser
  procedure :: del => end_profile_laser
  procedure :: launch => launch_profile_laser
  procedure :: norm_power => norm_power_laser 

end type profile_laser

interface
  subroutine get_prof_perp_intf( r, z, k, k0, prof_pars, mode, ar_re, ar_im, ai_re, ai_im )
    import kw_list
    implicit none
    real, intent(in) :: r, z, k, k0
    type(kw_list), intent(in) :: prof_pars
    integer, intent(in) :: mode
    real, intent(out) :: ar_re, ar_im, ai_re, ai_im
  end subroutine get_prof_perp_intf

  subroutine get_prof_lon_intf( z, prof_pars, envelope )
    import kw_list
    implicit none
    real, intent(in) :: z
    type(kw_list), intent(in) :: prof_pars
    real, intent(out) :: envelope
  end subroutine get_prof_lon_intf

  subroutine get_prof_array_intf( z, prof_pars, envelope)
    implicit none
    real, intent(in) :: z
    real, intent(in), dimension(:), pointer :: prof_pars
    real, intent(out) :: envelope
  end subroutine get_prof_array_intf

  subroutine set_prof_intf( input, sect_name, prof_pars )
    import input_json, kw_list
    implicit none
    type( input_json ), intent(inout) :: input
    character(len=*), intent(in) :: sect_name
    type(kw_list), intent(inout) :: prof_pars
  end subroutine set_prof_intf

  subroutine set_prof_array_intf( input, sect_name, prof_pars )
    import input_json
    implicit none
    type( input_json ), intent(inout) :: input
    character(len=*), intent(in) :: sect_name
    real, intent(inout), dimension(:), pointer :: prof_pars
  end subroutine set_prof_array_intf

  subroutine get_prof_perp_astrl_analytic_intf( r, z, t, k, k0, prof_pars, math_funcs, mode, ar_re, ar_im, ai_re, ai_im )
    import kw_list, t_fparser
    implicit none
    real, intent(in) :: r, z, t, k, k0
    type(kw_list), intent(in) :: prof_pars
    integer, intent(in) :: mode
    real, intent(out) :: ar_re, ar_im, ai_re, ai_im
    type(t_fparser), dimension(:), pointer, intent(inout) :: math_funcs
  end subroutine get_prof_perp_astrl_analytic_intf

  subroutine set_prof_perp_astrl_analytic_intf( input, sect_name, prof_pars, math_funcs )
    import input_json, kw_list, t_fparser
    implicit none
    type( input_json ), intent(inout) :: input
    character(len=*), intent(in) :: sect_name
    type(kw_list), intent(inout) :: prof_pars
    type(t_fparser), dimension(:), pointer, intent(inout) :: math_funcs
  end subroutine set_prof_perp_astrl_analytic_intf



end interface

integer, parameter, public :: p_prof_laser_gaussian = 0
integer, parameter, public :: p_prof_laser_laguerre = 1
integer, parameter, public :: p_prof_laser_astrl_analytic = 2
integer, parameter, public :: p_prof_laser_astrl_discrete = 3
integer, parameter, public :: p_prof_laser_sin2 = 100
integer, parameter, public :: p_prof_laser_poly = 101
integer, parameter, public :: p_prof_laser_const = 102
integer, parameter, public :: p_prof_laser_pw_linear = 1000
integer, parameter, public :: p_prof_laser_cubic_spline = 1001

character(len=32), save :: cls_name = 'profile_laser'
integer, save :: cls_level = 2

contains

subroutine init_profile_laser( this, input, opts, sect_id )

  implicit none
  class( profile_laser ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  integer, intent(in) :: sect_id

  character(len=:), allocatable :: read_str
  character(len=20) :: sect_name
  character(len=32), save :: sname = 'init_profile_laser'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%noff_r = opts%get_noff(1)
  this%noff_z = opts%get_noff(2)
  this%nr     = opts%get_nd(1)
  this%nz     = opts%get_nd(2)
  this%nrp    = opts%get_ndp(1)
  this%nzp    = opts%get_ndp(2)
  this%dr     = opts%get_dr()
  this%dz     = opts%get_dxi()

  sect_name = 'laser(' // num2str(sect_id) // ')'

  ! read and set profile types
  call input%get( trim(sect_name) // '.profile(1)', read_str )
  select case ( trim(read_str) )

    case ( 'gaussian' )
      this%prof_type(1)  = p_prof_laser_gaussian
      this%set_prof_perp => set_prof_perp_gaussian
      this%get_prof_perp => get_prof_perp_gaussian

    case ( 'laguerre' )
      this%prof_type(1)  = p_prof_laser_laguerre
      this%set_prof_perp => set_prof_perp_laguerre
      this%get_prof_perp => get_prof_perp_laguerre

    case ( 'astrl_analytic' )
      this%prof_type(1)  = p_prof_laser_astrl_analytic
      this%set_prof_perp_astrl_analytic => set_prof_perp_astrl_analytic
      this%get_prof_perp_astrl_analytic => get_prof_perp_astrl_analytic

    case ( 'astrl_discrete' )
      this%prof_type(1)  = p_prof_laser_astrl_discrete
      this%set_prof_perp_astrl_analytic => set_prof_perp_astrl_discrete
      this%get_prof_perp_astrl_analytic => get_prof_perp_astrl_discrete

    case default
      call write_err( 'Invalid intensity profile in direction 1! &
        &Currently available include "gaussian" and "laguerre".' )

  end select

  call input%get( trim(sect_name) // '.profile(2)', read_str )
  select case ( trim(read_str) )

    case ( 'sin2' )
      this%prof_type(2)  = p_prof_laser_sin2
      this%set_prof_lon => set_prof_lon_sin2
      this%get_prof_lon => get_prof_lon_sin2

    case ( 'polynomial' )
      this%prof_type(2)  = p_prof_laser_poly
      this%set_prof_lon => set_prof_lon_poly
      this%get_prof_lon => get_prof_lon_poly

    case ('piecewise_linear')
      this%prof_type(2)  = p_prof_laser_pw_linear
      this%set_prof_lon_lin => set_prof_lon_pw_linear
      this%get_prof_lon_lin => get_prof_lon_pw_linear
    
    case ('cubic_spline')
      this%prof_type(2)  = p_prof_laser_cubic_spline
      this%set_prof_lon_cub => set_prof_lon_cubic_spline
      this%get_prof_lon_cub => get_prof_lon_cubic_spline

   case ( 'const' )
      this%prof_type(2)  = p_prof_laser_const
      this%set_prof_lon => set_prof_const
      this%get_prof_lon => get_prof_const
    
    case default
      call write_err( 'Invalid intensity profile in direction 1! &
        &Currently available include "sin2", "polynomial", "piecewise_linear"and"cubic_spline".' )

  end select

  ! read and store the profile parameters into the parameter lists
  if(this%prof_type(1) == p_prof_laser_astrl_analytic) then
     call this%set_prof_perp_astrl_analytic( input, trim(sect_name), this%prof_perp_pars, this%math_funcs )
  else if( this%prof_type(1) == p_prof_laser_astrl_discrete ) then
     call this%set_prof_perp_astrl_analytic( input, trim(sect_name), this%prof_perp_pars, this%math_funcs )
  else
     call this%set_prof_perp( input, trim(sect_name), this%prof_perp_pars )
  endif
  
  if ( this%prof_type(2) == p_prof_laser_pw_linear ) then
    call this%set_prof_lon_lin( input, trim(sect_name), this%prof_lon_pars_lin )
  else if ( this%prof_type(2) ==  p_prof_laser_cubic_spline ) then
    call this%set_prof_lon_cub( input, trim(sect_name), this%prof_lon_pars_cub ) 
  else 
    call this%set_prof_lon( input, trim(sect_name), this%prof_lon_pars )
  endif

  if(this%prof_type(1) /= p_prof_laser_astrl_analytic .and. this%prof_type(1) /= p_prof_laser_astrl_discrete ) then
     call input%get( trim(sect_name) // '.a0', this%a0 )
  else 
     this%a0 = 1.0
  endif
  call input%get( trim(sect_name) // '.k0', this%k0 )
  call input%get( trim(sect_name) // '.lon_center', this%lon_center )
  call input%get( 'simulation.box.z(1)', this%z0 )

  this%chirp_coefs = [0.0]
  if ( input%found( trim(sect_name) // '.chirp_coefs' ) ) then
    call input%get( trim(sect_name) // '.chirp_coefs', this%chirp_coefs )
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_profile_laser

subroutine end_profile_laser( this )

  implicit none
  class( profile_laser ), intent(inout) :: this
  character(len=32), save :: sname = 'end_fdist3d_std'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  ! do nothing
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_profile_laser

subroutine launch_profile_laser( this, ar_re, ar_im, ai_re, ai_im )

  implicit none
  class( profile_laser ), intent(inout) :: this
  class( ufield ), intent(inout), dimension(:), pointer :: ar_re, ar_im, ai_re, ai_im

  integer :: i, j, m, l, max_mode
  real :: r, z, t, k, arr, ari, air, aii, env
  character(len=32), save :: sname = 'launch_profile_laser'
  logical :: if_norm_power
  real :: z_norm, power_norm 
  
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ! launch occurs at t=0 
  t = 0.0 
  
  max_mode = size(ar_im)
  do m = 0, max_mode
    do j = 1, this%nzp
      ! note that here "z" refers to xi = t - z
      z = ( this%noff_z + j - 1 ) * this%dz + this%z0 - this%lon_center

      ! longitudinal frequency chirp
      k = this%k0
      do l = 1, size(this%chirp_coefs)
        k = k + this%chirp_coefs(l) * z ** l
      enddo

      do i = 1, this%nrp
        r = ( this%noff_r + i - 1 ) * this%dr
        if(this%prof_type(1) == p_prof_laser_astrl_analytic) then
           call this%get_prof_perp_astrl_analytic( r, z, t, k, this%k0, &
                this%prof_perp_pars, this%math_funcs, m, arr, ari, air, aii )
        else if (this%prof_type(1) == p_prof_laser_astrl_discrete) then
           call this%get_prof_perp_astrl_analytic( r, z, t, k, this%k0, &
                this%prof_perp_pars, this%math_funcs, m, arr, ari, air, aii )           
        else
           call this%get_prof_perp( r, z, k, this%k0, this%prof_perp_pars, m, arr, ari, air, aii )
        endif

        if ( this%prof_type(2) == p_prof_laser_pw_linear ) then
          call this%get_prof_lon_lin( z, this%prof_lon_pars_lin, env )
        else if ( this%prof_type(2) ==  p_prof_laser_cubic_spline ) then 
          call this%get_prof_lon_cub( z, this%prof_lon_pars_cub, env )
        else 
          call this%get_prof_lon( z, this%prof_lon_pars, env )
       endif
        
        env = env * this%a0
        ar_re(m)%f2( 1, i, j ) = env * arr
        ai_re(m)%f2( 1, i, j ) = env * air
        if ( m == 0 ) cycle
        ar_im(m)%f2( 1, i, j ) = env * ari
        ai_im(m)%f2( 1, i, j ) = env * aii
      enddo
    enddo
  enddo

  ! normalize power to match the desired value 
  if( this%prof_type(1) == p_prof_laser_astrl_discrete ) then

     call this%prof_perp_pars%get( 'if_norm_power', if_norm_power )

     if( if_norm_power ) then 
        
        call this%prof_perp_pars%get( 'power_norm', power_norm )
        call this%norm_power( ar_re, ar_im, ai_re, ai_im, z_norm, power_norm )
        
        call this%prof_perp_pars%get( 'z_norm', z_norm )

     endif

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine launch_profile_laser

! perform a numerical integration to make the total power through the transverse slice
! of the nearest gridpoint to z_norm equal to power_norm
subroutine norm_power_laser( this, ar_re, ar_im, ai_re, ai_im, z_norm, power_norm )

  implicit none
  class( profile_laser ), intent(inout) :: this
  class( ufield ), intent(inout), dimension(:), pointer :: ar_re, ar_im, ai_re, ai_im
  real, intent(in) :: z_norm, power_norm 

  integer :: i, j_norm
  real :: r 
  real :: power_unnormalized, scale 
  
  ! get j index for integration, from z = ( this%noff_z + j - 1 ) * this%dz + this%z0
  j_norm = floor( ( z_norm - this%z0 ) / this%dz ) - this%noff_z + 1 
  
  ! numerically integrate the power (only implemented for mode 0)
  power_unnormalized = 0 

  do i = 1, this%nrp
     r = ( this%noff_r + i - 1 ) * this%dr     
     power_unnormalized = power_unnormalized + (1 / (8 * 3.1415) ) * (2 * 3.1415 * r) * &
          ( ar_re(0)%f2( 1, i, j_norm )**2 + ai_re(0)%f2( 1, i, j_norm )**2 )

  enddo

  scale = sqrt( power_norm / power_unnormalized ) 

  ar_re(0)%f2(:,:,:) = ar_re(0)%f2(:,:,:) * scale
  ar_im(0)%f2(:,:,:) = ar_im(0)%f2(:,:,:) * scale
  ai_re(0)%f2(:,:,:) = ai_re(0)%f2(:,:,:) * scale
  ai_im(0)%f2(:,:,:) = ai_im(0)%f2(:,:,:) * scale  
  
end subroutine norm_power_laser

end module profile_laser_class
