module profile_laser_class

! use parallel_module
use sysutil_module
use ufield_class
use kwargs_class
use input_class
use options_class
use profile_laser_lib

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

  ! Global intensity
  real :: a0

  ! central frequency
  real :: k0

  ! Parameter list of the intensity profiles in 1, 2, 3 directions
  type( kw_list ) :: prof_perp_pars, prof_lon_pars

  ! Procedure of setting the density profiles
  procedure( set_prof_intf ), nopass, pointer :: set_prof_perp => null()
  procedure( set_prof_intf ), nopass, pointer :: set_prof_lon => null()

  ! Procedure of getting the beam density
  procedure( get_prof_perp_intf ), nopass, pointer :: get_prof_perp => null()
  procedure( get_prof_lon_intf ), nopass, pointer :: get_prof_lon => null()

  contains

  procedure :: new => init_profile_laser
  procedure :: del => end_profile_laser
  procedure :: launch => launch_profile_laser

end type profile_laser

interface
  subroutine get_prof_perp_intf( r, z, k0, prof_pars, mode, ar_re, ar_im, ai_re, ai_im )
    import kw_list
    implicit none
    real, intent(in) :: r, z, k0
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

  subroutine set_prof_intf( input, sect_name, prof_pars )
    import input_json, kw_list
    implicit none
    type( input_json ), intent(inout) :: input
    character(len=*), intent(in) :: sect_name
    type(kw_list), intent(inout) :: prof_pars
  end subroutine set_prof_intf
end interface

integer, parameter, public :: p_prof_laser_gaussian = 0
integer, parameter, public :: p_prof_laser_laguerre = 1
integer, parameter, public :: p_prof_laser_sin2 = 100
integer, parameter, public :: p_prof_laser_poly = 101

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

    case default
      call write_err( 'Invalid intensity profile in direction 1! &
        &Currently available include "sin2" and "polynomial".' )

  end select

  ! read and store the profile parameters into the parameter lists
  call this%set_prof_perp( input, trim(sect_name), this%prof_perp_pars )
  call this%set_prof_lon( input, trim(sect_name), this%prof_lon_pars )

  call input%get( trim(sect_name) // '.a0', this%a0 )
  call input%get( trim(sect_name) // '.k0', this%k0 )
  call input%get( 'simulation.box.z(1)', this%z0 )

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

  integer :: i, j, m, max_mode
  real :: r, z, arr, ari, air, aii, env
  character(len=32), save :: sname = 'launch_profile_laser'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  max_mode = size(ar_im)
  do m = 0, max_mode
    do j = 1, this%nzp
      ! note that here "z" refers to xi = t - z
      z = ( this%noff_z + j - 1 ) * this%dz + this%z0
      do i = 1, this%nrp
        r = ( this%noff_r + i - 1 ) * this%dr
        call this%get_prof_perp( r, z, this%k0, this%prof_perp_pars, m, arr, ari, air, aii )
        call this%get_prof_lon( z, this%prof_lon_pars, env )
        env = env * this%a0
        ar_re(m)%f2( 1, i, j ) = env * arr
        ai_re(m)%f2( 1, i, j ) = env * air
        if ( m == 0 ) cycle
        ar_im(m)%f2( 1, i, j ) = env * ari
        ai_im(m)%f2( 1, i, j ) = env * aii
      enddo
    enddo
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine launch_profile_laser

end module profile_laser_class