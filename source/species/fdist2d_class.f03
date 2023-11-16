module fdist2d_class

use parallel_module
use options_class
use ufield_class
use input_class
use param
use sysutil_module
use fdist2d_lib
use math_module
use m_fparser
implicit none

private

public :: fdist2d

type :: fdist2d
  ! Density profiles in perpendicular and longitudinal directions
  integer, dimension(2) :: prof_type
  ! Number of particles per cell
  integer, dimension(2) :: ppc
  ! Number of azimuthal divisions
  integer :: num_theta
  ! Maximum number of particles in this partition
  integer :: npmax
  ! Frequency (number of 2D steps) of performing 2D particle sorting
  integer :: sort_freq
  real :: dr
  ! Maximum effective time step for subcycling pusher
  real :: dt_eff_max
  ! Clamped value of gamma/(1 + psi) for clamp pusher
  real :: fac_clamp
  integer :: noff, nr, nrp, max_mode
  ! Thermal velocity
  real, dimension(p_p_dim) :: uth
  ! Charge mass ratio
  real :: qm
  ! Global density
  real :: density
  ! Minimum density for particle injection
  real :: den_min
  ! Switch if setting neutralized background
  logical :: neutralized
  ! if distribute particles randomly in theta
  logical :: random_theta
  ! Parameter list of longitudinal density profile
  real, dimension(:), pointer :: prof_pars_lon => null()
  ! Parameter list of perpendicular density profile
  real, dimension(:), pointer :: prof_pars_perp => null()

  ! analytic density math function
  type(t_fparser), pointer :: math_func => null()


  procedure( set_prof_interface ), nopass, pointer :: set_prof_lon => null()
  procedure( set_prof_interface ), nopass, pointer :: set_prof_perp => null()

  procedure( set_prof_interface_analytic ), nopass, pointer :: set_prof_analytic => null()

  procedure( get_den_lon_interface ), nopass, pointer :: get_den_lon => null()
  procedure( get_den_perp_interface ), nopass, pointer :: get_den_perp => null()

  contains
  procedure :: new    => init_fdist2d
  procedure :: del    => end_fdist2d
  procedure :: inject => inject_fdist2d

end type fdist2d

interface
  subroutine get_den_perp_interface( x, s, prof_pars_perp, prof_pars_lon, den_value, math_func )
    import t_fparser
    implicit none
    real, intent(in), dimension(2) :: x
    real, intent(in) :: s
    real, intent(in), dimension(:), pointer :: prof_pars_perp, prof_pars_lon
    real, intent(out) :: den_value
    type(t_fparser), pointer, intent(inout) :: math_func
  end subroutine get_den_perp_interface
end interface

interface
  subroutine get_den_lon_interface( s, prof_pars_lon, den_value )
    implicit none
    real, intent(in) :: s
    real, intent(in), dimension(:), pointer :: prof_pars_lon
    real, intent(out) :: den_value
  end subroutine get_den_lon_interface
end interface

interface
  subroutine set_prof_interface( input, sect_name, prof_pars )
    import input_json
    implicit none
    type( input_json ), intent(inout) :: input
    character(len=*), intent(in) :: sect_name
    real, intent(inout), dimension(:), pointer :: prof_pars
  end subroutine set_prof_interface
end interface

interface
  subroutine set_prof_interface_analytic( input, sect_name, prof_pars, math_func )
    import input_json, t_fparser
    implicit none
    type( input_json ), intent(inout) :: input
    character(len=*), intent(in) :: sect_name
    real, intent(inout), dimension(:), pointer :: prof_pars
    type(t_fparser), pointer, intent(inout) :: math_func
  end subroutine set_prof_interface_analytic
end interface

character(len=20), parameter :: cls_name = "fdist2d"
integer, parameter :: cls_level = 2

contains

subroutine init_fdist2d( this, input, opts, sect, sect_id )

  implicit none

  class( fdist2d ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  character(len=*), intent(in) :: sect
  integer, intent(in) :: sect_id

  real :: xtra
  integer :: npmax_min
  character(len=20) :: sect_name
  character(len=:), allocatable :: prof_name
  character(len=18), save :: sname = 'init_fdist2d'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%noff = opts%get_noff(1)
  this%nr   = opts%get_nd(1)
  this%nrp  = opts%get_ndp(1)
  this%dr   = opts%get_dr()

  sect_name = trim(sect) // '(' // num2str(sect_id) // ')'
  call input%get( 'simulation.max_mode', this%max_mode )

  ! read profile types
  call input%get( trim(sect_name) // '.profile(1)', prof_name )
  select case ( trim(prof_name) )

    case ( 'uniform' )
      this%prof_type(1)  = p_prof_uniform
      this%set_prof_perp => set_prof_perp_uniform
      this%get_den_perp  => get_den_perp_uniform

    case ( 'parabolic-channel' )
      this%prof_type(1)  = p_prof_para_chl
      this%set_prof_perp => set_prof_perp_para_chl
      this%get_den_perp  => get_den_perp_para_chl

    case ( 'hollow-channel' )
      this%prof_type(1)  = p_prof_hllw_chl
      this%set_prof_perp => set_prof_perp_hllw_chl
      this%get_den_perp  => get_den_perp_hllw_chl

    case ( 'analytic' )
      this%prof_type(1)  = p_prof_analytic
      this%set_prof_analytic => set_prof_analytic
      this%get_den_perp  => get_den_perp_analytic

    case default
      call write_err( 'Invalid transverse density profile! Currently available &
        &include "uniform", "parabolic-channel" and "hollow-channel".' )

  end select

  call input%get( trim(sect_name) // '.profile(2)', prof_name )
  select case ( trim(prof_name) )

    case ( 'uniform' )
      this%prof_type(2) = p_prof_uniform
      this%set_prof_lon => set_prof_lon_uniform
      this%get_den_lon  => get_den_lon_uniform

    case ( 'piecewise-linear' )
      this%prof_type(2) = p_prof_pw_linear
      this%set_prof_lon => set_prof_lon_pw_linear
      this%get_den_lon  => get_den_lon_pw_linear

    case ( 'sine-series' )
      this%prof_type(2) = p_prof_sine
      this%set_prof_lon => set_prof_lon_sine
      this%get_den_lon  => get_den_lon_sine

    case ( 'cubic-spline' )
      this%prof_type(2) = p_prof_cubic_spline
      this%set_prof_lon => set_prof_lon_cubic_spline
      this%get_den_lon  => get_den_lon_cubic_spline

    case ( 'analytic' )
      this%prof_type(2)  = p_prof_analytic
      this%set_prof_analytic => set_prof_analytic
      this%get_den_lon  => get_den_lon_analytic

    case default
      call write_err( 'Invalid longitudinal density profile! Currently available &
        &include "uniform", "piecewise-linear", "sine-series" and "cubic-spline".' )

  end select


  ! read and store the profile parameters into the parameter lists
  if(this%prof_type(1) == p_prof_analytic .or. this%prof_type(2) == p_prof_analytic) then
    call this%set_prof_analytic( input, trim(sect_name), this%prof_pars_perp, this%math_func )
  else
    call this%set_prof_perp( input, trim(sect_name), this%prof_pars_perp )
    call this%set_prof_lon( input, trim(sect_name), this%prof_pars_lon )
  endif

  call input%get( trim(sect_name) // '.ppc(1)', this%ppc(1) )
  call input%get( trim(sect_name) // '.ppc(2)', this%ppc(2) )
  call input%get( trim(sect_name) // '.num_theta', this%num_theta )
  call input%get( trim(sect_name) // '.q', this%qm )
  call input%get( trim(sect_name) // '.density', this%density )

  this%den_min = 1.0d-10
  if ( input%found( trim(sect_name) // '.den_min' ) ) then
    call input%get( trim(sect_name) // '.den_min', this%den_min )
  endif

  this%dt_eff_max = 10.0
  if ( input%found( trim(sect_name) // '.dt_eff_max' ) ) then
    call input%get( trim(sect_name) // '.dt_eff_max', this%dt_eff_max )
  endif

  this%fac_clamp = 10.0
  if ( input%found( trim(sect_name) // '.fac_clamp' ) ) then
    call input%get( trim(sect_name) // '.fac_clamp', this%fac_clamp )
  endif

  this%uth = 0.0
  if ( input%found( trim(sect_name) // '.uth' ) ) then
    call input%get( trim(sect_name) // '.uth(1)', this%uth(1) )
    call input%get( trim(sect_name) // '.uth(2)', this%uth(2) )
    call input%get( trim(sect_name) // '.uth(3)', this%uth(3) )
  endif

  this%random_theta = .true.
  if ( input%found( trim(sect_name) // '.random_theta' ) ) then
    call input%get( trim(sect_name) // '.random_theta', this%random_theta )
  endif

  this%neutralized = .true.
  if ( input%found( trim(sect_name) // '.neutralized' ) ) then
    call input%get( trim(sect_name) // '.neutralized', this%neutralized )
  endif

  this%sort_freq = 0
  if ( input%found( trim(sect_name) // '.sort_freq' ) ) then
    call input%get( trim(sect_name) // '.sort_freq', this%sort_freq )
  endif

  ! calculate the maximum particles number allowed in this partition
  xtra = 2.0
  npmax_min = this%nrp * product(this%ppc) * this%num_theta
  this%npmax = int( npmax_min * xtra )
  if ( input%found( trim(sect_name) // '.npmax' ) ) then
    call input%get( trim(sect_name) // '.npmax', this%npmax )
    if ( this%npmax < npmax_min ) then
      call write_err( 'npmax is too small to initialize the 2D particles.' )
    endif
  endif

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist2d

subroutine end_fdist2d( this )

  implicit none

  class( fdist2d ), intent(inout) :: this
  character(len=18), save :: sname = 'end_fdist2d'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  if ( associated(this%prof_pars_lon) ) deallocate( this%prof_pars_lon )
  if ( associated(this%prof_pars_perp) ) deallocate( this%prof_pars_perp )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_fdist2d

subroutine inject_fdist2d( this, x, p, gamma, psi, q, npp, s )

  implicit none

  class( fdist2d ), intent(inout) :: this
  real, intent(inout), dimension(:,:) :: x, p
  real, intent(inout), dimension(:) :: gamma, psi, q
  integer(kind=LG), intent(inout) :: npp
  real, intent(in) :: s

  integer :: i, j, i1, i2, ppc_tot, n_theta, nrp, noff
  integer(kind=LG) :: ipart
  real :: dr, dtheta, den_lon, den_perp, rn, theta, coef
  real, dimension(2) :: x_tmp
  character(len=18), save :: sname = 'inject_fdist2d'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp     = this%nrp
  noff    = this%noff
  n_theta = this%num_theta
  dr      = this%dr
  dtheta  = 2.0 * pi / n_theta
  ppc_tot = product( this%ppc )
  call this%get_den_lon( s, this%prof_pars_lon, den_lon )
  coef = sign(1.0, this%qm) / ( real(ppc_tot) * real(n_theta) )

  ipart = 0
  do j = 1, n_theta
    do i = 1, nrp
      
      do i1 = 1, this%ppc(1)
        rn = (i1 - 0.5) / this%ppc(1) + real(i - 1 + noff)
        
        do i2 = 1, this%ppc(2)
          
          if (this%random_theta) then
            call random_number(theta) 
            theta = theta * 2.0 * pi
          else
            theta = ( (i2 - 0.5) / this%ppc(2) + j - 1.0 ) * dtheta
          endif

          x_tmp(1) = rn * dr * cos(theta)
          x_tmp(2) = rn * dr * sin(theta)

          call this%get_den_perp( x_tmp, s, this%prof_pars_perp, this%prof_pars_lon, den_perp, this%math_func )
          if ( den_lon * den_perp * this%density < this%den_min ) cycle

          ipart = ipart + 1
          x(1,ipart) = x_tmp(1)
          x(2,ipart) = x_tmp(2)
          q(ipart) = rn * den_perp * den_lon * this%density * coef
          p(1,ipart) = this%uth(1) * rand_norm()
          p(2,ipart) = this%uth(2) * rand_norm()
          p(3,ipart) = this%uth(3) * rand_norm()
          gamma(ipart) = sqrt( 1.0 + p(1,ipart)**2 + p(2,ipart)**2 + p(3,ipart)**2 )
          psi(ipart) = (1.0 - gamma(ipart) + p(3,ipart)) / this%qm
          
        enddo
      enddo
    enddo
  enddo

  npp = ipart

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine inject_fdist2d

end module fdist2d_class