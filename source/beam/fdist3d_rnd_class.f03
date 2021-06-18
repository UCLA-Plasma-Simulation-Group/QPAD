module fdist3d_rnd_class

use parallel_module
use options_class
use input_class
use iso_c_binding
use fdist3d_class
use fdist3d_rnd_lib
use part3d_class
use math_module
use sysutil_module

implicit none

private

type, extends(fdist3d) :: fdist3d_rnd

  ! private

  ! Density profiles in perpendicular and longitudinal directions
  integer, dimension(p_x_dim) :: prof_type

  ! Total number of macroparticles of this beam
  integer :: tot_num

  ! Initialization geometry
  integer :: geom

  ! Total charge in the unit of (e*n_p*w_p^-3*c^3)
  real :: tot_charge

  ! Thermal velocity
  real, dimension(p_p_dim) :: uth

  ! Gamma
  real :: gamma

  ! Twiss parameter alpha, beta in x and y directions
  real, dimension(2) :: alpha, beta, ctr

  ! Transverse position offset (as a function of xi)
  real, dimension(:), allocatable :: perp_offset_x, perp_offset_y

  ! Parameter list of the density profile
  real, dimension(:), pointer :: prof_pars1 => null()
  real, dimension(:), pointer :: prof_pars2 => null()
  real, dimension(:), pointer :: prof_pars3 => null()

  ! Ranges of particle injection. Particles beyond the ranges will be truncated
  real, dimension(2,p_x_dim) :: range

  ! Procedure of setting the density profiles
  procedure( set_prof_intf ), nopass, pointer :: set_prof1 => null()
  procedure( set_prof_intf ), nopass, pointer :: set_prof2 => null()
  procedure( set_prof_intf ), nopass, pointer :: set_prof3 => null()

  ! Procedure of getting the beam density
  procedure( get_rndpos_intf ), nopass, pointer :: get_rndpos1 => null()
  procedure( get_rndpos_intf ), nopass, pointer :: get_rndpos2 => null()
  procedure( get_rndpos_intf ), nopass, pointer :: get_rndpos3 => null()

  contains

  procedure :: new => init_fdist3d_rnd
  procedure :: del => end_fdist3d_rnd
  procedure :: inject => inject_fdist3d_rnd

end type fdist3d_rnd

interface
  subroutine get_rndpos_intf( prof_pars, pos )
    implicit none
    real, intent(in), dimension(:), pointer :: prof_pars
    real, intent(out) :: pos
  end subroutine get_rndpos_intf
end interface

interface
  subroutine set_prof_intf( input, sect_name, dim, prof_pars )
    import input_json
    implicit none
    type( input_json ), intent(inout) :: input
    character(len=*), intent(in) :: sect_name
    integer, intent(in) :: dim
    real, intent(inout), dimension(:), pointer :: prof_pars
  end subroutine set_prof_intf
end interface

character(len=32), save :: cls_name = 'fdist3d_rnd'
integer, save :: cls_level = 2

public :: fdist3d_rnd

contains

subroutine init_fdist3d_rnd( this, input, opts, sect_id )

  implicit none

  class( fdist3d_rnd ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  integer, intent(in) :: sect_id

  real :: xtra
  integer :: npmax_tmp
  integer(kind=LG) :: npmax_guess
  character(len=20) :: sect_name
  character(len=:), allocatable :: read_str
  character(len=18), save :: sname = 'init_fdist3d_rnd'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%noff_r = opts%get_noff(1)
  this%noff_z = opts%get_noff(2)
  this%nr     = opts%get_nd(1)
  this%nz     = opts%get_nd(2)
  this%nrp    = opts%get_ndp(1)
  this%nzp    = opts%get_ndp(2)
  this%dr     = opts%get_dr()
  this%dz     = opts%get_dxi()

  sect_name = 'beam(' // num2str(sect_id) // ')'

  ! read and set profile types
  call input%get( trim(sect_name) // '.profile(1)', read_str )
  select case ( trim(read_str) )

    case ( 'uniform' )
      this%prof_type(1) = p_prof_uniform
      this%set_prof1    => set_prof_uniform
      this%get_rndpos1  => get_rndpos_uniform

    case ( 'gaussian' )
      this%prof_type(1) = p_prof_gaussian
      this%set_prof1    => set_prof_gaussian
      this%get_rndpos1  => get_rndpos_gaussian

    case ( 'piecewise-linear' )
      this%prof_type(1) = p_prof_pw_linear
      this%set_prof1    => set_prof_pw_linear
      this%get_rndpos1  => get_rndpos_pw_linear

    case default
      call write_err( 'Invalid density profile in direction 1! Currently available &
        &include "uniform", "gaussian" and "piecewise-linear".' )

  end select

  call input%get( trim(sect_name) // '.profile(2)', read_str )
  select case ( trim(read_str) )

    case ( 'uniform' )
      this%prof_type(2) = p_prof_uniform
      this%set_prof2    => set_prof_uniform
      this%get_rndpos2  => get_rndpos_uniform

    case ( 'gaussian' )
      this%prof_type(2) = p_prof_gaussian
      this%set_prof2    => set_prof_gaussian
      this%get_rndpos2  => get_rndpos_gaussian

    case ( 'piecewise-linear' )
      this%prof_type(2) = p_prof_pw_linear
      this%set_prof2    => set_prof_pw_linear
      this%get_rndpos2  => get_rndpos_pw_linear

    case default
      call write_err( 'Invalid density profile in direction 2! Currently available &
        &include "uniform", "gaussian" and "piecewise-linear".' )

  end select

  call input%get( trim(sect_name) // '.profile(3)', read_str )
  select case ( trim(read_str) )

    case ( 'uniform' )
      this%prof_type(3) = p_prof_uniform
      this%set_prof3    => set_prof_uniform
      this%get_rndpos3  => get_rndpos_uniform

    case ( 'gaussian' )
      this%prof_type(3) = p_prof_gaussian
      this%set_prof3    => set_prof_gaussian
      this%get_rndpos3  => get_rndpos_gaussian

    case ( 'piecewise-linear' )
      this%prof_type(3) = p_prof_pw_linear
      this%set_prof3    => set_prof_pw_linear
      this%get_rndpos3  => get_rndpos_pw_linear

    case default
      call write_err( 'Invalid density profile in direction 3! Currently available &
        &include "uniform", "gaussian" and "piecewise-linear".' )

  end select

  ! read and store the profile parameters into the parameter lists
  call this%set_prof1( input, trim(sect_name), 1, this%prof_pars1 )
  call this%set_prof2( input, trim(sect_name), 2, this%prof_pars2 )
  call this%set_prof3( input, trim(sect_name), 3, this%prof_pars3 )

  call input%get( 'simulation.box.z(1)', this%z0 )
  call input%get( trim(sect_name) // '.total_num', this%tot_num )
  call input%get( trim(sect_name) // '.total_charge', this%tot_charge )
  call input%get( trim(sect_name) // '.q', this%qm )
  call input%get( trim(sect_name) // '.m', this%qbm )
  this%qbm = this%qm / this%qbm
  call input%get( trim(sect_name) // '.gamma', this%gamma )
  call input%get( trim(sect_name) // '.range1(1)', this%range(p_lower, 1) )
  call input%get( trim(sect_name) // '.range1(2)', this%range(p_upper, 1) )
  call input%get( trim(sect_name) // '.range2(1)', this%range(p_lower, 2) )
  call input%get( trim(sect_name) // '.range2(2)', this%range(p_upper, 2) )
  call input%get( trim(sect_name) // '.range3(1)', this%range(p_lower, 3) )
  call input%get( trim(sect_name) // '.range3(2)', this%range(p_upper, 3) )
  this%range(:, 3) = this%range(:, 3) - this%z0

  this%geom = p_geom_cart
  if ( input%found( trim(sect_name) // '.geometry' ) ) then
    call input%get( trim(sect_name) // '.geometry', read_str )
    select case ( trim(read_str) )
    case ( 'cartesian' )
      this%geom = p_geom_cart
    case ( 'cylindrical' )
      this%geom = p_geom_cyl
    end select
  endif

  this%evol = .true.
  if ( input%found( trim(sect_name) // '.evolution' ) ) then
    call input%get( trim(sect_name) // '.evolution', this%evol )
  endif

  this%quiet = .false.
  if ( input%found( trim(sect_name) // '.quiet_start' ) ) then
    call input%get( trim(sect_name) // '.quiet_start', this%quiet )
  endif

  this%has_spin = .false.
  if ( input%found( trim(sect_name) // '.has_spin' ) ) then
    call input%get( trim(sect_name) // '.has_spin', this%has_spin )
    if ( this%has_spin ) then
      call input%get( trim(sect_name) // '.anom_mag_moment', this%amm )
    endif
  endif

  this%uth = 0.0
  if ( input%found( trim(sect_name) // '.uth' ) ) then
    call input%get( trim(sect_name) // '.uth(1)', this%uth(1) )
    call input%get( trim(sect_name) // '.uth(2)', this%uth(2) )
    call input%get( trim(sect_name) // '.uth(3)', this%uth(3) )
  endif

  this%alpha = 0.0
  if ( input%found( trim(sect_name) // '.alpha' ) ) then
    if ( this%prof_type(1) == p_prof_gaussian .and. this%geom == p_geom_cart ) then
      call input%get( trim(sect_name) // '.alpha(1)', this%alpha(1) )
      call input%get( trim(sect_name) // '.gauss_center(1)', this%ctr(1) )
      call input%get( trim(sect_name) // '.gauss_sigma(1)', this%beta(1) )
      this%beta(1) = this%gamma * this%beta(1) / this%uth(1)
    else
      call write_err( 'Twiss parameter alpha is only available for Gaussian &
        &profile in the Cartesian geometry.' )
    endif
    if ( this%prof_type(2) == p_prof_gaussian .and. this%geom == p_geom_cart ) then
      call input%get( trim(sect_name) // '.alpha(2)', this%alpha(2) )
      call input%get( trim(sect_name) // '.gauss_center(2)', this%ctr(2) )
      call input%get( trim(sect_name) // '.gauss_sigma(2)', this%beta(2) )
      this%beta(2) = this%gamma * this%beta(2) / this%uth(2)
    else
      call write_err( 'Twiss parameter alpha is only available for Gaussian &
        &profile in the Cartesian geometry.' )
    endif
  endif

  if ( input%found( trim(sect_name) // '.perp_offset_x' ) ) then
    call input%get( trim(sect_name) // '.perp_offset_x', this%perp_offset_x )
    this%perp_offset_x(1) = this%perp_offset_x(1) - this%z0
  endif
  if ( input%found( trim(sect_name) // '.perp_offset_y' ) ) then
    call input%get( trim(sect_name) // '.perp_offset_y', this%perp_offset_y )
    this%perp_offset_y(1) = this%perp_offset_y(1) - this%z0
  endif

  ! if npmax is not given, guess a value for npmax
  xtra = 1.5
  ! this guess value must be smaller than the practical value, so the buffer reallocation
  ! will be invoked. Need a more smarter estimation algorithm.
  npmax_guess = this%tot_num / num_procs()
  this%npmax = int( npmax_guess * xtra, kind=LG )

  ! if npmax is given, set it as the maximum of the given value and p_npmax_min
  if ( input%found( trim(sect_name) // '.npmax' ) ) then
    call input%get( trim(sect_name) // '.npmax', npmax_tmp )
    this%npmax = int( npmax_tmp, kind=LG )
    if ( this%npmax < p_npmax_min ) then
      call write_wrn( 'The npmax may be too small for the initialization. It has &
        &been automatically changed to ' // num2str(p_npmax_min) // '.' )
      this%npmax = int( p_npmax_min, kind=LG )
    endif
  endif

  ! if using quiet start, double the buffer
  if ( this%quiet ) then
    this%npmax = this%npmax * 2
    call write_stdout( 'The number of particles will be doubled for quiet-start &
      &initialization.' )
  endif

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist3d_rnd

subroutine end_fdist3d_rnd( this )

  implicit none

  class( fdist3d_rnd ), intent(inout) :: this
  character(len=18), save :: sname = 'end_fdist3d_rnd'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  if ( associated(this%prof_pars1) ) deallocate( this%prof_pars1 )
  if ( associated(this%prof_pars2) ) deallocate( this%prof_pars2 )
  if ( associated(this%prof_pars3) ) deallocate( this%prof_pars3 )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_fdist3d_rnd

subroutine inject_fdist3d_rnd( this, part )

  implicit none

  class( fdist3d_rnd ), intent(inout) :: this
  class( part3d ), intent(inout) :: part

  integer :: i
  integer(kind=LG) :: ip
  real :: q_per_part, r_tmp, ratio, offset, coef
  real, dimension(3) :: x_tmp
  real, dimension(2) :: edge_r, edge_z
  character(len=18), save :: sname = 'inject_fdist3d_rnd'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  edge_r(p_lower) = this%noff_r * this%dr
  edge_z(p_lower) = this%noff_z * this%dz
  edge_r(p_upper) = edge_r(p_lower) + this%nrp * this%dr
  edge_z(p_upper) = edge_z(p_lower) + this%nzp * this%dz

  q_per_part = sign(1.0, this%qm) * this%tot_charge / ( 2.0 * pi * this%tot_num * this%dr**2 * this%dz )
  if ( this%geom == p_geom_cyl ) then
    q_per_part = q_per_part * 0.5 * pi**2
  endif

  do i = 1, this%tot_num

    ! generate particles until they are located within the range.
    ! Note that if a particle is generated beyond the physical boundary, it will
    ! not be handled here but dropped later by the particle manager.
    do

      ! generate transverse position and check if it is inside the initialization range.
      call this%get_rndpos1( this%prof_pars1, x_tmp(1) )
      if ( x_tmp(1) < this%range(p_lower,1) .or. x_tmp(1) >= this%range(p_upper,1) ) cycle
      call this%get_rndpos2( this%prof_pars2, x_tmp(2) )
      if ( x_tmp(2) < this%range(p_lower,2) .or. x_tmp(2) >= this%range(p_upper,2) ) cycle

      ! check if it is inside this partition
      select case ( this%geom )
      case ( p_geom_cart )
        r_tmp = sqrt( x_tmp(1) * x_tmp(1) + x_tmp(2) * x_tmp(2) )
      case ( p_geom_cyl )
        r_tmp = x_tmp(1)
      end select
      if ( r_tmp < edge_r(p_lower) .or. r_tmp >= edge_r(p_upper) ) cycle

      ! generate longitudinal position and check if it is inside the initialization range.
      call this%get_rndpos3( this%prof_pars3, x_tmp(3) )
      if ( x_tmp(3) < this%range(p_lower,3) .or. x_tmp(3) >= this%range(p_upper,3) ) cycle

      ! check if it is inside this partition
      if ( x_tmp(3) < edge_z(p_lower) .or. x_tmp(3) >= edge_z(p_upper) ) cycle

      exit

    enddo

    ! check if needs to reallocate the particle buffer
    if ( part%npp >= part%npmax ) then
      call part%realloc( ratio = p_buf_incr, buf_type = 'particle' )
    endif

    part%npp = part%npp + 1
    ip = part%npp

    select case ( this%geom )
    case ( p_geom_cart )
      part%x(1, ip) = x_tmp(1)
      part%x(2, ip) = x_tmp(2)
      part%x(3, ip) = x_tmp(3)
      part%q(ip) = q_per_part
    case ( p_geom_cyl )
      part%x(1, ip) = x_tmp(1) * cos(x_tmp(2))
      part%x(2, ip) = x_tmp(1) * sin(x_tmp(2))
      part%x(3, ip) = x_tmp(3)
      part%q(ip) = q_per_part * r_tmp
    end select

    ! momentum initialization uses Cartesian geometry
    part%p(1, ip) = this%uth(1) * ranorm()
    part%p(2, ip) = this%uth(2) * ranorm()
    part%p(3, ip) = this%uth(3) * ranorm() + this%gamma
    part%p(3, ip) = sqrt( part%p(3, ip)**2 - part%p(1, ip)**2 - part%p(2, ip)**2 - 1 )

    ! spin initialization has not yet implemented
    if ( this%has_spin ) then
      part%s(1, ip) = 0.0
      part%s(2, ip) = 0.0
      part%s(3, ip) = 1.0
    endif

  enddo

  ! if using quiet start. Note that the quiet start will double the total particle number
  if ( this%quiet ) then

    if ( part%npp * 2 > part%npmax ) then
      ratio = real(part%npp * 2) / part%npmax
      call part%realloc( ratio = p_buf_incr * ratio, buf_type = 'particle' )
    endif

    do i = 1, ip

      part%x(1, ip+i) = -part%x(1,i)
      part%x(2, ip+i) = -part%x(2,i)
      part%x(3, ip+i) =  part%x(3,i)

      part%p(1, ip+i) = -part%p(1,i)
      part%p(2, ip+i) = -part%p(2,i)
      part%p(3, ip+i) =  part%p(3,i)

      part%q(i) = 0.5 * part%q(i)
      part%q(ip+i) = part%q(i)

      ! spin initialization has not yet implemented
      if ( this%has_spin ) then
        part%s(1,ip+i) = 0.0
        part%s(2,ip+i) = 0.0
        part%s(3,ip+i) = 1.0
      endif

    enddo
    part%npp = part%npp * 2

  endif

  ! use Twiss parameter to initialize the tilt phase-space ellipse
  if ( this%alpha(1) /= 0.0 ) then
    coef = this%gamma * this%alpha(1) / this%beta(1)
    do i = 1, part%npp
      part%p(1,i) = part%p(1,i) - coef * ( part%x(1,i) - this%ctr(1) )
    enddo
  endif
  if ( this%alpha(2) /= 0.0 ) then
    coef = this%gamma * this%alpha(2) / this%beta(2)
    do i = 1, part%npp
      part%p(2,i) = part%p(2,i) - coef * ( part%x(2,i) - this%ctr(2) )
    enddo
  endif

  ! offset the beam particles according to xi
  if ( allocated(this%perp_offset_x) ) then
    do i = 1, part%npp
      offset = eval_polynomial( part%x(3,i), this%perp_offset_x(1), this%perp_offset_x(2:) )
      part%x(1,i) =  part%x(1,i) + offset
    enddo
  endif
  if ( allocated(this%perp_offset_y) ) then
    do i = 1, part%npp
      offset = eval_polynomial( part%x(3,i), this%perp_offset_y(1), this%perp_offset_y(2:) )
      part%x(2,i) =  part%x(2,i) + offset
    enddo
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine inject_fdist3d_rnd

end module fdist3d_rnd_class