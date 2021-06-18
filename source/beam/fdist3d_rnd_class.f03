module fdist3d_rnd_class

use parallel_module
use options_class
use input_class
use iso_c_binding
use fdist3d_class
use fdist3d_rnd_lib
use random
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
  integer(kind=LG) :: npmax_min
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

  ! TODO: buffer reallocation needs to be implemented here
  ! calculate the maximum particles number allowed in this partition
  xtra = 2.0
  npmax_min = this%tot_num
  if ( this%quiet ) then
    npmax_min = npmax_min * 2
    call write_stdout( 'The number of particles will be doubled for quiet-start &
      &initialization.' )
  endif
  this%npmax = int( npmax_min * xtra, kind=LG )
  if ( input%found( trim(sect_name) // '.npmax' ) ) then
    call input%get( trim(sect_name) // '.npmax', npmax_tmp )
    this%npmax = int( npmax_tmp, kind=LG )
    ! if ( this%npmax < npmax_min ) then
    !   call write_err( 'npmax is too small to initialize the 3D particles.' )
    ! endif
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

subroutine inject_fdist3d_rnd( this, x, p, s, q, npp )

  implicit none

  class( fdist3d_rnd ), intent(inout) :: this
  real, intent(inout), dimension(:,:) :: x, p, s
  real, intent(inout), dimension(:) :: q
  integer(kind=LG), intent(inout) :: npp

  integer :: i
  integer(kind=LG) :: ipart
  real :: q_per_part, r_tmp
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

  ipart = 0
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

    ipart = ipart + 1

    select case ( this%geom )
    case ( p_geom_cart )
      x(1, ipart) = x_tmp(1)
      x(2, ipart) = x_tmp(2)
      x(3, ipart) = x_tmp(3)
      q(ipart) = q_per_part
    case ( p_geom_cyl )
      x(1, ipart) = x_tmp(1) * cos(x_tmp(2))
      x(2, ipart) = x_tmp(1) * sin(x_tmp(2))
      x(3, ipart) = x_tmp(3)
      q(ipart) = q_per_part * r_tmp
    end select

    ! momentum initialization uses Cartesian geometry
    p(1, ipart) = this%uth(1) * ranorm()
    p(2, ipart) = this%uth(2) * ranorm()
    p(3, ipart) = this%uth(3) * ranorm() + this%gamma
    p(3, ipart) = sqrt( p(3, ipart)**2 - p(1, ipart)**2 - p(2, ipart)**2 - 1 )

    ! spin initialization has not yet implemented
    if ( this%has_spin ) then
      s(1, ipart) = 0.0
      s(2, ipart) = 0.0
      s(3, ipart) = 1.0
    endif

  enddo

  ! if using quiet start. Note that the quiet start will double the total particle number
  if ( this%quiet ) then

    do i = 1, ipart

      x(1, ipart+i) = -x(1,i)
      x(2, ipart+i) = -x(2,i)
      x(3, ipart+i) =  x(3,i)

      p(1, ipart+i) = -p(1,i)
      p(2, ipart+i) = -p(2,i)
      p(3, ipart+i) =  p(3,i)

      q(i) = 0.5 * q(i)
      q(ipart+i) = q(i)

      ! spin initialization has not yet implemented
      if ( this%has_spin ) then
        s(1,ipart+i) = 0.0
        s(2,ipart+i) = 0.0
        s(3,ipart+i) = 1.0
      endif

    enddo
    ipart = ipart * 2

  endif

  npp = ipart

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine inject_fdist3d_rnd

end module fdist3d_rnd_class