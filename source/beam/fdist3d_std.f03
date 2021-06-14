module fdist3d_std_class

use parallel_module
use options_class
use input_class
use iso_c_binding
use part3d_class
use fdist3d_class
use fdist3d_std_lib
use random

implicit none

private

type, extends(fdist3d) :: fdist3d_std

  private

  ! Density profiles in perpendicular and longitudinal directions
  integer, dimension(p_x_dim) :: prof_type

  ! Number of particles per (virtual) cell
  integer, dimension(3) :: ppc

  ! Number of particles in the azimuthal direction
  integer :: num_theta

  ! Global density
  real :: density

  ! Minimum density for particle injection
  real :: den_min

  ! Thermal velocity
  real, dimension(p_p_dim) :: uth

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
  procedure( get_den_intf ), nopass, pointer :: get_den1 => null()
  procedure( get_den_intf ), nopass, pointer :: get_den2 => null()
  procedure( get_den_intf ), nopass, pointer :: get_den3 => null()

  contains

  procedure :: init_fdist3d   => init_fdist3d_std
  procedure :: end_fdist3d    => end_fdist3d_std
  procedure :: inject_fdist3d => inject_fdist3d_std

end type fdist3d_std

interface
  subroutine get_den_intf( x, prof_pars, den_value )
    implicit none
    real, intent(in) :: x
    real, intent(in), dimension(:), pointer :: prof_pars
    real, intent(out) :: den_value
  end subroutine get_den_lon_intf
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

character(len=32), save :: cls_name = 'fdist3d_std'
integer, save :: cls_level = 2

public :: fdist3d_std

contains

subroutine init_fdist3d_std( this, input, opts, sect_id )

  implicit none

  class( fdist3d_std ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  integer, intent(in) :: sect_id

  real :: xtra
  integer :: npmax_min
  character(len=20) :: sect_name
  character(len=:), allocatable :: prof_name
  character(len=18), save :: sname = 'init_fdist3d_std'

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
  ! call input%get( 'simulation.max_mode', this%max_mode )

  ! read and set profile types
  call input%get( trim(sect_name) // '.profile(1)', prof_name )
  select case ( trim(prof_name) )

    case ( 'uniform' )
      this%prof_type(1) = p_prof_uniform
      this%set_prof1    => set_prof_uniform
      this%get_den1     => get_den_uniform

    case ( 'gaussian' )
      this%prof_type(1) = p_prof_gaussian
      this%set_prof1    => set_prof_gaussian
      this%get_den1     => get_den_gaussian

    case default
      call write_err( 'Invalid density profile in direction 1! Currently available &
        &include "gaussian".' )

  end select

  call input%get( trim(sect_name) // '.profile(2)', prof_name )
  select case ( trim(prof_name) )

    case ( 'uniform' )
      this%prof_type(2) = p_prof_uniform
      this%set_prof2    => set_prof_uniform
      this%get_den2     => get_den_uniform

    case ( 'gaussian' )
      this%prof_type(2) = p_prof_gaussian
      this%set_prof2    => set_prof_gaussian
      this%get_den2     => get_den_gaussian

    case default
      call write_err( 'Invalid density profile in direction 2! Currently available &
        &include "gaussian".' )

  end select

  call input%get( trim(sect_name) // '.profile(3)', prof_name )
  select case ( trim(prof_name) )

    case ( 'uniform' )
      this%prof_type(3) = p_prof_uniform
      this%set_prof3    => set_prof_uniform
      this%get_den3     => get_den_uniform

    case ( 'gaussian' )
      this%prof_type(3) = p_prof_gaussian
      this%set_prof3    => set_prof_gaussian
      this%get_den3     => get_den_gaussian

    case default
      call write_err( 'Invalid density profile in direction 3! Currently available &
        &include "gaussian".' )

  end select

  ! read and store the profile parameters into the parameter lists
  call this%set_prof1( input, trim(sect_name), 1, this%prof_pars1 )
  call this%set_prof2( input, trim(sect_name), 2, this%prof_pars2 )
  call this%set_prof3( input, trim(sect_name), 3, this%prof_pars3 )

  call input%get( 'simulation.box.z(1)', this%z0 )
  call input%get( trim(sect_name) // '.ppc(1)', this%ppc(1) )
  call input%get( trim(sect_name) // '.ppc(2)', this%ppc(2) )
  call input%get( trim(sect_name) // '.ppc(3)', this%ppc(3) )
  call input%get( trim(sect_name) // '.num_theta', this%num_theta )
  call input%get( trim(sect_name) // '.q', this%qm )
  call input%get( trim(sect_name) // '.m', this%qbm )
  this%qbm = this%qm / this%qbm
  call input%get( trim(sect_name) // '.density', this%density )
  call input%get( trim(sect_name) // '.range1(1)', this%range(p_lower, 1) )
  call input%get( trim(sect_name) // '.range1(2)', this%range(p_upper, 1) )
  call input%get( trim(sect_name) // '.range2(1)', this%range(p_lower, 2) )
  call input%get( trim(sect_name) // '.range2(2)', this%range(p_upper, 2) )
  call input%get( trim(sect_name) // '.range3(1)', this%range(p_lower, 3) )
  call input%get( trim(sect_name) // '.range3(2)', this%range(p_upper, 3) )
  this%range(:, 3) = this%range(:, 3) - this%z0

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
    call input%get( trim(sect_name) // '.anom_mag_moment', this%amm )
  endif

  this%den_min = 1.0d-10
  if ( input%found( trim(sect_name) // '.den_min' ) ) then
    call input%get( trim(sect_name) // '.den_min', this%den_min )
  endif

  this%uth = 0.0
  if ( input%found( trim(sect_name) // '.uth' ) ) then
    call input%get( trim(sect_name) // '.uth(1)', this%uth(1) )
    call input%get( trim(sect_name) // '.uth(2)', this%uth(2) )
    call input%get( trim(sect_name) // '.uth(3)', this%uth(3) )
  endif

  ! calculate the maximum particles number allowed in this partition
  xtra = 2.0
  npmax_min = this%nrp * this%nzp * product(this%ppc) * this%num_theta
  if ( this%quiet ) then
    npmax_min = npmax_min * 2
    call write_stdout( 'The number of particles will be doubled for quiet-start &
      &initialization in "standard" beam profile type.' )
  endif
  this%npmax = int( npmax_min * xtra )
  if ( input%found( trim(sect_name) // '.npmax' ) ) then
    call input%get( trim(sect_name) // '.npmax', this%npmax )
    if ( this%npmax < npmax_min ) then
      call write_err( 'npmax is too small to initialize the 3D particles.' )
    endif
  endif

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist3d_std

subroutine end_fdist3d_std( this )

  implicit none

  class( fdist3d_std ), intent(inout) :: this
  character(len=18), save :: sname = 'end_fdist3d_std'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  if ( associated(this%prof_pars1) ) deallocate( this%prof_pars1 )
  if ( associated(this%prof_pars2) ) deallocate( this%prof_pars2 )
  if ( associated(this%prof_pars3) ) deallocate( this%prof_pars3 )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_fdist3d_std

subroutine inject_fdist3d_std( this, x, p, s, q, npp )

  implicit none

  class( fdist3d_std ), intent(inout) :: this
  real, intent(inout), dimension(:,:) :: x, p, s
  real, intent(inout), dimension(:) :: q
  integer(kind=LG), intent(inout) :: npp

  integer :: i, j, k, i1, i2, i3, ppc_tot, n_theta, nrp, nzp, noff_r, noff_z
  integer(kind=LG) :: ipart
  real :: dr, dz, dtheta, rn, zn, theta, coef, den_loc
  real, dimension(3) :: den_val
  real, dimension(3) :: x_tmp
  character(len=18), save :: sname = 'inject_fdist3d_std'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp     = this%nrp
  nzp     = this%nzp
  noff_r  = this%noff_r
  noff_z  = this%noff_z
  n_theta = this%num_theta
  dr      = this%dr
  dz      = this%dz
  dtheta  = 2.0 * pi / n_theta
  ppc_tot = product( this%ppc )

  call this%get_den_lon( s, this%prof_pars_lon, den_lon )
  coef = sign(1.0, this%qm) / ( real(ppc_tot) * real(n_theta) )

  ipart = 0
  do k = 1, nzp
    do j = 1, n_theta
      do i = 1, nrp
        
        do i3 = 1, this%ppc(3)
          zn = (i3 - 0.5) / this%ppc(3) + real(k - 1 + noff_z)
          do i2 = 1, this%ppc(2)
            theta = ( (i2 - 0.5) / this%ppc(2) + j - 1.5 ) * dtheta
            do i1 = 1, this%ppc(1)
              rn = (i1 - 0.5) / this%ppc(1) + real(i - 1 + noff_r)

              x_tmp(1) = rn * dr * cos(theta)
              x_tmp(2) = rn * dr * sin(theta)
              x_tmp(3) = zn * dz

              ! check if the particle to be injected is within the range
              if ( x_tmp(1) < this%range(p_lower,1) .or. x_tmp(1) > this%range(p_upper,1) ) cycle
              if ( x_tmp(2) < this%range(p_lower,2) .or. x_tmp(2) > this%range(p_upper,2) ) cycle
              if ( x_tmp(3) < this%range(p_lower,3) .or. x_tmp(3) > this%range(p_upper,3) ) cycle

              call this%get_den1( x_tmp(1), this%prof_pars1, den_val(1) )
              call this%get_den2( x_tmp(2), this%prof_pars2, den_val(2) )
              call this%get_den3( x_tmp(3), this%prof_pars3, den_val(3) )
              den_loc = product(den_val) * this%density

              if ( den_loc < this%den_min ) cycle

              ipart = ipart + 1

              x(1,ipart) = x_tmp(1)
              x(2,ipart) = x_tmp(2)
              x(3,ipart) = x_tmp(3)

              p(1,ipart) = this%uth(1) * ranorm()
              p(2,ipart) = this%uth(2) * ranorm()
              p(3,ipart) = this%uth(3) * ranorm()

              q(ipart) = rn * den_loc * coef

              ! spin initialization has not yet implemented
              if ( this%has_spin ) then
                s(1,ipart) = 0.0
                s(2,ipart) = 0.0
                s(3,ipart) = 1.0
              endif
              
            enddo  
          enddo
        enddo

      enddo
    enddo
  enddo

  if ( this%quiet ) then

    do i = 1, ipart

      x(1,ipart+i) = -x(1,i)
      x(2,ipart+i) = -x(2,i)
      x(3,ipart+i) =  x(3,i)

      p(1,ipart+i) = -p(1,i)
      p(2,ipart+i) = -p(2,i)
      p(3,ipart+i) =  p(3,i)

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

end subroutine inject_fdist3d_std

end module fdist3d_std_class