module field_src_class

use field_class
use parallel_pipe_class
use grid_class
use ufield_class
use param
use system

implicit none

private

public :: field_rho, field_jay, field_djdxi

type, extends( field ) :: field_rho

  contains

  generic :: new => init_field_rho
  procedure, private :: init_field_rho

end type field_rho

type, extends( field ) :: field_jay

  contains

  generic :: new => init_field_jay
  procedure, private :: init_field_jay

end type field_jay

type, extends( field ) :: field_djdxi

  contains

  generic :: new => init_field_djdxi
  generic :: solve => solve_field_djdxi
  procedure, private :: init_field_djdxi, solve_field_djdxi

end type field_djdxi

contains

subroutine init_field_rho( this, pp, gp, dr, dxi, num_modes, part_shape, &
  smooth_type, smooth_order )

  implicit none

  class( field_rho ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape
  real, intent(in) :: dr, dxi
  integer, intent(in), optional :: smooth_type, smooth_order

  integer, dimension(2,2) :: gc_num
  integer :: dim, i
  character(len=20), save :: cls_name = "field_rho"
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = "init_field_rho"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/1, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 1
  ! call initialization routine of the parent class
  if ( present(smooth_type) .and. present(smooth_order) ) then

    gc_num(1,1) = max( gc_num(1,1), smooth_order )
    gc_num(2,1) = max( gc_num(2,1), smooth_order )

    call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num, &
      smooth_type=smooth_type, smooth_order=smooth_order )

  else

    call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_rho

subroutine init_field_jay( this, pp, gp, dr, dxi, num_modes, part_shape, &
  smooth_type, smooth_order )

  implicit none

  class( field_jay ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape
  real, intent(in) :: dr, dxi
  integer, intent(in), optional :: smooth_type, smooth_order

  integer, dimension(2,2) :: gc_num
  integer :: dim, i
  character(len=20), save :: cls_name = "field_jay"
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = "init_field_jay"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/1, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 3
  ! call initialization routine of the parent class
  if ( present(smooth_type) .and. present(smooth_order) ) then

    gc_num(1,1) = max( gc_num(1,1), smooth_order )
    gc_num(2,1) = max( gc_num(2,1), smooth_order )

    call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num, &
      smooth_type=smooth_type, smooth_order=smooth_order )

  else

    call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_jay

subroutine init_field_djdxi( this, pp, gp, dr, dxi, num_modes, part_shape, &
  smooth_type, smooth_order )

  implicit none

  class( field_djdxi ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape
  real, intent(in) :: dr, dxi
  integer, intent(in), optional :: smooth_type, smooth_order

  integer, dimension(2,2) :: gc_num
  integer :: dim, i
  character(len=20), save :: cls_name = "field_djdxi"
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = "init_field_djdxi"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/1, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 2
  ! call initialization routine of the parent class
  if ( present(smooth_type) .and. present(smooth_order) ) then

    gc_num(1,1) = max( gc_num(1,1), smooth_order )
    gc_num(2,1) = max( gc_num(2,1), smooth_order )

    call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num, &
      smooth_type=smooth_type, smooth_order=smooth_order )

  else

    call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_djdxi

subroutine solve_field_djdxi( this, acu, amu )

  implicit none

  class( field_djdxi ), intent(inout) :: this
  class( field_djdxi ), intent(in) :: acu
  class( field_jay ), intent(inout) :: amu

  type( ufield ), dimension(:), pointer :: acu_re => null(), acu_im => null()
  type( ufield ), dimension(:), pointer :: amu_re => null(), amu_im => null()
  real, dimension(:,:), pointer :: uacu_re => null(), uacu_im => null()
  real, dimension(:,:), pointer :: uamu_re => null(), uamu_im => null()
  real, dimension(:,:), pointer :: udcu_re => null(), udcu_im => null()
  integer :: mode, i, nrp, noff, idproc, nvp
  real :: idr, idrh, ir, k0
  character(len=20), save :: cls_name = "field_djdxi"  
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = 'solve_field_djdxi'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'set source' )

  idr = 1.0 / this%dr
  idrh = idr * 0.5
  nrp = this%rf_re(0)%get_ndp(1)

  acu_re => acu%get_rf_re()
  acu_im => acu%get_rf_im()
  amu_re => amu%get_rf_re()
  amu_im => amu%get_rf_im()

  noff = this%rf_re(0)%get_noff(1)
  nvp = this%rf_re(0)%pp%getlnvp()
  idproc = this%rf_re(0)%pp%getlidproc()

  do mode = 0, this%num_modes

    uacu_re => acu_re(mode)%get_f1()
    uamu_re => amu_re(mode)%get_f1()
    udcu_re => this%rf_re(mode)%get_f1()
    if ( mode == 0 ) then
      do i = 1, nrp
        k0 = real(i+noff) - 0.5
        ir = idr / k0
        udcu_re(1,i) = uacu_re(1,i) - idrh * ( uamu_re(1,i+1) - uamu_re(1,i-1) ) - ir*uamu_re(1,i)
        udcu_re(2,i) = uacu_re(2,i) - idrh * ( uamu_re(2,i+1) - uamu_re(2,i-1) ) - ir*uamu_re(2,i)
      enddo
      if ( idproc == 0 ) then
        udcu_re(1,1) = uacu_re(1,1) - idrh * ( 4.0 * uamu_re(1,2) - &
        &uamu_re(1,3) - 3.0 * uamu_re(1,1) ) - 2.0*idr * uamu_re(1,1)
        udcu_re(2,1) = uacu_re(2,1) - idrh * ( 4.0 * uamu_re(2,2) - &
        &uamu_re(2,3) - 3.0 * uamu_re(2,1) ) - 2.0*idr * uamu_re(2,1)
      endif
      if ( idproc == nvp-1 ) then
        k0 = real(nrp+noff)-0.5
        ir = idr / k0
        udcu_re(1,nrp) = uacu_re(1,nrp) + idrh * ( 4.0 * uamu_re(1,nrp-1) - &
        &uamu_re(1,nrp-2) - 3.0 * uamu_re(1,nrp) ) - ir * uamu_re(1,nrp)
        udcu_re(2,nrp) = uacu_re(2,nrp) + idrh * ( 4.0 * uamu_re(2,nrp-1) - &
        &uamu_re(2,nrp-2) - 3.0 * uamu_re(2,nrp) ) - ir * uamu_re(2,nrp)
      endif
      cycle
    endif

    uacu_im => acu_im(mode)%get_f1()
    uamu_im => amu_im(mode)%get_f1()
    udcu_im => this%rf_im(mode)%get_f1()

    do i = 1, nrp
      k0 = real(i+noff) - 0.5
      ir = idr / k0
      udcu_re(1,i) = uacu_re(1,i) - idrh * ( uamu_re(1,i+1) - uamu_re(1,i-1) ) - ir*uamu_re(1,i)&
      & + mode*ir*uamu_im(2,i)
      udcu_re(2,i) = uacu_re(2,i) - idrh * ( uamu_re(2,i+1) - uamu_re(2,i-1) ) - ir*uamu_re(2,i)&
      & + mode*ir*uamu_im(3,i)
      udcu_im(1,i) = uacu_im(1,i) - idrh * ( uamu_im(1,i+1) - uamu_im(1,i-1) ) - ir*uamu_im(1,i)&
      & - mode*ir*uamu_re(2,i)
      udcu_im(2,i) = uacu_im(2,i) - idrh * ( uamu_im(2,i+1) - uamu_im(2,i-1) ) - ir*uamu_im(2,i)&
      & - mode*ir*uamu_re(3,i)
    enddo
    if ( idproc == 0 ) then
      ! ir = 2.0 * idr
      udcu_re(1,1) = uacu_re(1,1) - idrh * ( 4.0 * uamu_re(1,2) - &
      &uamu_re(1,3) - 3.0 * uamu_re(1,1) ) - 2.0*idr * uamu_re(1,1) + mode*2.0*idr*uamu_im(2,1)
      udcu_re(2,1) = uacu_re(2,1) - idrh * ( 4.0 * uamu_re(2,2) - &
      &uamu_re(2,3) - 3.0 * uamu_re(2,1) ) - 2.0*idr * uamu_re(2,1) + mode*2.0*idr*uamu_im(3,1)
      udcu_im(1,1) = uacu_im(1,1) - idrh * ( 4.0 * uamu_im(1,2) - &
      &uamu_im(1,3) - 3.0 * uamu_im(1,1) ) - 2.0*idr * uamu_im(1,1) - mode*2.0*idr*uamu_re(2,1)
      udcu_im(2,1) = uacu_im(2,1) - idrh * ( 4.0 * uamu_im(2,2) - &
      &uamu_im(2,3) - 3.0 * uamu_im(2,1) ) - 2.0*idr * uamu_im(2,1) - mode*2.0*idr*uamu_re(3,1)
    endif
    if ( idproc == nvp-1 ) then
      k0 = real(nrp+noff)-0.5
      ir = idr / k0
      udcu_re(1,nrp) = uacu_re(1,nrp) + idrh * ( 4.0 * uamu_re(1,nrp-1) - &
      &uamu_re(1,nrp-2) - 3.0 * uamu_re(1,nrp) ) - ir * uamu_re(1,nrp) + mode*ir*uamu_im(2,nrp)
      udcu_re(2,nrp) = uacu_re(2,nrp) + idrh * ( 4.0 * uamu_re(2,nrp-1) - &
      &uamu_re(2,nrp-2) - 3.0 * uamu_re(2,nrp) ) - ir * uamu_re(2,nrp) + mode*ir*uamu_im(3,nrp)
      udcu_im(1,nrp) = uacu_im(1,nrp) + idrh * ( 4.0 * uamu_im(1,nrp-1) - &
      &uamu_im(1,nrp-2) - 3.0 * uamu_im(1,nrp) ) - ir * uamu_im(1,nrp) - mode*ir*uamu_re(2,nrp)
      udcu_im(2,nrp) = uacu_im(2,nrp) + idrh * ( 4.0 * uamu_im(2,nrp-1) - &
      &uamu_im(2,nrp-2) - 3.0 * uamu_im(2,nrp) ) - ir * uamu_im(2,nrp) - mode*ir*uamu_re(3,nrp)
    endif

  enddo

  call this%copy_gc_f1()

  call stop_tprof( 'set source' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_djdxi

end module field_src_class