module field_src_class

use field_class
use parallel_module
use options_class
use ufield_class
use param
use sysutil_module
use mpi

implicit none

private

public :: field_rho, field_jay, field_djdxi

type, extends( field ) :: field_rho

  real :: q_ax ! on-axis charge
  integer :: smooth_order = 0

  contains

  generic :: new => init_field_rho
  procedure :: init_field_rho
  procedure :: smooth => smooth_field_rho

end type field_rho

type, extends( field ) :: field_jay

  integer :: smooth_order = 0

  contains

  generic :: new => init_field_jay
  procedure :: init_field_jay
  procedure :: smooth => smooth_field_jay

end type field_jay

type, extends( field ) :: field_djdxi

  integer :: smooth_order = 0

  contains

  generic :: new => init_field_djdxi
  procedure :: init_field_djdxi
  procedure :: solve => solve_field_djdxi
  procedure :: smooth => smooth_field_djdxi

end type field_djdxi

contains

subroutine init_field_rho(this, opts, max_mode, part_shape, smooth_order, has_2d)

  implicit none

  class( field_rho ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: max_mode, part_shape
  integer, intent(in) :: smooth_order
  logical, intent(in), optional :: has_2d

  integer, dimension(2,2) :: gc_num
  integer :: dim
  logical :: has_2d_
  character(len=20), save :: cls_name = "field_rho"
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = "init_field_rho"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( part_shape )

    case ( p_ps_linear )

      gc_num(:,1) = (/1, 1/)
      gc_num(:,2) = (/0, 1/)

    case ( p_ps_quadratic )
        
    !   call write_err( "Quadratic particle shape not implemented." )
        gc_num(:,1) = (/1, 2/)
        gc_num(:,2) = (/0, 2/)
        
    case default
      call write_err( "Invalid particle shape." )

  end select

  dim = 1
  has_2d_ = .true.
  if ( present(has_2d) ) has_2d_ = has_2d

  call this%field%new( opts, dim, max_mode, gc_num, has_2d=has_2d_ )
  this%smooth_order = smooth_order

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_rho

subroutine smooth_field_rho(this)

  implicit none
  class( field_rho ), intent(inout) :: this

  integer :: mode, k
  character(len=20), save :: cls_name = "field_rho"
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = "smooth_field_rho"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do k = 1, this%smooth_order
    call this%rf_re(0)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.true.])
    do mode = 1, this%max_mode
      call this%rf_re(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.false.])
      call this%rf_im(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.false.])
    enddo
    call this%copy_gc_f1()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine smooth_field_rho

subroutine init_field_jay(this, opts, max_mode, part_shape, smooth_order)

  implicit none

  class( field_jay ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: max_mode, part_shape
  integer, intent(in) :: smooth_order

  integer, dimension(2,2) :: gc_num
  integer :: dim
  character(len=20), save :: cls_name = "field_jay"
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = "init_field_jay"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( part_shape )

    case ( p_ps_linear )

      gc_num(:,1) = (/1, 1/)
      gc_num(:,2) = (/0, 1/)

    case ( p_ps_quadratic )
        
      gc_num(:,1) = (/1, 2/)
      gc_num(:,2) = (/0, 2/)
        
    !   call write_err( "Quadratic particle shape not implemented." )

    case default

      call write_err( "Invalid particle shape." )

  end select

  dim = 3
  call this%field%new( opts, dim, max_mode, gc_num )
  this%smooth_order = smooth_order

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_jay

subroutine smooth_field_jay(this)

  implicit none
  class( field_jay ), intent(inout) :: this

  integer :: mode, k
  character(len=20), save :: cls_name = "field_jay"
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = "smooth_field_jay"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ! TODO: note that the amu is also the instance of field_jay class, but amu corrects the
  ! on-axis value for m=2 mode. Therefore, this routine may not work for amu smooth.
  do k = 1, this%smooth_order
    call this%rf_re(0)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.false., .false., .true.])
    do mode = 1, this%max_mode
      if (mode == 1) then
        call this%rf_re(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.true., .true., .false.])
        call this%rf_im(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.true., .true., .false.])
      else
        call this%rf_re(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.false., .false., .false.])
        call this%rf_im(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.false., .false., .false.])
      endif
    enddo
    call this%copy_gc_f1()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine smooth_field_jay

subroutine init_field_djdxi( this, opts, max_mode, part_shape, smooth_order )

  implicit none

  class( field_djdxi ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: max_mode, part_shape
  integer, intent(in) :: smooth_order

  integer, dimension(2,2) :: gc_num
  integer :: dim
  character(len=20), save :: cls_name = "field_djdxi"
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = "init_field_djdxi"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( part_shape )

  case ( p_ps_linear )

    gc_num(:,1) = (/1, 1/)
    gc_num(:,2) = (/0, 1/)

  case ( p_ps_quadratic )
    
    gc_num(:,1) = (/1, 2/)
    gc_num(:,2) = (/0, 1/)
      
    ! call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 2
  call this%field%new( opts, dim, max_mode, gc_num )
  this%smooth_order = smooth_order

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_djdxi

subroutine smooth_field_djdxi(this)

  implicit none
  class( field_djdxi ), intent(inout) :: this

  integer :: mode, k
  character(len=20), save :: cls_name = "field_djdxi"
  integer, parameter :: cls_level = 3
  character(len=20), save :: sname = "smooth_field_djdxi"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do k = 1, this%smooth_order
    call this%rf_re(0)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.false., .false.])
    do mode = 1, this%max_mode
      if (mode == 1) then
        call this%rf_re(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.true., .true.])
        call this%rf_im(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.true., .true.])
      else
        call this%rf_re(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.false., .false.])
        call this%rf_im(mode)%smooth_f1(stencil=[1.0, 2.0, 1.0], ax_smooth=[.false., .false.])
      endif
    enddo
    call this%copy_gc_f1()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine smooth_field_djdxi

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
  real :: idr, idrh, ir
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

  noff   = this%rf_re(0)%get_noff(1)
  nvp    = num_procs_loc()
  idproc = id_proc_loc()

  do mode = 0, this%max_mode

    uacu_re => acu_re(mode)%get_f1()
    uamu_re => amu_re(mode)%get_f1()
    udcu_re => this%rf_re(mode)%get_f1()

    if ( mode == 0 ) then
      do i = 2, nrp
        ir = idr / real(i+noff-1)
        udcu_re(1,i) = uacu_re(1,i) - idrh * ( uamu_re(1,i+1) - uamu_re(1,i-1) ) - ir * uamu_re(1,i)
        udcu_re(2,i) = uacu_re(2,i) - idrh * ( uamu_re(2,i+1) - uamu_re(2,i-1) ) - ir * uamu_re(2,i)
      enddo

      ! calculate the first cell
      if ( idproc == 0 ) then
        udcu_re(1,1) = 0.0
        udcu_re(2,1) = 0.0
      else
        ir = idr / real(noff)
        udcu_re(1,1) = uacu_re(1,1) - idrh * ( uamu_re(1,2) - uamu_re(1,0) ) - ir * uamu_re(1,1)
        udcu_re(2,1) = uacu_re(2,1) - idrh * ( uamu_re(2,2) - uamu_re(2,0) ) - ir * uamu_re(2,1)
      endif

      ! outer boundary
      if ( idproc == nvp-1 ) then
        ir = idr / real(nrp+noff-1)
        udcu_re(1,nrp) = uacu_re(1,nrp) + idrh * ( 4.0 * uamu_re(1,nrp-1) - &
          uamu_re(1,nrp-2) - 3.0 * uamu_re(1,nrp) ) - ir * uamu_re(1,nrp)
        udcu_re(2,nrp) = uacu_re(2,nrp) + idrh * ( 4.0 * uamu_re(2,nrp-1) - &
          uamu_re(2,nrp-2) - 3.0 * uamu_re(2,nrp) ) - ir * uamu_re(2,nrp)
      endif
      cycle
    endif

    uacu_im => acu_im(mode)%get_f1()
    uamu_im => amu_im(mode)%get_f1()
    udcu_im => this%rf_im(mode)%get_f1()

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      udcu_re(1,i) = uacu_re(1,i) - idrh * ( uamu_re(1,i+1) - uamu_re(1,i-1) ) - ir * uamu_re(1,i) + mode * ir * uamu_im(2,i)
      udcu_re(2,i) = uacu_re(2,i) - idrh * ( uamu_re(2,i+1) - uamu_re(2,i-1) ) - ir * uamu_re(2,i) + mode * ir * uamu_im(3,i)
      udcu_im(1,i) = uacu_im(1,i) - idrh * ( uamu_im(1,i+1) - uamu_im(1,i-1) ) - ir * uamu_im(1,i) - mode * ir * uamu_re(2,i)
      udcu_im(2,i) = uacu_im(2,i) - idrh * ( uamu_im(2,i+1) - uamu_im(2,i-1) ) - ir * uamu_im(2,i) - mode * ir * uamu_re(3,i)
    enddo

    ! calculate the first cell
    if ( idproc == 0 ) then
      if ( mode == 1 ) then
        udcu_re(1,1) = uacu_re(1,1) - 2.0 * idr * uamu_re(1,2) + mode * idr * uamu_im(2,2)
        udcu_re(2,1) = uacu_re(2,1) - 2.0 * idr * uamu_re(2,2) + mode * idr * uamu_im(3,2)
        udcu_im(1,1) = uacu_im(1,1) - 2.0 * idr * uamu_im(1,2) - mode * idr * uamu_re(2,2)
        udcu_im(2,1) = uacu_im(2,1) - 2.0 * idr * uamu_im(2,2) - mode * idr * uamu_re(3,2)
      else
        udcu_re(1,1) = 0.0
        udcu_re(2,1) = 0.0
        udcu_im(1,1) = 0.0
        udcu_im(2,1) = 0.0

        ! since amu(m=2) is multiplied by factor 8 on axis, the derivative on index=2 is
        ! calculated using forward difference
        if ( mode == 2 ) then
          ir = idr
          udcu_re(1,2) = uacu_re(1,2) - idr * ( uamu_re(1,3) - uamu_re(1,2) ) - ir * uamu_re(1,2) + mode * ir * uamu_im(2,2)
          udcu_re(2,2) = uacu_re(2,2) - idr * ( uamu_re(2,3) - uamu_re(2,2) ) - ir * uamu_re(2,2) + mode * ir * uamu_im(3,2)
          udcu_im(1,2) = uacu_im(1,2) - idr * ( uamu_im(1,3) - uamu_im(1,2) ) - ir * uamu_im(1,2) - mode * ir * uamu_re(2,2)
          udcu_im(2,2) = uacu_im(2,2) - idr * ( uamu_im(2,3) - uamu_im(2,2) ) - ir * uamu_im(2,2) - mode * ir * uamu_re(3,2)
        endif
      endif
    else
      ir = idr / real(noff)
      udcu_re(1,1) = uacu_re(1,1) - idrh * ( uamu_re(1,2) - uamu_re(1,0) ) - ir * uamu_re(1,1) + mode * ir * uamu_im(2,1)
      udcu_re(2,1) = uacu_re(2,1) - idrh * ( uamu_re(2,2) - uamu_re(2,0) ) - ir * uamu_re(2,1) + mode * ir * uamu_im(3,1)
      udcu_im(1,1) = uacu_im(1,1) - idrh * ( uamu_im(1,2) - uamu_im(1,0) ) - ir * uamu_im(1,1) - mode * ir * uamu_re(2,1)
      udcu_im(2,1) = uacu_im(2,1) - idrh * ( uamu_im(2,2) - uamu_im(2,0) ) - ir * uamu_im(2,1) - mode * ir * uamu_re(3,1)
    endif

    ! outer boundary
    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-1)
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