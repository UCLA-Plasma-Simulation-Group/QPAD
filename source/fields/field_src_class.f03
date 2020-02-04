module field_src_class

use field_class
use parallel_pipe_class
use grid_class
use ufield_class
use ufield_smooth_class
use param
use sysutil
use mpi

implicit none

private

public :: field_rho, field_jay, field_djdxi

type, extends( field ) :: field_rho

  real :: q_ax ! on-axis charge

  contains

  generic :: new => init_field_rho
  ! procedure :: get_q_ax1, get_q_ax2
  ! procedure :: get_q_ax ! get the on-axis charge according to charge conservation
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

subroutine init_field_rho( this, pp, gp, num_modes, part_shape, &
  smooth_type, smooth_order )

  implicit none

  class( field_rho ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape
  integer, intent(in), optional :: smooth_type, smooth_order

  integer, dimension(2,2) :: gc_num
  integer :: dim
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

    call this%field%new( pp, gp, dim, num_modes, gc_num, &
      smooth_type=smooth_type, smooth_order=smooth_order )

  else

    call this%field%new( pp, gp, dim, num_modes, gc_num )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_rho

! subroutine get_q_ax( this )

!   implicit none
!   class( field_rho ), intent(inout) :: this

!   real, dimension(:,:), pointer :: f1_re => null()
!   integer :: dtype, comm, idproc, ierr, i, nrp, noff
!   real :: r, q_local

!   idproc = this%rf_re(0)%pp%getlidproc()
!   comm   = this%rf_re(0)%pp%getlgrp()
!   dtype  = this%rf_re(0)%pp%getmreal()
!   nrp    = this%rf_re(0)%get_ndp(1)
!   noff   = this%rf_re(0)%get_noff(1)

!   ! if ( idproc == 0 ) then
!   !   f1_re => this%rf_re(0)%get_f1()
!   !   this%q_ax = -1.0 * this%dr**2 * f1_re(1,0)
!   ! endif

!   ! call MPI_BCAST( this%q_ax, 1, dtype, 0, comm, ierr )

!   f1_re => this%rf_re(0)%get_f1()
!   q_local = 0.0
!   do i = 1, nrp
!     r = real( noff + i - 0.5 )
!     q_local = q_local + f1_re(1,i) * r
!   enddo
!   ! q_local = q_local * this%dr**2

!   call MPI_ALLREDUCE( q_local, this%q_ax, 1, dtype, MPI_SUM, comm, ierr )

! end subroutine get_q_ax

! subroutine get_q_ax2( this )

!   implicit none
!   class( field_rho ), intent(inout) :: this
!   ! class( ufield_smooth ), intent(in) :: smooth

!   real, dimension(:,:), pointer :: f1_re => null()
!   integer :: dtype, comm, idproc, nsm, i, j, ierr
!   real, dimension(:), pointer :: scoef
!   real :: dr2

!   idproc = this%rf_re(0)%pp%getlidproc()
!   comm   = this%rf_re(0)%pp%getlgrp()
!   dtype  = this%rf_re(0)%pp%getmreal()
!   nsm    = this%smooth%get_order()
!   dr2    = this%dr**2
!   scoef  => this%smooth%get_scoef()

!   if ( nsm <= 0 ) return

!   if ( idproc == 0 ) then

!     f1_re => this%rf_re(0)%get_f1()

!     ! add the charge smoothed out
!     do j = 1, nsm
!       do i = -nsm, -j
!         this%q_ax = this%q_ax - dr2 * f1_re(1,j) * scoef(i) * real(i+j-0.5)
!       enddo
!     enddo

!     ! deduct the charge smoothed in
!     do j = 1, nsm
!       do i = 1, nsm-j+1
!         this%q_ax = this%q_ax + dr2 * real(j-0.5) * f1_re(1,i) * scoef(1-j-i)
!       enddo
!     enddo

!   endif

!   call MPI_BCAST( this%q_ax, 1, dtype, 0, comm, ierr )

! end subroutine get_q_ax2

subroutine init_field_jay( this, pp, gp, num_modes, part_shape, &
  smooth_type, smooth_order )

  implicit none

  class( field_jay ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape
  integer, intent(in), optional :: smooth_type, smooth_order

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

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 3
  ! call initialization routine of the parent class
  if ( present(smooth_type) .and. present(smooth_order) ) then

    gc_num(1,1) = max( gc_num(1,1), smooth_order )
    gc_num(2,1) = max( gc_num(2,1), smooth_order )

    call this%field%new( pp, gp, dim, num_modes, gc_num, &
      smooth_type=smooth_type, smooth_order=smooth_order )

  else

    call this%field%new( pp, gp, dim, num_modes, gc_num )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_jay

subroutine init_field_djdxi( this, pp, gp, num_modes, part_shape, &
  smooth_type, smooth_order )

  implicit none

  class( field_djdxi ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape
  integer, intent(in), optional :: smooth_type, smooth_order

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

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 2
  ! call initialization routine of the parent class
  if ( present(smooth_type) .and. present(smooth_order) ) then

    gc_num(1,1) = max( gc_num(1,1), smooth_order )
    gc_num(2,1) = max( gc_num(2,1), smooth_order )

    call this%field%new( pp, gp, dim, num_modes, gc_num, &
      smooth_type=smooth_type, smooth_order=smooth_order )

  else

    call this%field%new( pp, gp, dim, num_modes, gc_num )

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

  noff = this%rf_re(0)%get_noff(1)
  nvp = this%rf_re(0)%pp%getlnvp()
  idproc = this%rf_re(0)%pp%getlidproc()

  do mode = 0, this%num_modes

    uacu_re => acu_re(mode)%get_f1()
    uamu_re => amu_re(mode)%get_f1()
    udcu_re => this%rf_re(mode)%get_f1()

    if ( mode == 0 ) then
      do i = 2, nrp
        ir = idr / real(i+noff-1)
        udcu_re(1,i) = uacu_re(1,i) - idrh * ( uamu_re(1,i+1) - uamu_re(1,i-1) ) - ir * uamu_re(1,i)
        udcu_re(2,i) = uacu_re(2,i) - idrh * ( uamu_re(2,i+1) - uamu_re(2,i-1) ) - ir * uamu_re(2,i)
      enddo

      ! on axis
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
    endif

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