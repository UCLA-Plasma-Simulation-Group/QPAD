#ifndef __TEMPLATE__

module part2d_class

use param
use sysutil_module
use parallel_module
use sort_module
use options_class
use field_class
use field_laser_class
use ufield_class
use fdist2d_class
use hdf5io_class
use mpi
use interpolation
use interp_part2d
use debug_tool

implicit none

private

public :: part2d

type part2d

  ! qbm = particle charge/mass ratio
  ! dt = time interval between successive calculations
  ! dr = radial cell size
  real :: qbm, dt, dr

  ! npp = number of particles in current partition
  ! npmax = maximum number of particles in each partition
  integer(kind=LG) :: npmax, npp

  ! dimension of particle coordinates
  integer :: part_dim

  ! array for particle position
  real, dimension(:,:), allocatable :: x
  ! array for particle momenta
  real, dimension(:,:), allocatable :: p
  ! array for time-centered gamma
  real, dimension(:), allocatable :: gamma
  ! array for particle charge
  real, dimension(:), allocatable :: q
  ! array for psi
  real, dimension(:), allocatable :: psi
  ! particle upper boundaries
  real :: edge
  ! particle buffer
  real, dimension(:), allocatable :: pbuf
  ! on-axis deposition correction factor
  real, dimension(1,2) :: deposit_ax_corr
  integer :: part_shape

  ! temporary arrays used for buffer reallocation
  real, private, dimension(:), allocatable :: tmp1
  real, private, dimension(:,:), allocatable :: tmp2

  contains

  procedure :: new                   => init_part2d
  procedure :: renew                 => renew_part2d
  procedure :: del                   => end_part2d
  procedure :: qdeposit              => qdeposit_part2d
  procedure :: deposit_chi           => deposit_chi_part2d
  procedure :: amjdeposit_std        => amjdeposit_std_part2d
  procedure :: amjdeposit_robust     => amjdeposit_robust_part2d
  procedure :: amjdeposit_std_pgc    => amjdeposit_std_pgc_part2d
  procedure :: amjdeposit_robust_pgc => amjdeposit_robust_pgc_part2d
  procedure :: push_u_std            => push_u_std_part2d
  procedure :: push_u_robust         => push_u_robust_part2d
  procedure :: push_u_std_pgc        => push_u_std_pgc_part2d
  procedure :: push_u_robust_pgc     => push_u_robust_pgc_part2d
  procedure :: push_x                => push_x_part2d
  procedure :: update_bound          => update_bound_part2d
  procedure :: pipesend              => pipesend_part2d
  procedure :: piperecv              => piperecv_part2d
  procedure :: wr                    => writehdf5_part2d
  procedure :: realloc               => realloc_part2d
  procedure :: sort                  => sort_part2d
  procedure :: interp_psi            => interp_psi_part2d

end type

save

character(len=20), parameter :: cls_name = "part2d"
integer, parameter :: cls_level = 2

real, dimension(:), allocatable :: recv_buf
integer :: recv_buf_size = 0

integer, dimension(:), allocatable :: sort_idx

integer, parameter :: p_max_subcyc = 1024
real, parameter :: p_buf_incr = 1.5

contains

subroutine init_part2d( this, opts, pf, qbm, dt, s, part_shape, if_empty )

   implicit none

   class(part2d), intent(inout) :: this
   type(options), intent(in) :: opts
   class(fdist2d), intent(inout) :: pf
   real, intent(in) :: qbm, dt, s
   integer, intent(in) :: part_shape
   logical, intent(in), optional :: if_empty
   ! local data
   character(len=18), save :: sname = 'init_part2d'
   integer :: npmax
   logical :: empty = .false.

   call write_dbg( cls_name, sname, cls_level, 'starts' )

   this%qbm = qbm
   this%dt  = dt
   this%part_dim = 2 + p_p_dim + 3

   npmax      = pf%npmax
   this%npmax = npmax
   this%npp   = 0

   this%dr   = opts%get_dr()
   this%edge = opts%get_nd(1) * this%dr
   this%part_shape = part_shape
   this%deposit_ax_corr = get_deposit_ax_corr( pf%ppc(1), part_shape )
   
   allocate( this%x( 2, npmax ) )
   allocate( this%p( p_p_dim, npmax ) )
   allocate( this%gamma( npmax ), this%q( npmax ), this%psi( npmax ) )
   allocate( this%pbuf( this%part_dim * npmax ) )

   recv_buf_size = max( recv_buf_size, npmax )

   if ( present( if_empty ) ) empty = if_empty

   ! initialize particle coordinates according to specified profile
   if ( .not. empty ) call pf%inject( this%x, this%p, this%gamma, this%psi, this%q, this%npp, s )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_part2d

subroutine end_part2d(this)

   implicit none

   class(part2d), intent(inout) :: this
   character(len=18), save :: sname = 'end_part2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   deallocate( this%x, this%p, this%gamma, this%q, this%psi )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_part2d

subroutine realloc_part2d( this, ratio )

  implicit none

  class(part2d), intent(inout) :: this
  real, intent(in) :: ratio

  integer :: i, npmax
  character(len=18), save :: sname = 'realloc_part2d'

  call write_dbg(cls_name, sname, cls_level, 'starts')

  this%npmax = int( this%npmax * ratio )
  npmax = this%npmax

  allocate( this%tmp2( 2, npmax ) )
  this%tmp2 = 0.0
  this%tmp2( 1:2, 1:this%npp ) = this%x( 1:2, 1:this%npp )
  call move_alloc( this%tmp2, this%x )

  allocate( this%tmp2( p_p_dim, npmax ) )
  this%tmp2 = 0.0
  this%tmp2( 1:p_p_dim, 1:this%npp ) = this%p( 1:p_p_dim, 1:this%npp )
  call move_alloc( this%tmp2, this%p )

  allocate( this%tmp1( npmax ) )
  this%tmp1 = 0.0
  this%tmp1( 1:this%npp ) = this%gamma( 1:this%npp )
  call move_alloc( this%tmp1, this%gamma )

  allocate( this%tmp1( npmax ) )
  this%tmp1 = 0.0
  this%tmp1( 1:this%npp ) = this%psi( 1:this%npp )
  call move_alloc( this%tmp1, this%psi )

  allocate( this%tmp1( npmax ) )
  this%tmp1 = 0.0
  this%tmp1( 1:this%npp ) = this%q( 1:this%npp )
  call move_alloc( this%tmp1, this%q )

  deallocate( this%pbuf )
  allocate( this%pbuf( this%part_dim * npmax ) )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine realloc_part2d

subroutine renew_part2d( this, pf, s, if_empty )

   implicit none

   class(part2d), intent(inout) :: this
   class(fdist2d), intent(inout) :: pf
   real, intent(in) :: s
   logical, intent(in), optional :: if_empty

   ! local data
   logical :: empty
   character(len=18), save :: sname = 'renew_part2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%npp = 0
   empty = .false.   
   if ( present( if_empty ) ) empty = if_empty
   if ( .not. empty ) call pf%inject( this%x, this%p, this%gamma, this%psi, this%q, this%npp, s )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine renew_part2d

subroutine qdeposit_part2d( this, q )
    implicit none
    class(part2d), intent(in) :: this
    class(field), intent(in) :: q
    
    select case (this%part_shape)
    case(p_ps_linear)
        call qdeposit_part2d_s1(this, q)
    case(p_ps_quadratic)
        call qdeposit_part2d_s2(this, q)
    case default
        call write_err('Not implemented yet')
    end select
    
end subroutine qdeposit_part2d

subroutine deposit_chi_part2d( this, chi)
    implicit none
    class(part2d), intent(in) :: this
    class(field), intent(inout) :: chi
    
    select case (this%part_shape)
    case(p_ps_linear)
        call deposit_chi_part2d_s1(this, chi)
    case(p_ps_quadratic)
        call deposit_chi_part2d_s2(this, chi)
    case default
        call write_err('Not implemented yet')
    end select
    
end subroutine deposit_chi_part2d

subroutine amjdeposit_std_part2d(this, ef, bf, cu, amu, dcu, dt_)
    ! deposit the current, acceleration and momentum flux
    
    implicit none
    
    class(part2d), intent(inout) :: this
    class(field), intent(in) :: cu, amu, dcu
    class(field), intent(in) :: ef, bf
    real, intent(in), optional :: dt_

    select case(this%part_shape)
    case(p_ps_linear)
      call amjdeposit_std_part2d_s1(this, ef, bf, cu, amu, dcu, dt_)
    case(p_ps_quadratic)
      call amjdeposit_std_part2d_s2(this, ef, bf, cu, amu, dcu, dt_)
    case default
      call write_err('Not implemented yet')
    end select
  
end subroutine amjdeposit_std_part2d

subroutine amjdeposit_robust_part2d(this, ef, bf, cu, amu, dcu, dt_)
    ! deposit the current, acceleration and momentum flux
    
    implicit none
    
    class(part2d), intent(inout) :: this
    class(field), intent(in) :: cu, amu, dcu
    class(field), intent(in) :: ef, bf
    real, intent(in), optional :: dt_

    select case(this%part_shape)
    case(p_ps_linear)
      call amjdeposit_robust_part2d_s1(this, ef, bf, cu, amu, dcu, dt_)
    case(p_ps_quadratic)
      call amjdeposit_robust_part2d_s2(this, ef, bf, cu, amu, dcu, dt_)
    case default
      call write_err('Not implemented yet')
    end select
  
end subroutine amjdeposit_robust_part2d


subroutine amjdeposit_std_pgc_part2d(this, ef, bf, af, cu, amu, dcu, dt_)
    ! deposit the current, acceleration and momentum flux
    
      implicit none
    
      class(part2d), intent(inout) :: this
      class(field), intent(in) :: cu, amu, dcu
      class(field), intent(in) :: ef, bf
      class(field_laser), intent(in) :: af
      real, intent(in), optional :: dt_
      
      select case(this%part_shape)
      case(p_ps_linear)
        call amjdeposit_std_pgc_part2d_s1(this, ef, bf, af, cu, amu, dcu, dt_)
      case(p_ps_quadratic)
        call amjdeposit_std_pgc_part2d_s2(this, ef, bf, af, cu, amu, dcu, dt_)
      case default
        call write_err('Not implemented yet')
      end select
      
    end subroutine amjdeposit_std_pgc_part2d
    
subroutine amjdeposit_robust_pgc_part2d(this, ef, bf, af, cu, amu, dcu, dt_)
! deposit the current, acceleration and momentum flux

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: cu, amu, dcu
  class(field), intent(in) :: ef, bf
  class(field_laser), intent(in) :: af
  real, intent(in), optional :: dt_
  
  select case(this%part_shape)
  case(p_ps_linear)
    call amjdeposit_robust_pgc_part2d_s1(this, ef, bf, af, cu, amu, dcu, dt_)
  case(p_ps_quadratic)
    call amjdeposit_robust_pgc_part2d_s2(this, ef, bf, af, cu, amu, dcu, dt_)
  case default
    call write_err('Not implemented yet')
  end select
  
end subroutine amjdeposit_robust_pgc_part2d

! subroutine interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, x, dr, bp, ep, np, ptrcur, &
!   geom, weight, ix, pcos, psin )

!   implicit none

!   type(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im, bf_re, bf_im
!   integer, intent(in) :: max_mode, np, geom
!   real, intent(in) :: dr
!   real, dimension(:,:), intent(in) :: x
!   real, dimension(:,:), intent(inout) :: bp, ep
!   integer(kind=LG), intent(in) :: ptrcur
!   real, dimension(:,:), intent(out), optional :: weight
!   integer, dimension(:), intent(out), optional :: ix
!   real, dimension(:), intent(out), optional :: pcos, psin

!   real, dimension(:,:), pointer :: e0, b0, er, ei, br, bi
!   integer :: noff, i, j, nn, mode
!   integer(kind=LG) :: pp
!   real :: pos, idr, ph_r, ph_i, cc, ss
!   real, dimension(0:1) :: wt
!   complex(kind=DB) :: phase0, phase

!   idr = 1.0 / dr
!   noff = ef_re(0)%get_noff(1)

!   e0 => ef_re(0)%get_f1()
!   b0 => bf_re(0)%get_f1()

!   ep = 0.0
!   bp = 0.0

!   pp = ptrcur
!   do i = 1, np
!     pos = sqrt( x(1,pp)**2 + x(2,pp)**2 ) * idr
!     ! cosine and sine
!     cc = x(1,pp) / pos * idr
!     ss = x(2,pp) / pos * idr
!     phase0 = cmplx( cc, ss )

!     ! in-cell position
!     nn  = int( pos )
!     pos = pos - real(nn)

!     ! cell index
!     nn = nn - noff + 1

!     ! get interpolation weight factor
!     call spline_linear( pos, wt )

!     if ( present(weight) ) weight(:,i) = wt
!     if ( present(ix) ) ix(i) = nn
!     if ( present(pcos) ) then
!       pcos(i) = cc
!       psin(i) = ss
!     endif

!     ! interpolate m=0 mode
!     do j = 0, 1
!       ep(:,i) = ep(:,i) + e0(:,nn+j) * wt(j)
!       bp(:,i) = bp(:,i) + b0(:,nn+j) * wt(j)
!     enddo

!     ! interpolate m>0 modes
!     phase = cmplx( 1.0, 0.0 )
!     do mode = 1, max_mode
!         phase = phase * phase0
!         ph_r = 2.0 * real(phase)
!         ph_i = 2.0 * aimag(phase)

!         er => ef_re(mode)%get_f1()
!         ei => ef_im(mode)%get_f1()
!         br => bf_re(mode)%get_f1()
!         bi => bf_im(mode)%get_f1()

!         do j = 0, 1
!           ep(:,i) = ep(:,i) + ( er(:,nn+j) * ph_r - ei(:,nn+j) * ph_i ) * wt(j)
!           bp(:,i) = bp(:,i) + ( br(:,nn+j) * ph_r - bi(:,nn+j) * ph_i ) * wt(j)
!         enddo
!     enddo

!     ! transform from cylindrical geometry to Cartesian geometry
!     if ( geom == p_cartesian ) then
!       ! ph_r, ph_i are temporary variables here
!       ph_r = ep(1,i) * cc - ep(2,i) * ss
!       ph_i = ep(1,i) * ss + ep(2,i) * cc
!       ep(1,i) = ph_r
!       ep(2,i) = ph_i

!       ph_r = bp(1,i) * cc - bp(2,i) * ss
!       ph_i = bp(1,i) * ss + bp(2,i) * cc
!       bp(1,i) = ph_r
!       bp(2,i) = ph_i
!     endif

!     pp = pp + 1
!   enddo

! end subroutine interp_emf_part2d

! subroutine interp_laser_part2d( ar_re, ar_im, ai_re, ai_im, ar_grad_re, ar_grad_im, &
!   ai_grad_re, ai_grad_im, max_mode, x, dr, apr, api, apr_grad, api_grad, np, ptrcur, geom )

!    implicit none

!    type(ufield), dimension(:), pointer, intent(in) :: ar_re, ar_im, ai_re, ai_im
!    type(ufield), dimension(:), pointer, intent(in) :: ar_grad_re, ar_grad_im, ai_grad_re, ai_grad_im
!    integer, intent(in) :: max_mode, np, geom
!    real, intent(in) :: dr
!    real, dimension(:,:), intent(in) :: x
!    real, dimension(:), intent(inout) :: apr, api
!    real, dimension(:,:), intent(inout) :: apr_grad, api_grad
!    integer(kind=LG), intent(in) :: ptrcur

!    integer :: noff, i, j, nn, m
!    integer(kind=LG) :: pp
!    real :: pos, idr, ph_r, ph_i, cc, ss
!    real, dimension(0:1) :: wt
!    complex :: phase0, phase

!    idr = 1.0 / dr
!    noff = ar_re(0)%get_noff(1)
!    apr = 0.0
!    api = 0.0
!    apr_grad = 0.0
!    api_grad = 0.0

!    pp = ptrcur
!    do i = 1, np
!       pos = sqrt( x(1,pp)**2 + x(2,pp)**2 ) * idr
!       ! cosine and sine
!       cc = x(1,pp) / pos * idr
!       ss = x(2,pp) / pos * idr
!       phase0 = cmplx( cc, ss )

!       ! in-cell position
!       nn  = int( pos )
!       pos = pos - real(nn)

!       ! cell index
!       nn = nn - noff + 1

!       ! get interpolation weight factor
!       call spline_linear( pos, wt )

!       ! interpolate m = 0 mode
!       do j = 0, 1
!         apr(i) = apr(i) + ar_re(0)%f1(1,nn+j) * wt(j)
!         api(i) = api(i) + ai_re(0)%f1(1,nn+j) * wt(j)
!         apr_grad(:,i) = apr_grad(:,i) + ar_grad_re(0)%f1(:,nn+j) * wt(j)
!         api_grad(:,i) = api_grad(:,i) + ai_grad_re(0)%f1(:,nn+j) * wt(j)
!       enddo

!       ! interpolate m > 0 modes
!       phase = cmplx( 1.0, 0.0 )
!       do m = 1, max_mode
!          phase = phase * phase0
!          ph_r = 2.0 * real(phase)
!          ph_i = 2.0 * aimag(phase)

!          do j = 0, 1
!             apr(i) = apr(i) + ( ar_re(m)%f1(1,nn+j) * ph_r - ar_im(m)%f1(1,nn+j) * ph_i ) * wt(j)
!             api(i) = api(i) + ( ai_re(m)%f1(1,nn+j) * ph_r - ai_im(m)%f1(1,nn+j) * ph_i ) * wt(j)
!             apr_grad(:,i) = apr_grad(:,i) + ( ar_grad_re(m)%f1(:,nn+j) * ph_r - ar_grad_im(m)%f1(:,nn+j) * ph_i ) * wt(j)
!             api_grad(:,i) = api_grad(:,i) + ( ai_grad_re(m)%f1(:,nn+j) * ph_r - ai_grad_im(m)%f1(:,nn+j) * ph_i ) * wt(j)
!          enddo
!       enddo

!       ! transform from cylindrical geometry to Cartesian geometry
!       if ( geom == p_cartesian ) then
!         ! ph_r, ph_i are temporary variables here
!         ph_r = apr_grad(1,i) * cc - apr_grad(2,i) * ss
!         ph_i = apr_grad(1,i) * ss + apr_grad(2,i) * cc
!         apr_grad(1,i) = ph_r
!         apr_grad(2,i) = ph_i

!         ph_r = api_grad(1,i) * cc - api_grad(2,i) * ss
!         ph_i = api_grad(1,i) * ss + api_grad(2,i) * cc
!         api_grad(1,i) = ph_r
!         api_grad(2,i) = ph_i
!       endif

!       pp = pp + 1
!     enddo

! end subroutine interp_laser_part2d

subroutine push_u_std_part2d(this, ef, bf, dt_)
    implicit none
    
    class(part2d), intent(inout) :: this
    class(field), intent(in) :: ef, bf
    real, intent(in), optional :: dt_
  
    select case(this%part_shape)
    case(p_ps_linear)
      call push_u_std_part2d_s1(this, ef, bf, dt_)
    case(p_ps_quadratic)
      call push_u_std_part2d_s2(this, ef, bf, dt_)
    case default
      call write_err('Not implemented yet')
    end select
    
end subroutine push_u_std_part2d

subroutine push_u_robust_part2d(this, ef, bf, dt_)
    implicit none
    
    class(part2d), intent(inout) :: this
    class(field), intent(in) :: ef, bf
    real, intent(in), optional :: dt_
  
    select case(this%part_shape)
    case(p_ps_linear)
      call push_u_robust_part2d_s1(this, ef, bf, dt_)
    case(p_ps_quadratic)
      call push_u_robust_part2d_s2(this, ef, bf, dt_)
    case default
      call write_err('Not implemented yet')
    end select
    
end subroutine push_u_robust_part2d

subroutine push_u_robust_pgc_part2d(this, ef, bf, af, dt_)
  implicit none
  
  class(part2d), intent(inout) :: this
  class(field), intent(in) :: ef, bf
  class(field_laser), intent(in) :: af
  real, intent(in), optional :: dt_

  select case(this%part_shape)
  case(p_ps_linear)
    call push_u_robust_pgc_part2d_s1(this, ef, bf, af, dt_)
  case(p_ps_quadratic)
    call push_u_robust_pgc_part2d_s2(this, ef, bf, af, dt_)
  case default
    call write_err('Not implemented yet')
  end select
  
end subroutine push_u_robust_pgc_part2d

subroutine push_u_std_pgc_part2d(this, ef, bf, af, dt_)

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: ef, bf
  class(field_laser), intent(in) :: af
  real, intent(in), optional :: dt_
  
  select case(this%part_shape)
  case(p_ps_linear)
    call push_u_std_pgc_part2d_s1(this, ef, bf, af, dt_)
  case(p_ps_quadratic)
    call push_u_std_pgc_part2d_s2(this, ef, bf, af, dt_)
  case default
    call write_err('Not implemented yet')
  end select
  
end subroutine push_u_std_pgc_part2d

subroutine push_x_part2d(this, dt_)

  implicit none

  class(part2d), intent(inout) :: this
  real, intent(in), optional :: dt_
  ! local data
  character(len=18), save :: sname = 'push_x_part2d'

  integer :: i, np
  real :: gam, dtc, dt
  integer(kind=LG) :: ptrcur, pp

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'push 2D particles' )

  dt = this%dt
  if (present(dt_)) dt = dt_

  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    pp = ptrcur
    do i = 1, np
      dtc = dt / (this%gamma(pp) - this%p(3, pp))
      this%x(1, pp) = this%x(1, pp) + this%p(1, pp) * dtc
      this%x(2, pp) = this%x(2, pp) + this%p(2, pp) * dtc
      pp = pp + 1
    enddo

  enddo

  call stop_tprof( 'push 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_x_part2d

subroutine interp_psi_part2d(this, psi_re, psi_im, max_mode)

  implicit none
  class(part2d), intent(inout) :: this
  type(ufield), dimension(:), pointer, intent(in) :: psi_re, psi_im
  integer, intent(in) :: max_mode

  ! real, dimension(:, :), pointer :: psi_re_ptr, psi_im_ptr
  integer :: noff, i, np
  integer(kind=LG) :: pp, ptrcur
  ! real :: pos, idr, ph_r, ph_i, cc, ss
  ! real, dimension(0:1) :: wt
  ! complex(kind=DB) :: phase0, phase

  real, dimension(:, :), allocatable :: weight
  real, dimension(p_cache_size) :: pcos, psin, psip
  integer, dimension(p_cache_size) :: idx

  ! idr = 1.0 / this%dr
  noff = psi_re(0)%get_noff(1)
  ! psi_re_ptr => psi_re(0)%get_f1()
  select case(this%part_shape)
  case(p_ps_linear)
    allocate(weight(0:1, p_cache_size))
  case(p_ps_quadratic)
    allocate(weight(-1:1, p_cache_size))
case default
    call write_err('Not implemented yet')
end select

  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    call gen_interp_info(this%x, this%dr, noff, np, ptrcur, weight, idx, pcos, psin)
    call interp_field(psi_re, psi_im, max_mode, weight, idx, pcos, psin, np, psip)

    pp = ptrcur
    do i = 1, np
      this%psi(pp) = psip(i)
    enddo

  enddo
  
end subroutine interp_psi_part2d

subroutine update_bound_part2d( this )

   implicit none

   class(part2d), intent(inout) :: this
   ! local
   integer(kind=LG) :: i
   real :: pos
   character(len=32), save :: sname = "update_bound_part2d"

   call write_dbg(cls_name, sname, cls_level, 'starts')

   if ( this%npp == 0 ) return

   call start_tprof( 'push 2D particles' )

   i = 1

   do while ( i < this%npp )

      pos = sqrt( this%x(1,i)**2 + this%x(2,i)**2 )

      ! check if particle goes out of the physical edge
      if ( pos >= this%edge ) then
         this%x(:,i)   = this%x(:, this%npp)
         this%p(:,i)   = this%p(:, this%npp)
         this%gamma(i) = this%gamma(this%npp)
         this%psi(i)   = this%psi(this%npp)
         this%q(i)     = this%q(this%npp)
         this%npp = this%npp - 1
         cycle
      endif

      i = i + 1
   enddo

   ! deal with the last particle
   pos = sqrt( this%x(1,this%npp)**2 + this%x(2,this%npp)**2 )

   if ( pos >= this%edge ) then
      this%npp = this%npp - 1
   endif

   call stop_tprof( 'push 2D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine update_bound_part2d

subroutine pipesend_part2d(this, tag, id)

  implicit none

  class(part2d), intent(inout) :: this
  integer, intent(in) :: tag
  integer, intent(inout) :: id

  ! local data
  character(len=18), save :: sname = 'pipesend_part2d'
  integer :: des, ierr, i, stride

  call write_dbg(cls_name, sname, cls_level, 'starts')

  des = id_proc() + num_procs_loc()

  if ( des >= num_procs() ) then
    id = MPI_REQUEST_NULL
    call write_dbg(cls_name, sname, cls_level, 'ends')
    return
  endif

  do i = 1, this%npp

    stride = (i-1) * this%part_dim

    this%pbuf( 1 + stride ) = this%x(1,i)
    this%pbuf( 2 + stride ) = this%x(2,i)
    this%pbuf( 3 + stride ) = this%p(1,i)
    this%pbuf( 4 + stride ) = this%p(2,i)
    this%pbuf( 5 + stride ) = this%p(3,i)
    this%pbuf( 6 + stride ) = this%gamma(i)
    this%pbuf( 7 + stride ) = this%psi(i)
    this%pbuf( 8 + stride ) = this%q(i)

  enddo

  ! NOTE: npp*xdim might be larger than MAX_INT32
  call mpi_isend( this%pbuf, int(this%npp * this%part_dim), p_dtype_real, des, tag, &
    comm_world(), id, ierr )

  ! check for errors
  if (ierr /= 0) then
    call write_err('MPI_ISEND failed')
  endif

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine pipesend_part2d

subroutine piperecv_part2d(this, tag)

  implicit none

  class(part2d), intent(inout) :: this
  integer, intent(in) :: tag
  ! local data
  character(len=18), save :: sname = 'piperecv_part2d'
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer :: recv_cnt, src, ierr, i, stride
  real :: ratio

  call write_dbg(cls_name, sname, cls_level, 'starts')

  if ( .not. allocated(recv_buf) ) then
    allocate( recv_buf( this%part_dim * recv_buf_size ) )
  endif

  src = id_proc() - num_procs_loc()

  if (src < 0) then
    call write_dbg(cls_name, sname, cls_level, 'ends')
    return
  endif

  ! probe for incoming message
  call mpi_probe( src, tag, comm_world(), stat, ierr )
  call mpi_get_count( stat, p_dtype_real, recv_cnt, ierr )
  recv_cnt = recv_cnt / this%part_dim

  ! check if the receiving buffer needs to be reallocated
  if ( recv_cnt > recv_buf_size ) then
    call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
          &particle pipeline receiving buffer!', only_root = .false. )
    deallocate( recv_buf )
    recv_buf_size = int( recv_cnt * p_buf_incr )
    allocate( recv_buf( recv_buf_size * this%part_dim ) )
  endif

  ! NOTE: npp*xdim might be larger than MAX_INT32
  call mpi_recv( recv_buf, this%part_dim * recv_buf_size, p_dtype_real, &
    src, tag, comm_world(), MPI_STATUS_IGNORE, ierr )

  ! check if the particle buffer needs to be reallocated
  ratio = real( recv_cnt ) / this%npmax
  if ( ratio > 1.0 ) then
    call write_stdout( '[process ' // num2str(id_proc()) // ']: Resizing 2D &
          &particle buffer!', only_root = .false. )
    call this%realloc( ratio = ratio * p_buf_incr )
  endif

  this%npp = recv_cnt

  do i = 1, this%npp

    stride = (i-1) * this%part_dim

    this%x(1,i)   = recv_buf( 1 + stride )
    this%x(2,i)   = recv_buf( 2 + stride )
    this%p(1,i)   = recv_buf( 3 + stride )
    this%p(2,i)   = recv_buf( 4 + stride )
    this%p(3,i)   = recv_buf( 5 + stride )
    this%gamma(i) = recv_buf( 6 + stride )
    this%psi(i)   = recv_buf( 7 + stride )
    this%q(i)     = recv_buf( 8 + stride )

  enddo

  ! check for errors
  if (ierr /= 0) then
    call write_err('MPI failed')
  endif

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine piperecv_part2d

subroutine writehdf5_part2d(this,file)

  implicit none

  class(part2d), intent(inout) :: this
  class(hdf5file), intent(in) :: file
  ! local data
  character(len=18), save :: sname = 'writehdf5_part2d'
  integer :: ierr

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call pwpart( file, this%x, this%p, this%q, this%npp, 1, ierr )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_part2d

subroutine sort_part2d( this, nrp, noff )

  implicit none
  class(part2d), intent(inout) :: this
  integer, intent(in) :: nrp, noff

  ! local data
  integer :: i
  integer, dimension(:), allocatable :: ix
  real :: idr, pos
  character(len=18), save :: sname = 'sort_part2d'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'sort 2D particles' )

  if ( .not. allocated(sort_idx) ) allocate( sort_idx(this%npmax) )
  if ( this%npmax > size(sort_idx) ) then
    deallocate(sort_idx)
    allocate( sort_idx(this%npmax) )
  endif

  allocate( ix(this%npp) )
  ix = 0

  idr = 1.0 / this%dr
  ! calculate the grid index of particles
  do i = 1, this%npp
    pos = sqrt( this%x(1,i)**2 + this%x(2,i)**2 ) * idr
    ix(i) = floor(pos) - noff + 1
  enddo

  ! generate the sorted indices
  call generate_sort_idx_1d( ix, sort_idx, int(this%npp, kind=4), nrp )

  ! rearrange the particle position
  this%pbuf(1:this%npp) = this%x( 1, 1:this%npp )
  do i = 1, this%npp
    this%x( 1, sort_idx(i) ) = this%pbuf(i)
  enddo
  this%pbuf(1:this%npp) = this%x( 2, 1:this%npp )
  do i = 1, this%npp
    this%x( 2, sort_idx(i) ) = this%pbuf(i)
  enddo

  ! rearrange the particle momentum
  this%pbuf(1:this%npp) = this%p( 1, 1:this%npp )
  do i = 1, this%npp
    this%p( 1, sort_idx(i) ) = this%pbuf(i)
  enddo
  this%pbuf(1:this%npp) = this%p( 2, 1:this%npp )
  do i = 1, this%npp
    this%p( 2, sort_idx(i) ) = this%pbuf(i)
  enddo
  this%pbuf(1:this%npp) = this%p( 3, 1:this%npp )
  do i = 1, this%npp
    this%p( 3, sort_idx(i) ) = this%pbuf(i)
  enddo

  ! rearrange the particle gamma
  this%pbuf(1:this%npp) = this%gamma( 1:this%npp )
  do i = 1, this%npp
    this%gamma( sort_idx(i) ) = this%pbuf(i)
  enddo

  ! rearrange the particle charge
  this%pbuf(1:this%npp) = this%q( 1:this%npp )
  do i = 1, this%npp
    this%q( sort_idx(i) ) = this%pbuf(i)
  enddo

  ! rearrange the particle psi
  this%pbuf(1:this%npp) = this%psi( 1:this%npp )
  do i = 1, this%npp
    this%psi( sort_idx(i) ) = this%pbuf(i)
  enddo

  deallocate( ix )

  call stop_tprof( 'sort 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')
    
end subroutine sort_part2d

function get_deposit_ax_corr( ppc_r, part_shape ) result(res)
    implicit none
    integer, intent(in) :: ppc_r
    real, dimension(1,2) :: res
    integer, intent(in) :: part_shape
    
    res = 0.0
    select case(part_shape)
    case(p_ps_linear)
        res(1,1) = (12.0 * ppc_r**2) / ( 1.0 + 2.0 * ppc_r**2 ) 
    case(p_ps_quadratic)
        if (mod(ppc_r,2)==0) then
            res(1,1) = (384.0 * ppc_r**2) / (14.0 + 79.0 * ppc_r**2)
            res(1,2) = (384.0 * ppc_r**2) / (2.0 + 385.0 * ppc_r**2)
        else
            res(1,1) = (384.0 * ppc_r**4) / (3.0 + 14.0 * ppc_r**2 + 79.0 * ppc_r**4)
            res(1,2) = (384.0 * ppc_r**4) / (-3.0 + 2.0 * ppc_r**2 + 385.0 * ppc_r**4)
        endif
  case default
      call write_err('Not implemented yet')
  end select
  end function get_deposit_ax_corr
  
#define __TEMPLATE__

#define SPLINE spline_linear
#define QDEPOSIT qdeposit_part2d_s1
#define CHIDEPOSIT deposit_chi_part2d_s1
! #define AX_COOR get_deposit_ax_corr_s1
#define LINEAR
! #define CORRNUM 1
#define AMJ_STD amjdeposit_std_part2d_s1
#define AMJ_ROB amjdeposit_robust_part2d_s1
#define AMJ_STD_PGC amjdeposit_std_pgc_part2d_s1
#define AMJ_ROB_PGC amjdeposit_robust_pgc_part2d_s1
#define PUSH_STD push_u_std_part2d_s1
#define PUSH_ROB push_u_robust_part2d_s1
#define PUSH_ROB_PGC push_u_robust_pgc_part2d_s1
#define PUSH_STD_PGC push_u_std_pgc_part2d_s1
#define LP 0
#define UP 1

#include __FILE__

#define SPLINE spline_quadratic
#define QDEPOSIT qdeposit_part2d_s2
#define CHIDEPOSIT deposit_chi_part2d_s2
! #define AX_COOR get_deposit_ax_corr_s2
#define QUADRATIC
! #define CORRNUM 2
#define AMJ_STD amjdeposit_std_part2d_s2
#define AMJ_ROB amjdeposit_robust_part2d_s2
#define AMJ_STD_PGC amjdeposit_std_pgc_part2d_s2
#define AMJ_ROB_PGC amjdeposit_robust_pgc_part2d_s2
#define PUSH_STD push_u_std_part2d_s2
#define PUSH_ROB push_u_robust_part2d_s2
#define PUSH_ROB_PGC push_u_robust_pgc_part2d_s2
#define PUSH_STD_PGC push_u_std_pgc_part2d_s2
#define LP -1
#define UP 1

#include __FILE__

end module part2d_class

#else
    
subroutine QDEPOSIT( this, q )
! deposit the charge density (rho - Jz)

  implicit none

  class(part2d), intent(in) :: this
  class(field), intent(in) :: q
  ! local data
  type(ufield), dimension(:), pointer :: q_re => null(), q_im => null()
  complex(kind=DB) :: phase, phase0
  real, dimension(:,:), pointer :: q0 => null(), qr => null(), qi => null()
  real, dimension(LP:UP) :: wt ! interpolation weight
  real :: pos_norm ! normalized position
  real :: idr, ir
  integer(kind=LG) :: ptrcur, pp
  integer :: i, j, nn, noff, nrp, np, mode, max_mode

  character(len=18), save :: sname = 'qdeposit_part2d'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'deposit 2D particles' )

  idr = 1.0 / this%dr

  q_re => q%get_rf_re()
  q_im => q%get_rf_im()

  max_mode = q%get_max_mode()

  noff = q_re(0)%get_noff(1)
  nrp  = q_re(0)%get_ndp(1)
  q0   => q_re(0)%get_f1()

  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    pp = ptrcur
    do i = 1, np
       ! get normalized position
       pos_norm= sqrt( this%x(1, pp)**2 + this%x(2, pp)**2 ) * idr
       phase0 = cmplx( this%x(1, pp), -this%x(2, pp) ) / pos_norm * idr

        nn = floor( pos_norm )
        ! in-cell position
        pos_norm = pos_norm - real(nn)
        nn = nn - noff + 1

#ifdef QUADRATIC
        nn = nn + floor(pos_norm-0.5) + 1
        pos_norm = pos_norm - floor(pos_norm-0.5) - 1
#endif

        call SPLINE( pos_norm, wt )

        phase = cmplx( 1.0, 0.0 ) * this%q(pp)
        ! deposit m=0 mode
        do j = LP, UP
            q0( 1, nn+j ) = q0( 1, nn+j ) + wt(j) * real(phase)
        enddo
        
        ! deposit m>0 mode
        do mode = 1, max_mode
            qr => q_re(mode)%get_f1()
            qi => q_im(mode)%get_f1()
            phase = phase * phase0
    
            do j = LP, UP
            qr( 1, nn+j ) = qr( 1, nn+j ) + wt(j) * real(phase)
            qi( 1, nn+j ) = qi( 1, nn+j ) + wt(j) * aimag(phase)
            enddo

        enddo
        pp = pp + 1
    enddo

  enddo

  if ( id_proc_loc() == 0 ) then
    
! Correct the charge on axis by adding the value on the guard cell 0
! TODO: the data exchanged between partitions should be 2!!
#ifdef QUADRATIC
    q0(1,1) = q0(1,1) + q0(1,0)
#endif

    q0(1,0) = 0.0 ! guard cell is useless on axis
    q0(1,1) = 8.0 * q0(1,1)
    ! q0(1,1) = q0(1,1) * this%deposit_ax_corr
    do j = 2, nrp + 1
      ir = 1.0 / ( j + noff - 1 )
      q0(1,j) = q0(1,j) * ir
    enddo

    do i = 1, max_mode
      qr => q_re(i)%get_f1()
      qi => q_im(i)%get_f1()
      qr(1,0) = 0.0
      qi(1,0) = 0.0
      qr(1,1) = 0.0
      qi(1,1) = 0.0

      do j = 2, nrp + 1
        ir = 1.0 / ( j + noff - 1 )
        qr(1,j) = qr(1,j) * ir
        qi(1,j) = qi(1,j) * ir
      enddo
    enddo

  else
! #ifdef QUADRATIC
! ! Correct the charge on the boundary by adding the value on the guard cell nrp+2
!     if (id_proc_loc() == num_procs_loc()-1) then
!         q0(1, nrp) = q0(1, nrp) + q0(1, nrp+2)
!         do i = 1, max_mode
!             qr => q_re(i)%get_f1()
!             qi => q_im(i)%get_f1()
!             qr(1, nrp) = qr(1, nrp) + qr(1, nrp+2)
!             qi(1, nrp) = qi(1, nrp) + qi(1, nrp+2)
!         enddo
!     endif
! #endif
    do j = 0, nrp + 1
      ir = 1.0 / ( j + noff - 1 )
      q0(1,j) = q0(1,j) * ir
    enddo

    do i = 1, max_mode
      qr => q_re(i)%get_f1()
      qi => q_im(i)%get_f1()
      do j = 0, nrp + 1
        ir = 1.0 / ( j + noff - 1 )
        qr(1,j) = qr(1,j) * ir
        qi(1,j) = qi(1,j) * ir
      enddo
    enddo

  endif

  call stop_tprof( 'deposit 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine QDEPOSIT


subroutine CHIDEPOSIT( this, chi )
! deposit the plasma susceptibility

  implicit none
  class(part2d), intent(in) :: this
  class(field), intent(inout) :: chi
  ! local data
  type(ufield), dimension(:), pointer :: chi_re => null(), chi_im => null()
  complex(kind=DB) :: phase, phase0
  real, dimension(LP:UP) :: wt ! interpolation weight
  real :: pos_norm ! normalized position
  real :: idr, ir
  integer(kind=LG) :: ptrcur, pp
  integer :: i, j, nn, noff, nrp, np, mode, max_mode
  
  character(len=18), save :: sname = 'deposit_chi_part2d'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'deposit 2D particles' )

  chi_re => chi%get_rf_re()
  chi_im => chi%get_rf_im()
  idr      = 1.0 / this%dr
  max_mode = chi%get_max_mode()
  noff     = chi_re(0)%get_noff(1)
  nrp      = chi_re(0)%get_ndp(1)

  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    pp = ptrcur
    do i = 1, np

      ! get normalized position
      pos_norm = sqrt( this%x(1, pp)**2 + this%x(2, pp)**2 ) * idr
      phase0 = cmplx( this%x(1, pp), -this%x(2, pp) ) / pos_norm * idr

      ! calculate in-cell position
      nn = floor( pos_norm )
      pos_norm = pos_norm - real(nn)
      nn = nn - noff + 1

#ifdef QUADRATIC
      nn = nn + floor(pos_norm-0.5) + 1
      pos_norm = pos_norm - floor(pos_norm-0.5) - 1
#endif

      call SPLINE( pos_norm, wt )

      ! phase = -1.0 * cmplx( 1.0, 0.0 ) * this%qbm * this%q(pp) / this%psi(pp)
      phase = -1.0 * cmplx( 1.0, 0.0 ) * this%qbm * this%q(pp) / (1.0 - this%qbm * this%psi(pp))
      ! deposit m = 0 mode
      do j = LP, UP
        chi_re(0)%f1(1,nn+j) = chi_re(0)%f1(1,nn+j) + wt(j) * real(phase)
      enddo

      ! deposit m > 0 modes
      do mode = 1, max_mode
        phase = phase * phase0
        do j = LP, UP
          chi_re(mode)%f1(1,nn+j) = chi_re(mode)%f1(1,nn+j) + wt(j) * real(phase)
          chi_im(mode)%f1(1,nn+j) = chi_im(mode)%f1(1,nn+j) + wt(j) * aimag(phase)
        enddo
      enddo

      pp = pp + 1
    enddo

  enddo

  if ( id_proc_loc() == 0 ) then
    
! Correct the charge on axis by adding the value on the guard cell 0
#ifdef QUADRATIC
    chi_re(0)%f1(1,1) = chi_re(0)%f1(1,1) + chi_re(0)%f1(1,0)
#endif

    chi_re(0)%f1(1,0) = 0.0 ! guard cell is useless on axis
    ! chi_re(0)%f1(1,1) = 8.0 * chi_re(0)%f1(1,1) ! is this correct???
#ifdef LINEAR
    chi_re(0)%f1(1,1) = chi_re(0)%f1(1,1) * this%deposit_ax_corr(1,1)
#endif

#ifdef QUADRATIC
    chi_re(0)%f1(1,1) = chi_re(0)%f1(1,1) * this%deposit_ax_corr(1,1)
    chi_re(0)%f1(1,2) = chi_re(0)%f1(1,2) * this%deposit_ax_corr(1,2)
#endif

    do j = 2, nrp + 1
      ir = 1.0 / ( j + noff - 1 )
      chi_re(0)%f1(1,j) = chi_re(0)%f1(1,j) * ir
    enddo

    do mode = 1, max_mode
      ! guard cell is useless on axis
      chi_re(mode)%f1(1,0) = 0.0
      chi_im(mode)%f1(1,0) = 0.0
      ! on-axis values are all zeros
      chi_re(mode)%f1(1,1) = 0.0
      chi_im(mode)%f1(1,1) = 0.0
#ifdef QUADRATIC
      chi_re(mode)%f1(1,2) = chi_re(mode)%f1(1,2) * this%deposit_ax_corr(1,2)
      chi_im(mode)%f1(1,2) = chi_im(mode)%f1(1,2) * this%deposit_ax_corr(1,2)
#endif

      do j = 2, nrp + 1
        ir = 1.0 / ( j + noff - 1 )
        chi_re(mode)%f1(1,j) = chi_re(mode)%f1(1,j) * ir
        chi_im(mode)%f1(1,j) = chi_im(mode)%f1(1,j) * ir
      enddo
    enddo

  else
    ! #ifdef QUADRATIC
    ! ! Correct the chi on the boundary by adding the value on the guard cell nrp+2
    !     if (id_proc_loc() == num_procs_loc() - 1) then
    !         chi_re(0)%f1(1,nrp) = chi_re(0)%f1(1,nrp) + chi_re(0)%f1(1,nrp+2)            
    !         do mode = 1, max_mode
    !             chi_re(mode)%f1(1,nrp) = chi_re(mode)%f1(1,nrp) + chi_re(mode)%f1(1,nrp+2)
    !             chi_im(mode)%f1(1,nrp) = chi_im(mode)%f1(1,nrp) + chi_im(mode)%f1(1,nrp+2)
    !         enddo
    !     endif
    ! #endif
    do j = 0, nrp + 1
      ir = 1.0 / ( j + noff - 1 )
      chi_re(0)%f1(1,j) = chi_re(0)%f1(1,j) * ir
    enddo

    do mode = 1, max_mode
      do j = 0, nrp + 1
        ir = 1.0 / ( j + noff - 1 )
        chi_re(mode)%f1(1,j) = chi_re(mode)%f1(1,j) * ir
        chi_im(mode)%f1(1,j) = chi_im(mode)%f1(1,j) * ir
      enddo
    enddo

  endif

  call stop_tprof( 'deposit 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine CHIDEPOSIT


subroutine AMJ_STD(this, ef, bf, cu, amu, dcu, dt_)
    ! deposit the current, acceleration and momentum flux
    
      implicit none
    
      class(part2d), intent(inout) :: this
      class(field), intent(in) :: cu, amu, dcu
      class(field), intent(in) :: ef, bf
      real, intent(in), optional :: dt_
      ! local data
      character(len=32), save :: sname = 'amjdeposit_std_part2d'
      type(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
      type(ufield), dimension(:), pointer :: bf_re => null(), bf_im => null()
      type(ufield), dimension(:), pointer :: cu_re => null(), cu_im => null()
      type(ufield), dimension(:), pointer :: dcu_re => null(), dcu_im => null()
      type(ufield), dimension(:), pointer :: amu_re => null(), amu_im => null()
    
      real, dimension(:,:), pointer :: cu0 => null(), dcu0 => null(), amu0 => null()
      real, dimension(:,:), pointer :: cur => null(), dcur => null(), amur => null()
      real, dimension(:,:), pointer :: cui => null(), dcui => null(), amui => null()
    
      integer(kind=LG) :: ptrcur, pp
      integer :: i, j, noff, nrp, np, mode, max_mode
      integer, dimension(p_cache_size) :: ix
      real, dimension(p_p_dim, p_cache_size) :: bp, ep, wp, u0, u
      real, dimension(:, :), allocatable :: weight
      real, dimension(p_cache_size) :: pcos, psin
      real, dimension(p_p_dim) :: du, u2, utmp
      real :: qtmh, qtmh1, qtmh2, idt, gam, ostq, ipsi, dpsi, w, ir, dt
      complex(kind=DB) :: phase, phase0
    
      call write_dbg(cls_name, sname, cls_level, 'starts')
      call start_tprof( 'deposit 2D particles' )
      allocate(weight(LP:UP, p_cache_size))
    
      ef_re  => ef%get_rf_re();  ef_im  => ef%get_rf_im()
      bf_re  => bf%get_rf_re();  bf_im  => bf%get_rf_im()
      cu_re  => cu%get_rf_re();  cu_im  => cu%get_rf_im()
      dcu_re => dcu%get_rf_re(); dcu_im => dcu%get_rf_im()
      amu_re => amu%get_rf_re(); amu_im => amu%get_rf_im()
    
      dt = this%dt
      if (present(dt_)) dt = dt_
    
      idt = 1.0 / dt
      qtmh = 0.5 * this%qbm * dt
      max_mode = ef%get_max_mode()
    
      noff = cu_re(0)%get_noff(1)
      nrp  = cu_re(0)%get_ndp(1)
    
      cu0  => cu_re(0)%get_f1()
      dcu0 => dcu_re(0)%get_f1()
      amu0 => amu_re(0)%get_f1()
    
      do ptrcur = 1, this%npp, p_cache_size
    
        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > this%npp ) then
          np = this%npp - ptrcur + 1
        else
          np = p_cache_size
        endif
    
        ! interpolate fields to particles
        call gen_interp_info(this%x, this%dr, noff, np, ptrcur, weight, ix, pcos, psin)
        call interp_field(ef_re, ef_im, max_mode, weight, ix, pcos, psin, np, ep)
        call interp_field(bf_re, bf_im, max_mode, weight, ix, pcos, psin, np, bp)
    
        pp = ptrcur
        do i = 1, np
    
          ! calculate wake field
          wp(1, i) = ep(1, i) - bp(2, i)
          wp(2, i) = ep(2, i) + bp(1, i)
          wp(3, i) = ep(3, i)
    
          ! transform momentum from Cartesian to cylindrical coordinates
          u0(1, i) = this%p(1, pp) * pcos(i) + this%p(2, pp) * psin(i)
          u0(2, i) = this%p(2, pp) * pcos(i) - this%p(1, pp) * psin(i)
          u0(3, i) = this%p(3, pp)
    
          ! TODO: There is a strict time-centered scheme and it is left for future work.
    
          ! normalize the E and B fields
          gam = sqrt(1.0 + u0(1, i)**2 + u0(2, i)**2 + u0(3, i)**2)
          qtmh2 = qtmh / (1.0 - this%qbm * this%psi(pp))
          qtmh1 = qtmh2 * gam
          bp(:, i) = bp(:, i) * qtmh2
    
          ! half electric acceleration
          utmp(:) = u0(:, i) + ep(:, i) * qtmh1
    
          ! magnetic rotation
          u(1, i) = utmp(1) + utmp(2) * bp(3, i) - utmp(3) * bp(2, i)
          u(2, i) = utmp(2) + utmp(3) * bp(1, i) - utmp(1) * bp(3, i)
          u(3, i) = utmp(3) + utmp(1) * bp(2, i) - utmp(2) * bp(1, i)
    
          ostq = 2.0 / (1.0 + bp(1, i)**2 + bp(2, i)**2 + bp(3, i)**2)
          bp(:,i) = bp(:,i) * ostq
    
          utmp(1) = utmp(1) + u(2, i) * bp(3, i) - u(3, i) * bp(2, i)
          utmp(2) = utmp(2) + u(3, i) * bp(1, i) - u(1, i) * bp(3, i)
          utmp(3) = utmp(3) + u(1, i) * bp(2, i) - u(2, i) * bp(1, i)
    
          ! half electric acceleration
          gam = sqrt(1.0 + utmp(1)**2 + utmp(2)**2 + utmp(3)**2)
          qtmh1 = qtmh2 * gam
          u(:, i) = utmp(:) + ep(:, i) * qtmh1
    
          pp = pp + 1
        enddo
    
        ! calculate and store time-centered values
        ! deposit momentum flux, acceleration density, and current density
        pp = ptrcur
        do i = 1, np
    
          du(1) = idt * ( u(1,i) - u0(1,i) )
          du(2) = idt * ( u(2,i) - u0(2,i) )
          
          u(:, i)  = 0.5 * ( u(:, i) + u0(:, i) )
          this%gamma(pp) = sqrt(1.0 + u(1, i)**2 + u(2, i)**2 + u(3, i)**2)
    
          ipsi = 1.0 / (1.0 - this%qbm * this%psi(pp))
          dpsi = this%qbm * (wp(3, i) - (wp(1, i) * u(1, i) + wp(2, i) * u(2, i)) * ipsi)
    
          du(1) = du(1) + u(1, i) * dpsi * ipsi
          du(2) = du(2) + u(2, i) * dpsi * ipsi
    
          u2(1) = u(1, i) * u(1, i) * ipsi
          u2(2) = u(1, i) * u(2, i) * ipsi
          u2(3) = u(2, i) * u(2, i) * ipsi
    
          phase0 = cmplx(pcos(i), -psin(i))
          phase  = cmplx(1.0, 0.0) * this%q(pp) * ipsi
    
          ! deposit m = 0 mode
          do j = LP, UP
            w = weight(j, i) * real(phase)
            cu0(1:3, ix(i) + j)  = cu0(1:3, ix(i) + j)  + w * u(1:3, i)
            dcu0(1:2, ix(i) + j) = dcu0(1:2, ix(i) + j) + w * du(1:2)
            amu0(1:3, ix(i) + j) = amu0(1:3, ix(i) + j) + w * u2(1:3)
          enddo
    
          ! deposit m > 0 mode
          do mode = 1, max_mode
    
            cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
            dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
            amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
            phase = phase * phase0
    
            do j = LP, UP
              w = weight(j,i) * real(phase)
              cur( 1:3, ix(i)+j )  = cur( 1:3, ix(i)+j )  + w * u(1:3,i)
              dcur( 1:2, ix(i)+j ) = dcur( 1:2, ix(i)+j ) + w * du(1:2)
              amur( 1:3, ix(i)+j ) = amur( 1:3, ix(i)+j ) + w * u2(1:3)
    
              w = weight(j,i) * aimag(phase)
              cui( 1:3, ix(i)+j )  = cui( 1:3, ix(i)+j )  + w * u(1:3,i)
              dcui( 1:2, ix(i)+j ) = dcui( 1:2, ix(i)+j ) + w * du(1:2)
              amui( 1:3, ix(i)+j ) = amui( 1:3, ix(i)+j ) + w * u2(1:3)
            enddo
    
          enddo
    
          pp = pp + 1
        enddo
    
      enddo
    
      if ( id_proc_loc() ==  0 ) then
! Correct the used values on axis by adding the value on the guard cell 0
#ifdef QUADRATIC
    cu0(3,1) = cu0(3,1) + cu0(3,0)
#endif
        ! guard cells on the axis are useless
        cu0(1:3,0)  = 0.0
        dcu0(1:2,0) = 0.0
        amu0(1:3,0) = 0.0
    
        cu0(1:2,1)  = 0.0; cu0(3,1) = 8.0 * cu0(3,1)
        dcu0(1:2,1) = 0.0
        amu0(1:3,1) = 0.0
    
        do mode = 1, max_mode
    
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          ! guard cells on the axis are useless
          cur(1:3,0)  = 0.0; cui(1:3,0)  = 0.0
          dcur(1:2,0) = 0.0; dcui(1:2,0) = 0.0
          amur(1:3,0) = 0.0; amui(1:3,0) = 0.0
           
          if ( mode == 1 ) then
#ifdef QUADRATIC
          cur(1:2,1) = cur(1:2,1) + cur(1:2,0)
          dcur(1:2,1) = dcur(1:2,1) + dcur(1:2,0)
          cui(1:2,1) = cui(1:2,1) + cui(1:2,0)
          dcui(1:2,1) = dcui(1:2,1) + dcui(1:2,0)
#endif
            cur(1:2,1)  = 8.0 * cur(1:2,1); cur(3,1) = 0.0
            dcur(1:2,1) = 8.0 * dcur(1:2,1)
            amur(1:3,1) = 0.0
            cui(1:2,1)  = 8.0 * cui(1:2,1); cui(3,1) = 0.0
            dcui(1:2,1) = 8.0 * dcui(1:2,1)
            amui(1:3,1) = 0.0
          elseif ( mode == 2 ) then
#ifdef QUADRATIC
            amur(1:3,1) = amur(1:3,1) + amur(1:3,0)
            amui(1:3,1) = amui(1:3,1) + amui(1:3,0)
#endif
            cur(1:3,1)  = 0.0
            dcur(1:2,1) = 0.0
            amur(1:3,1) = 8.0 * amur(1:3,1)
            cui(1:3,1)  = 0.0
            dcui(1:2,1) = 0.0
            amui(1:3,1) = 8.0 * amui(1:3,1)
          else
            cur(1:3,1)  = 0.0
            dcur(1:2,1) = 0.0
            amur(1:3,1) = 0.0
            cui(1:3,1)  = 0.0
            dcui(1:2,1) = 0.0
            amui(1:3,1) = 0.0
           endif
        enddo
    
        do j = 2, nrp + 1
          ir = 1.0 / ( j + noff - 1 )
          cu0(1:3,j)  = cu0(1:3,j)  * ir
          dcu0(1:2,j) = dcu0(1:2,j) * ir
          amu0(1:3,j) = amu0(1:3,j) * ir
        enddo
    
        do mode = 1, max_mode
    
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          do j = 2, nrp + 1
            ir = 1.0 / ( j + noff - 1 )
            cur(1:3,j)  = cur(1:3,j)  * ir; cui(1:3,j)  = cui(1:3,j)  * ir
            dcur(1:2,j) = dcur(1:2,j) * ir; dcui(1:2,j) = dcui(1:2,j) * ir
            amur(1:3,j) = amur(1:3,j) * ir; amui(1:3,j) = amui(1:3,j) * ir
          enddo
        enddo
    
      else
! #ifdef QUADRATIC
!         if (id_proc_loc() == num_procs_loc() - 1) then
!         ! Correct the values on the boundary by adding the value on the guard cell nrp+2
!             cu0(1:3, nrp) = cu0(1:3, nrp) + cu0(1:3, nrp+2)
!             dcu0(1:2, nrp) = dcu0(1:2, nrp) + dcu0(1:2, nrp+2)
!             amu0(1:3, nrp) = amu0(1:3, nrp) + amu0(1:3, nrp+2)
            
!             do mode = 1, max_mode
!                 cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
!                 dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
!                 amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()

!                 cur(1:3, nrp) = cur(1:3, nrp) + cur(1:3, nrp+2);  cui(1:3, nrp) = cui(1:3, nrp) + cui(1:3, nrp+2)
!                 dcur(1:2, nrp) = dcur(1:2, nrp) + dcur(1:2, nrp+2);  dcui(1:2, nrp) = dcui(1:2, nrp) + dcui(1:2, nrp+2)
!                 amur(1:3, nrp) = amur(1:3, nrp) + amur(1:3, nrp+2);  amui(1:3, nrp) = amui(1:3, nrp) + amui(1:3, nrp+2)
!             enddo
!         endif
! #endif
        do j = 0, nrp + 1
           ir = 1.0 / ( j + noff - 1 )
           cu0(1:3,j)  = cu0(1:3,j)  * ir
           dcu0(1:2,j) = dcu0(1:2,j) * ir
           amu0(1:3,j) = amu0(1:3,j) * ir
        enddo
    
        do mode = 1, max_mode
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          do j = 0, nrp+1
            ir = 1.0 / ( j + noff - 1 )
            cur(1:3,j)  = cur(1:3,j)  * ir; cui(1:3,j)  = cui(1:3,j)  * ir
            dcur(1:2,j) = dcur(1:2,j) * ir; dcui(1:2,j) = dcui(1:2,j) * ir
            amur(1:3,j) = amur(1:3,j) * ir; amui(1:3,j) = amui(1:3,j) * ir
          enddo
        enddo
    
      endif
      
      call stop_tprof( 'deposit 2D particles' )
      call write_dbg(cls_name, sname, cls_level, 'ends')
    
end subroutine AMJ_STD


subroutine AMJ_ROB(this, ef, bf, cu, amu, dcu, dt_)
    ! deposit the current, acceleration and momentum flux
    
      implicit none
    
      class(part2d), intent(inout) :: this
      class(field), intent(in) :: cu, amu, dcu
      class(field), intent(in) :: ef, bf
      real, intent(in), optional :: dt_
      ! local data
      character(len=32), save :: sname = 'amjdeposit_robust_part2d'
      type(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
      type(ufield), dimension(:), pointer :: bf_re => null(), bf_im => null()
      type(ufield), dimension(:), pointer :: cu_re => null(), cu_im => null()
      type(ufield), dimension(:), pointer :: dcu_re => null(), dcu_im => null()
      type(ufield), dimension(:), pointer :: amu_re => null(), amu_im => null()
    
      real, dimension(:, :), pointer :: cu0 => null(), dcu0 => null(), amu0 => null()
      real, dimension(:, :), pointer :: cur => null(), dcur => null(), amur => null()
      real, dimension(:, :), pointer :: cui => null(), dcui => null(), amui => null()
    
      integer(kind=LG) :: ptrcur, pp
      integer :: i, j, noff, nrp, np, mode, max_mode
      integer, dimension(p_cache_size) :: ix
      real, dimension(p_p_dim, p_cache_size) :: bp, ep, wp, u0, u, utmp
      real, dimension(:, :), allocatable :: weight
      real, dimension(p_cache_size) :: pcos, psin
      real, dimension(p_p_dim) :: du, u2
      real :: qtmh, qtmh1, qtmh2, idt, gam, ostq, ipsi, dpsi, w, ir, dt
      complex(kind=DB) :: phase, phase0
    
      call write_dbg(cls_name, sname, cls_level, 'starts')
      call start_tprof( 'deposit 2D particles' )
      allocate(weight(LP:UP, p_cache_size))
    
      ef_re  => ef%get_rf_re();  ef_im  => ef%get_rf_im()
      bf_re  => bf%get_rf_re();  bf_im  => bf%get_rf_im()
      cu_re  => cu%get_rf_re();  cu_im  => cu%get_rf_im()
      dcu_re => dcu%get_rf_re(); dcu_im => dcu%get_rf_im()
      amu_re => amu%get_rf_re(); amu_im => amu%get_rf_im()
    
      dt = this%dt
      if (present(dt_)) dt = dt_
    
      idt = 1.0 / dt
      qtmh = 0.5 * this%qbm * dt
      max_mode = ef%get_max_mode()
    
      noff = cu_re(0)%get_noff(1)
      nrp  = cu_re(0)%get_ndp(1)
    
      cu0  => cu_re(0)%get_f1()
      dcu0 => dcu_re(0)%get_f1()
      amu0 => amu_re(0)%get_f1()
    
      do ptrcur = 1, this%npp, p_cache_size
    
        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > this%npp ) then
          np = this%npp - ptrcur + 1
        else
          np = p_cache_size
        endif
    
        ! interpolate fields to particles
        call gen_interp_info(this%x, this%dr, noff, np, ptrcur, weight, ix, pcos, psin)
        call interp_field(ef_re, ef_im, max_mode, weight, ix, pcos, psin, np, ep)
        call interp_field(bf_re, bf_im, max_mode, weight, ix, pcos, psin, np, bp)
    
        ! calculate wake field
        pp = ptrcur
        do i = 1, np
          wp(1, i) = ep(1, i) - bp(2, i)
          wp(2, i) = ep(2, i) + bp(1, i)
          wp(3, i) = ep(3, i)
    
          ! transform momentum from Cartesian to cylindrical coordinates
          u0(1, i) = this%p(1, pp) * pcos(i) + this%p(2, pp) * psin(i)
          u0(2, i) = this%p(2, pp) * pcos(i) - this%p(1, pp) * psin(i)
          u0(3, i) = this%p(3, pp)
    
          ! half electric acceleration
          gam = sqrt(1.0 + u0(1, i)**2 + u0(2, i)**2 + u0(3, i)**2)
          qtmh1 = qtmh * gam / (gam - u0(3, i))
          ep(:, i) = ep(:, i) * qtmh1
          utmp(:, i) = u0(:, i) + ep(:, i)
    
          ! scale magnetic field
          gam = sqrt(1.0 + utmp(1, i)**2 + utmp(2, i)**2 + utmp(3, i)**2)
          qtmh2 = qtmh / (gam - utmp(3, i))
          bp(:, i) = bp(:, i) * qtmh2
    
          ! magnetic rotation
          u(1, i) = utmp(1, i) + utmp(2, i) * bp(3, i) - utmp(3, i) * bp(2, i)
          u(2, i) = utmp(2, i) + utmp(3, i) * bp(1, i) - utmp(1, i) * bp(3, i)
          u(3, i) = utmp(3, i) + utmp(1, i) * bp(2, i) - utmp(2, i) * bp(1, i)
    
          ostq = 2.0 / (1.0 + bp(1, i)**2 + bp(2, i)**2 + bp(3, i)**2)
          bp(:, i) = bp(:, i) * ostq
    
          utmp(1, i) = utmp(1, i) + u(2, i) * bp(3, i) - u(3, i) * bp(2, i)
          utmp(2, i) = utmp(2, i) + u(3, i) * bp(1, i) - u(1, i) * bp(3, i)
          utmp(3, i) = utmp(3, i) + u(1, i) * bp(2, i) - u(2, i) * bp(1, i)
    
          ! half electric acceleration
          u(:, i) = utmp(:, i) + ep(:, i)
    
          pp = pp + 1
        enddo
    
        ! calculate and store time-centered values
        ! deposit momentum flux, acceleration density, and current density
        pp = ptrcur
        do i = 1, np
    
          du(1) = idt * ( u(1,i) - u0(1,i) )
          du(2) = idt * ( u(2,i) - u0(2,i) )
          
          u(:,i)  = 0.5 * ( u(:,i) + u0(:,i) )
          this%gamma(pp) = sqrt( 1.0 + u(1,i)**2 + u(2,i)**2 + u(3,i)**2 )
          ipsi = 1.0 / (this%gamma(pp) - u(3,i))
          this%psi(pp) = (1.0 - 1.0 / ipsi) / this%qbm
    
          dpsi = this%qbm * (wp(3, i) - (wp(1, i) * u(1, i) + wp(2, i) * u(2, i)) * ipsi)
    
          du(1) = du(1) + u(1, i) * dpsi * ipsi
          du(2) = du(2) + u(2, i) * dpsi * ipsi
    
          u2(1) = u(1, i) * u(1, i) * ipsi
          u2(2) = u(1, i) * u(2, i) * ipsi
          u2(3) = u(2, i) * u(2, i) * ipsi
    
          phase0 = cmplx( pcos(i), -psin(i) )
          phase  = cmplx( 1.0, 0.0 ) * this%q(pp) * ipsi
    
          ! deposit m = 0 mode
          do j = LP, UP
            w = weight(j,i) * real(phase)
            cu0( 1:3, ix(i)+j )  = cu0( 1:3, ix(i)+j )  + w * u(1:3,i)
            dcu0( 1:2, ix(i)+j ) = dcu0( 1:2, ix(i)+j ) + w * du(1:2)
            amu0( 1:3, ix(i)+j ) = amu0( 1:3, ix(i)+j ) + w * u2(1:3)
          enddo
    
          ! deposit m > 0 mode
          do mode = 1, max_mode
    
            cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
            dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
            amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
            phase = phase * phase0
    
            do j = LP, UP
              w = weight(j,i) * real(phase)
              cur( 1:3, ix(i)+j )  = cur( 1:3, ix(i)+j )  + w * u(1:3,i)
              dcur( 1:2, ix(i)+j ) = dcur( 1:2, ix(i)+j ) + w * du(1:2)
              amur( 1:3, ix(i)+j ) = amur( 1:3, ix(i)+j ) + w * u2(1:3)
    
              w = weight(j,i) * aimag(phase)
              cui( 1:3, ix(i)+j )  = cui( 1:3, ix(i)+j )  + w * u(1:3,i)
              dcui( 1:2, ix(i)+j ) = dcui( 1:2, ix(i)+j ) + w * du(1:2)
              amui( 1:3, ix(i)+j ) = amui( 1:3, ix(i)+j ) + w * u2(1:3)
            enddo
    
          enddo
    
          pp = pp + 1
        enddo
    
      enddo
    
      if ( id_proc_loc() ==  0 ) then
! Correct the used values on axis by adding the value on the guard cell 0
#ifdef QUADRATIC
        cu0(3,1) = cu0(3,1) + cu0(3,0)
#endif
        ! guard cells on the axis are useless
        cu0(1:3,0)  = 0.0
        dcu0(1:2,0) = 0.0
        amu0(1:3,0) = 0.0
    
        cu0(1:2,1)  = 0.0; cu0(3,1) = 8.0 * cu0(3,1)
        dcu0(1:2,1) = 0.0
        amu0(1:3,1) = 0.0
    
        do mode = 1, max_mode
    
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          ! guard cells on the axis are useless
          cur(1:3,0)  = 0.0; cui(1:3,0)  = 0.0
          dcur(1:2,0) = 0.0; dcui(1:2,0) = 0.0
          amur(1:3,0) = 0.0; amui(1:3,0) = 0.0
           
          if ( mode == 1 ) then
#ifdef QUADRATIC
            cur(1:2,1) = cur(1:2,1) + cur(1:2,0)
            dcur(1:2,1) = dcur(1:2,1) + dcur(1:2,0)
            cui(1:2,1) = cui(1:2,1) + cui(1:2,0)
            dcui(1:2,1) = dcui(1:2,1) + dcui(1:2,0)
  #endif
            cur(1:2,1)  = 8.0 * cur(1:2,1); cur(3,1) = 0.0
            dcur(1:2,1) = 8.0 * dcur(1:2,1)
            amur(1:3,1) = 0.0
            cui(1:2,1)  = 8.0 * cui(1:2,1); cui(3,1) = 0.0
            dcui(1:2,1) = 8.0 * dcui(1:2,1)
            amui(1:3,1) = 0.0
          elseif ( mode == 2 ) then
#ifdef QUADRATIC
            amur(1:3,1) = amur(1:3,1) + amur(1:3,0)
            amui(1:3,1) = amui(1:3,1) + amui(1:3,0)
#endif
            cur(1:3,1)  = 0.0
            dcur(1:2,1) = 0.0
            amur(1:3,1) = 8.0 * amur(1:3,1)
            cui(1:3,1)  = 0.0
            dcui(1:2,1) = 0.0
            amui(1:3,1) = 8.0 * amui(1:3,1)
          else
            cur(1:3,1)  = 0.0
            dcur(1:2,1) = 0.0
            amur(1:3,1) = 0.0
            cui(1:3,1)  = 0.0
            dcui(1:2,1) = 0.0
            amui(1:3,1) = 0.0
           endif
        enddo
    
        do j = 2, nrp + 1
          ir = 1.0 / ( j + noff - 1 )
          cu0(1:3,j)  = cu0(1:3,j)  * ir
          dcu0(1:2,j) = dcu0(1:2,j) * ir
          amu0(1:3,j) = amu0(1:3,j) * ir
        enddo
    
        do mode = 1, max_mode
    
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          do j = 2, nrp + 1
            ir = 1.0 / ( j + noff - 1 )
            cur(1:3,j)  = cur(1:3,j)  * ir; cui(1:3,j)  = cui(1:3,j)  * ir
            dcur(1:2,j) = dcur(1:2,j) * ir; dcui(1:2,j) = dcui(1:2,j) * ir
            amur(1:3,j) = amur(1:3,j) * ir; amui(1:3,j) = amui(1:3,j) * ir
          enddo
        enddo
    
      else
! #ifdef QUADRATIC
!         if (id_proc_loc() == num_procs_loc() - 1) then
!         ! Correct the values on the boundary by adding the value on the guard cell nrp+2
!             cu0(1:3, nrp) = cu0(1:3, nrp) + cu0(1:3, nrp+2)
!             dcu0(1:2, nrp) = dcu0(1:2, nrp) + dcu0(1:2, nrp+2)
!             amu0(1:3, nrp) = amu0(1:3, nrp) + amu0(1:3, nrp+2)
            
!             do mode = 1, max_mode
!                 cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
!                 dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
!                 amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()

!                 cur(1:3, nrp) = cur(1:3, nrp) + cur(1:3, nrp+2);  cui(1:3, nrp) = cui(1:3, nrp) + cui(1:3, nrp+2)
!                 dcur(1:2, nrp) = dcur(1:2, nrp) + dcur(1:2, nrp+2);  dcui(1:2, nrp) = dcui(1:2, nrp) + dcui(1:2, nrp+2)
!                 amur(1:3, nrp) = amur(1:3, nrp) + amur(1:3, nrp+2);  amui(1:3, nrp) = amui(1:3, nrp) + amui(1:3, nrp+2)
!             enddo
!         endif
! #endif
        do j = 0, nrp + 1
           ir = 1.0 / ( j + noff - 1 )
           cu0(1:3,j)  = cu0(1:3,j)  * ir
           dcu0(1:2,j) = dcu0(1:2,j) * ir
           amu0(1:3,j) = amu0(1:3,j) * ir
        enddo
    
        do mode = 1, max_mode
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          do j = 0, nrp+1
            ir = 1.0 / ( j + noff - 1 )
            cur(1:3,j)  = cur(1:3,j)  * ir; cui(1:3,j)  = cui(1:3,j)  * ir
            dcur(1:2,j) = dcur(1:2,j) * ir; dcui(1:2,j) = dcui(1:2,j) * ir
            amur(1:3,j) = amur(1:3,j) * ir; amui(1:3,j) = amui(1:3,j) * ir
          enddo
        enddo
    
      endif
      
      call stop_tprof( 'deposit 2D particles' )
      call write_dbg(cls_name, sname, cls_level, 'ends')
    
end subroutine AMJ_ROB

subroutine AMJ_STD_PGC(this, ef, bf, af, cu, amu, dcu, dt_)
    ! deposit the current, acceleration and momentum flux
    
      implicit none
    
      class(part2d), intent(inout) :: this
      class(field), intent(in) :: cu, amu, dcu
      class(field), intent(in) :: ef, bf
      class(field_laser), intent(in) :: af
      real, intent(in), optional :: dt_
      ! local data
      character(len=32), save :: sname = 'amjdeposit_std_pgc_part2d'
      type(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
      type(ufield), dimension(:), pointer :: bf_re => null(), bf_im => null()
      type(ufield), dimension(:), pointer :: cu_re => null(), cu_im => null()
      type(ufield), dimension(:), pointer :: dcu_re => null(), dcu_im => null()
      type(ufield), dimension(:), pointer :: amu_re => null(), amu_im => null()
      type(ufield), dimension(:), pointer :: ar_re, ar_im, ai_re, ai_im
      type(ufield), dimension(:), pointer :: ar_grad_re, ar_grad_im, ai_grad_re, ai_grad_im
    
      real, dimension(:,:), pointer :: cu0 => null(), dcu0 => null(), amu0 => null()
      real, dimension(:,:), pointer :: cur => null(), dcur => null(), amur => null()
      real, dimension(:,:), pointer :: cui => null(), dcui => null(), amui => null()
    
      integer(kind=LG) :: ptrcur, pp
      integer :: i, j, noff, nrp, np, mode, max_mode
      integer, dimension(p_cache_size) :: ix
      real, dimension(p_p_dim, p_cache_size) :: bp, ep, wp, u0, u
      real, dimension(p_cache_size) :: apr, api
      real, dimension(3,p_cache_size) :: apr_grad, api_grad
      real, dimension(:, :), allocatable :: weight
      real, dimension(p_cache_size) :: pcos, psin
      real, dimension(p_p_dim) :: du, u2, utmp
      real :: qtmh, qtmh_e, qtmh_b, idt, gam, ostq, ipsi, dpsi, w, ir, gam_corr, tmp, dt
      complex(kind=DB) :: phase, phase0
    
      call write_dbg(cls_name, sname, cls_level, 'starts')
      call start_tprof( 'deposit 2D particles' )
      allocate(weight(LP:UP, p_cache_size))
    
      ef_re  => ef%get_rf_re();  ef_im  => ef%get_rf_im()
      bf_re  => bf%get_rf_re();  bf_im  => bf%get_rf_im()
      cu_re  => cu%get_rf_re();  cu_im  => cu%get_rf_im()
      dcu_re => dcu%get_rf_re(); dcu_im => dcu%get_rf_im()
      amu_re => amu%get_rf_re(); amu_im => amu%get_rf_im()
    
      dt = this%dt
      if (present(dt_)) dt = dt_
    
      idt = 1.0 / dt
      qtmh = 0.5 * this%qbm * dt
      max_mode = ef%get_max_mode()
    
      noff = cu_re(0)%get_noff(1)
      nrp  = cu_re(0)%get_ndp(1)
    
      cu0  => cu_re(0)%get_f1()
      dcu0 => dcu_re(0)%get_f1()
      amu0 => amu_re(0)%get_f1()
    
      ar_re => af%get_cfr_re()
      ar_im => af%get_cfr_im()
      ai_re => af%get_cfi_re()
      ai_im => af%get_cfi_im()
    
      ar_grad_re => af%ar_grad_re
      ar_grad_im => af%ar_grad_im
      ai_grad_re => af%ai_grad_re
      ai_grad_im => af%ai_grad_im
    
      do ptrcur = 1, this%npp, p_cache_size
    
        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > this%npp ) then
          np = this%npp - ptrcur + 1
        else
          np = p_cache_size
        endif
    
        ! interpolate fields to particles
        ! call interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, this%x, this%dr, &
        !   bp, ep, np, ptrcur, p_cylindrical, weight = wt, ix = ix, pcos = cc, psin = ss )
        ! call interp_laser_part2d( ar_re, ar_im, ai_re, ai_im, &
        !   ar_grad_re, ar_grad_im, ai_grad_re, ai_grad_im, &
        !   max_mode, this%x, this%dr, apr, api, apr_grad, api_grad, np, ptrcur, p_cylindrical )
        call gen_interp_info(this%x, this%dr, noff, np, ptrcur, weight, ix, pcos, psin)
        call interp_field(ef_re, ef_im, max_mode, weight, ix, pcos, psin, np, ep)
        call interp_field(bf_re, bf_im, max_mode, weight, ix, pcos, psin, np, bp)
        call interp_field(ar_re, ar_im, max_mode, weight, ix, pcos, psin, np, apr)
        call interp_field(ai_re, ai_im, max_mode, weight, ix, pcos, psin, np, api)
        call interp_field(ar_grad_re, ar_grad_im, max_mode, weight, ix, pcos, psin, np, apr_grad)
        call interp_field(ai_grad_re, ai_grad_im, max_mode, weight, ix, pcos, psin, np, api_grad)
    
        ! advance the particle momentum
        pp = ptrcur
        do i = 1, np
    
          ! calculate the correction factor for gamma due to ponderomotive force
          gam_corr = 0.5 * this%qbm**2 * (apr(i)**2 + api(i)**2)
    
          ! transform momentum from Cartesian to cylindrical coordinates
          u0(1, i) = this%p(1, pp) * pcos(i) + this%p(2, pp) * psin(i)
          u0(2, i) = this%p(2, pp) * pcos(i) - this%p(1, pp) * psin(i)
          u0(3, i) = this%p(3, pp)
    
          ! calculate the averaged gamma factor
          gam = sqrt(1.0 + u0(1, i)**2 + u0(2, i)**2 + u0(3, i)**2 + gam_corr)
    
          ! calculate wake field
          wp(1, i) = ep(1, i) - bp(2, i)
          wp(2, i) = ep(2, i) + bp(1, i)
          wp(3, i) = ep(3, i)
    
          ! TODO: There is a strict time-centered scheme for the standard-PGC pusher, which is left for future work.
    
          ! calculate the effective electric fields
          ! note that here the time-centered gamma is unknown, so we approximate it using the initial one.
          tmp = 0.5 * this%qbm / gam
          ep(1, i) = ep(1, i) - tmp * (apr(i) * apr_grad(1, i) + api(i) * api_grad(1, i))
          ep(2, i) = ep(2, i) - tmp * (apr(i) * apr_grad(2, i) + api(i) * api_grad(2, i))
          ep(3, i) = ep(3, i) + tmp * (apr(i) * apr_grad(3, i) + api(i) * api_grad(3, i))
    
          ! half electric acceleration
          qtmh_b = qtmh / (1.0 - this%qbm * this%psi(pp))
          qtmh_e = qtmh_b * gam
          utmp(:) = u0(:, i) + ep(:, i) * qtmh_e
    
          ! magnetic rotation
          bp(:, i) = bp(:, i) * qtmh_b
          u(1, i) = utmp(1) + utmp(2) * bp(3, i) - utmp(3) * bp(2, i)
          u(2, i) = utmp(2) + utmp(3) * bp(1, i) - utmp(1) * bp(3, i)
          u(3, i) = utmp(3) + utmp(1) * bp(2, i) - utmp(2) * bp(1, i)
    
          ostq = 2.0 / (1.0 + bp(1, i)**2 + bp(2, i)**2 + bp(3, i)**2)
          bp(:, i) = bp(:, i) * ostq
    
          utmp(1) = utmp(1) + u(2, i) * bp(3, i) - u(3, i) * bp(2, i)
          utmp(2) = utmp(2) + u(3, i) * bp(1, i) - u(1, i) * bp(3, i)
          utmp(3) = utmp(3) + u(1, i) * bp(2, i) - u(2, i) * bp(1, i)
    
          ! half electric acceleration
          gam = sqrt(1.0 + utmp(1)**2 + utmp(2)**2 + utmp(3)**2 + gam_corr)
          qtmh_e = qtmh_b * gam
          u(:, i) = utmp(:) + ep(:, i) * qtmh_e
    
          ! calculate and store time-centered values
          ! deposit momentum flux, acceleration density, and current density
          du(1) = idt * (u(1, i) - u0(1, i))
          du(2) = idt * (u(2, i) - u0(2, i))
          
          u(:, i)  = 0.5 * (u(:, i) + u0(:, i))
          ! store the time-centered gamma
          this%gamma(pp) = sqrt(1.0 + u(1, i)**2 + u(2, i)**2 + u(3, i)**2 + gam_corr)
    
          ipsi = 1.0 / (1.0 - this%qbm * this%psi(pp))
          dpsi = this%qbm * (wp(3, i) - (wp(1, i) * u(1, i) + wp(2, i) * u(2, i)) * ipsi)
    
          du(1) = du(1) + u(1, i) * dpsi * ipsi
          du(2) = du(2) + u(2, i) * dpsi * ipsi
    
          u2(1) = u(1, i) * u(1, i) * ipsi
          u2(2) = u(1, i) * u(2, i) * ipsi
          u2(3) = u(2, i) * u(2, i) * ipsi
    
          phase0 = cmplx(pcos(i), -psin(i))
          phase  = cmplx(1.0, 0.0) * this%q(pp) * ipsi
    
          ! deposit m = 0 mode
          do j = LP, UP
            w = weight(j, i) * real(phase)
            cu0(1:3, ix(i)+j)  = cu0(1:3, ix(i) + j)  + w * u(1:3, i)
            dcu0(1:2, ix(i)+j) = dcu0(1:2, ix(i) + j) + w * du(1:2)
            amu0(1:3, ix(i)+j) = amu0(1:3, ix(i) + j) + w * u2(1:3)
          enddo
    
          ! deposit m > 0 mode
          do mode = 1, max_mode
    
            cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
            dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
            amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
            phase = phase * phase0
    
            do j = LP, UP
              w = weight(j,i) * real(phase)
              cur( 1:3, ix(i)+j )  = cur( 1:3, ix(i)+j )  + w * u(1:3,i)
              dcur( 1:2, ix(i)+j ) = dcur( 1:2, ix(i)+j ) + w * du(1:2)
              amur( 1:3, ix(i)+j ) = amur( 1:3, ix(i)+j ) + w * u2(1:3)
    
              w = weight(j,i) * aimag(phase)
              cui( 1:3, ix(i)+j )  = cui( 1:3, ix(i)+j )  + w * u(1:3,i)
              dcui( 1:2, ix(i)+j ) = dcui( 1:2, ix(i)+j ) + w * du(1:2)
              amui( 1:3, ix(i)+j ) = amui( 1:3, ix(i)+j ) + w * u2(1:3)
            enddo
    
          enddo
    
          pp = pp + 1
        enddo
    
      enddo
    
      if ( id_proc_loc() ==  0 ) then
! Correct the used values on axis by adding the value on the guard cell 0
#ifdef QUADRATIC
        cu0(3,1) = cu0(3,1) + cu0(3,0)
#endif
        ! guard cells on the axis are useless
        cu0(1:3,0)  = 0.0
        dcu0(1:2,0) = 0.0
        amu0(1:3,0) = 0.0
    
        cu0(1:2,1)  = 0.0; cu0(3,1) = 8.0 * cu0(3,1)
        dcu0(1:2,1) = 0.0
        amu0(1:3,1) = 0.0
    
        do mode = 1, max_mode
    
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          ! guard cells on the axis are useless
          cur(1:3,0)  = 0.0; cui(1:3,0)  = 0.0
          dcur(1:2,0) = 0.0; dcui(1:2,0) = 0.0
          amur(1:3,0) = 0.0; amui(1:3,0) = 0.0
            
          if ( mode == 1 ) then
#ifdef QUADRATIC
            cur(1:2,1) = cur(1:2,1) + cur(1:2,0)
            dcur(1:2,1) = dcur(1:2,1) + dcur(1:2,0)
            cui(1:2,1) = cui(1:2,1) + cui(1:2,0)
            dcui(1:2,1) = dcui(1:2,1) + dcui(1:2,0)
#endif
            cur(1:2,1)  = 8.0 * cur(1:2,1); cur(3,1) = 0.0
            dcur(1:2,1) = 8.0 * dcur(1:2,1)
            amur(1:3,1) = 0.0
            cui(1:2,1)  = 8.0 * cui(1:2,1); cui(3,1) = 0.0
            dcui(1:2,1) = 8.0 * dcui(1:2,1)
            amui(1:3,1) = 0.0
          elseif ( mode == 2 ) then
#ifdef QUADRATIC
            amur(1:3,1) = amur(1:3,1) + amur(1:3,0)
            amui(1:3,1) = amui(1:3,1) + amui(1:3,0)
#endif
            cur(1:3,1)  = 0.0
            dcur(1:2,1) = 0.0
            amur(1:3,1) = 8.0 * amur(1:3,1)
            cui(1:3,1)  = 0.0
            dcui(1:2,1) = 0.0
            amui(1:3,1) = 8.0 * amui(1:3,1)
          else
            cur(1:3,1)  = 0.0
            dcur(1:2,1) = 0.0
            amur(1:3,1) = 0.0
            cui(1:3,1)  = 0.0
            dcui(1:2,1) = 0.0
            amui(1:3,1) = 0.0
            endif
        enddo
    
        do j = 2, nrp + 1
          ir = 1.0 / ( j + noff - 1 )
          cu0(1:3,j)  = cu0(1:3,j)  * ir
          dcu0(1:2,j) = dcu0(1:2,j) * ir
          amu0(1:3,j) = amu0(1:3,j) * ir
        enddo
    
        do mode = 1, max_mode
    
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          do j = 2, nrp + 1
            ir = 1.0 / ( j + noff - 1 )
            cur(1:3,j)  = cur(1:3,j)  * ir; cui(1:3,j)  = cui(1:3,j)  * ir
            dcur(1:2,j) = dcur(1:2,j) * ir; dcui(1:2,j) = dcui(1:2,j) * ir
            amur(1:3,j) = amur(1:3,j) * ir; amui(1:3,j) = amui(1:3,j) * ir
          enddo
        enddo
    
      else
! #ifdef QUADRATIC
!         if (id_proc_loc() == num_procs_loc() - 1) then
!         ! Correct the values on the boundary by adding the value on the guard cell nrp+2
!             cu0(1:3, nrp) = cu0(1:3, nrp) + cu0(1:3, nrp+2)
!             dcu0(1:2, nrp) = dcu0(1:2, nrp) + dcu0(1:2, nrp+2)
!             amu0(1:3, nrp) = amu0(1:3, nrp) + amu0(1:3, nrp+2)
            
!             do mode = 1, max_mode
!                 cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
!                 dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
!                 amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()

!                 cur(1:3, nrp) = cur(1:3, nrp) + cur(1:3, nrp+2);  cui(1:3, nrp) = cui(1:3, nrp) + cui(1:3, nrp+2)
!                 dcur(1:2, nrp) = dcur(1:2, nrp) + dcur(1:2, nrp+2);  dcui(1:2, nrp) = dcui(1:2, nrp) + dcui(1:2, nrp+2)
!                 amur(1:3, nrp) = amur(1:3, nrp) + amur(1:3, nrp+2);  amui(1:3, nrp) = amui(1:3, nrp) + amui(1:3, nrp+2)
!             enddo
!         endif
! #endif
        do j = 0, nrp + 1
            ir = 1.0 / ( j + noff - 1 )
            cu0(1:3,j)  = cu0(1:3,j)  * ir
            dcu0(1:2,j) = dcu0(1:2,j) * ir
            amu0(1:3,j) = amu0(1:3,j) * ir
        enddo
    
        do mode = 1, max_mode
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          do j = 0, nrp+1
            ir = 1.0 / ( j + noff - 1 )
            cur(1:3,j)  = cur(1:3,j)  * ir; cui(1:3,j)  = cui(1:3,j)  * ir
            dcur(1:2,j) = dcur(1:2,j) * ir; dcui(1:2,j) = dcui(1:2,j) * ir
            amur(1:3,j) = amur(1:3,j) * ir; amui(1:3,j) = amui(1:3,j) * ir
          enddo
        enddo
    
      endif
      
      call stop_tprof( 'deposit 2D particles' )
      call write_dbg(cls_name, sname, cls_level, 'ends')
      
end subroutine AMJ_STD_PGC

subroutine AMJ_ROB_PGC(this, ef, bf, af, cu, amu, dcu, dt_)
    ! deposit the current, acceleration and momentum flux
    
      implicit none
    
      class(part2d), intent(inout) :: this
      class(field), intent(in) :: cu, amu, dcu
      class(field), intent(in) :: ef, bf
      class(field_laser), intent(in) :: af
      real, intent(in), optional :: dt_
      ! local data
      character(len=32), save :: sname = 'amjdeposit_robust_pgc_part2d'
      type(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
      type(ufield), dimension(:), pointer :: bf_re => null(), bf_im => null()
      type(ufield), dimension(:), pointer :: cu_re => null(), cu_im => null()
      type(ufield), dimension(:), pointer :: dcu_re => null(), dcu_im => null()
      type(ufield), dimension(:), pointer :: amu_re => null(), amu_im => null()
      type(ufield), dimension(:), pointer :: ar_re, ar_im, ai_re, ai_im
      type(ufield), dimension(:), pointer :: ar_grad_re, ar_grad_im, ai_grad_re, ai_grad_im
    
      real, dimension(:,:), pointer :: cu0 => null(), dcu0 => null(), amu0 => null()
      real, dimension(:,:), pointer :: cur => null(), dcur => null(), amur => null()
      real, dimension(:,:), pointer :: cui => null(), dcui => null(), amui => null()
    
      integer(kind=LG) :: ptrcur, pp
      integer :: i, j, noff, nrp, np, mode, max_mode
      integer, dimension(p_cache_size) :: ix
      real, dimension(p_p_dim, p_cache_size) :: bp, ep, wp, u0, u
      real, dimension(p_cache_size) :: apr, api
      real, dimension(3,p_cache_size) :: apr_grad, api_grad
      real, dimension(:, :), allocatable :: weight
      real, dimension(p_cache_size) :: pcos, psin
      real, dimension(p_p_dim) :: du, u2, utmp
      real :: qtmh, qtmh_e, qtmh_b, idt, gam, ostq, ipsi, dpsi, w, ir, gam_corr, tmp, dt
      complex(kind=DB) :: phase, phase0
    
      call write_dbg(cls_name, sname, cls_level, 'starts')
      call start_tprof( 'deposit 2D particles' )
      
      allocate(weight(LP:UP, p_cache_size))
      
      ef_re  => ef%get_rf_re();  ef_im  => ef%get_rf_im()
      bf_re  => bf%get_rf_re();  bf_im  => bf%get_rf_im()
      cu_re  => cu%get_rf_re();  cu_im  => cu%get_rf_im()
      dcu_re => dcu%get_rf_re(); dcu_im => dcu%get_rf_im()
      amu_re => amu%get_rf_re(); amu_im => amu%get_rf_im()
    
      dt = this%dt
      if (present(dt_)) dt = dt_
    
      idt = 1.0 / dt
      qtmh = 0.5 * this%qbm * dt
      max_mode = ef%get_max_mode()
    
      noff = cu_re(0)%get_noff(1)
      nrp  = cu_re(0)%get_ndp(1)
    
      cu0  => cu_re(0)%get_f1()
      dcu0 => dcu_re(0)%get_f1()
      amu0 => amu_re(0)%get_f1()
    
      ar_re => af%get_cfr_re()
      ar_im => af%get_cfr_im()
      ai_re => af%get_cfi_re()
      ai_im => af%get_cfi_im()
    
      ar_grad_re => af%ar_grad_re
      ar_grad_im => af%ar_grad_im
      ai_grad_re => af%ai_grad_re
      ai_grad_im => af%ai_grad_im
    
      do ptrcur = 1, this%npp, p_cache_size
    
        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > this%npp ) then
          np = this%npp - ptrcur + 1
        else
          np = p_cache_size
        endif
    
        ! interpolate fields to particles
        call gen_interp_info(this%x, this%dr, noff, np, ptrcur, weight, ix, pcos, psin)
        call interp_field(ef_re, ef_im, max_mode, weight, ix, pcos, psin, np, ep)
        call interp_field(bf_re, bf_im, max_mode, weight, ix, pcos, psin, np, bp)
        call interp_field(ar_re, ar_im, max_mode, weight, ix, pcos, psin, np, apr)
        call interp_field(ai_re, ai_im, max_mode, weight, ix, pcos, psin, np, api)
        call interp_field(ar_grad_re, ar_grad_im, max_mode, weight, ix, pcos, psin, np, apr_grad)
        call interp_field(ai_grad_re, ai_grad_im, max_mode, weight, ix, pcos, psin, np, api_grad)
    
        ! advance the particle momentum
        pp = ptrcur
        do i = 1, np
    
          ! calculate the correction factor for gamma due to ponderomotive force
          gam_corr = 0.5 * this%qbm**2 * ( apr(i)**2 + api(i)**2 )
    
          ! transform momentum from Cartesian to cylindrical coordinates
          u0(1,i) = this%p(1,pp) * pcos(i) + this%p(2,pp) * psin(i)
          u0(2,i) = this%p(2,pp) * pcos(i) - this%p(1,pp) * psin(i)
          u0(3,i) = this%p(3,pp)
    
          ! calculate the averaged gamma factor
          gam = sqrt( 1.0 + u0(1,i)**2 + u0(2,i)**2 + u0(3,i)**2 + gam_corr )
    
          ! calculate wake field
          wp(1,i) = ep(1,i) - bp(2,i)
          wp(2,i) = ep(2,i) + bp(1,i)
          wp(3,i) = ep(3,i)
    
          ! calculate the effective electric fields
          ! note that here the time-centered gamma is unknown, so we approximate it using the initial one.
          tmp = 0.5 * this%qbm / gam
          ep(1,i) = ep(1,i) - tmp * ( apr(i) * apr_grad(1,i) + api(i) * api_grad(1,i) )
          ep(2,i) = ep(2,i) - tmp * ( apr(i) * apr_grad(2,i) + api(i) * api_grad(2,i) )
          ep(3,i) = ep(3,i) + tmp * ( apr(i) * apr_grad(3,i) + api(i) * api_grad(3,i) )
    
          ! half electric acceleration
          qtmh_e = qtmh * gam / ( gam - u0(3,i) )
          utmp(:) = u0(:,i) + ep(:,i) * qtmh_e
    
          ! magnetic rotation
          gam = sqrt( 1.0 + utmp(1)**2 + utmp(2)**2 + utmp(3)**2 + gam_corr )
          qtmh_b = qtmh / ( gam - utmp(3) )
          bp(:,i) = bp(:,i) * qtmh_b
    
          u(1,i) = utmp(1) + utmp(2) * bp(3,i) - utmp(3) * bp(2,i)
          u(2,i) = utmp(2) + utmp(3) * bp(1,i) - utmp(1) * bp(3,i)
          u(3,i) = utmp(3) + utmp(1) * bp(2,i) - utmp(2) * bp(1,i)
    
          ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
          bp(:,i) = bp(:,i) * ostq
    
          utmp(1) = utmp(1) + u(2,i) * bp(3,i) - u(3,i) * bp(2,i)
          utmp(2) = utmp(2) + u(3,i) * bp(1,i) - u(1,i) * bp(3,i)
          utmp(3) = utmp(3) + u(1,i) * bp(2,i) - u(2,i) * bp(1,i)
    
          ! half electric acceleration
          gam = sqrt( 1.0 + utmp(1)**2 + utmp(2)**2 + utmp(3)**2 + gam_corr )
          qtmh_e = qtmh * gam / ( gam - utmp(3) )
          u(:,i) = utmp(:) + ep(:,i) * qtmh_e
    
          ! calculate and store time-centered values
          ! deposit momentum flux, acceleration density, and current density
          du(1) = idt * ( u(1,i) - u0(1,i) )
          du(2) = idt * ( u(2,i) - u0(2,i) )
          
          u(:,i)  = 0.5 * ( u(:,i) + u0(:,i) )
          ! store the time-centered gamma and (1-(q/m)*psi)
          this%gamma(pp) = sqrt( 1.0 + u(1,i)**2 + u(2,i)**2 + u(3,i)**2 + gam_corr )
          ipsi = 1.0 / (this%gamma(pp) - u(3,i))
          this%psi(pp) = (1.0 - 1.0 / ipsi) / this%qbm
    
          dpsi = this%qbm * ( wp(3,i) - ( wp(1,i) * u(1,i) + wp(2,i) * u(2,i) ) * ipsi )
    
          du(1) = du(1) + u(1,i) * dpsi * ipsi
          du(2) = du(2) + u(2,i) * dpsi * ipsi
    
          u2(1) = u(1,i) * u(1,i) * ipsi
          u2(2) = u(1,i) * u(2,i) * ipsi
          u2(3) = u(2,i) * u(2,i) * ipsi
    
          phase0 = cmplx( pcos(i), -psin(i) )
          phase  = cmplx( 1.0, 0.0 ) * this%q(pp) * ipsi
    
          ! deposit m = 0 mode
          do j = LP, UP
            w = weight(j,i) * real(phase)
            cu0( 1:3, ix(i)+j )  = cu0( 1:3, ix(i)+j )  + w * u(1:3,i)
            dcu0( 1:2, ix(i)+j ) = dcu0( 1:2, ix(i)+j ) + w * du(1:2)
            amu0( 1:3, ix(i)+j ) = amu0( 1:3, ix(i)+j ) + w * u2(1:3)
          enddo
    
          ! deposit m > 0 mode
          do mode = 1, max_mode
    
            cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
            dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
            amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
            phase = phase * phase0
    
            do j = LP, UP
              w = weight(j,i) * real(phase)
              cur( 1:3, ix(i)+j )  = cur( 1:3, ix(i)+j )  + w * u(1:3,i)
              dcur( 1:2, ix(i)+j ) = dcur( 1:2, ix(i)+j ) + w * du(1:2)
              amur( 1:3, ix(i)+j ) = amur( 1:3, ix(i)+j ) + w * u2(1:3)
    
              w = weight(j,i) * aimag(phase)
              cui( 1:3, ix(i)+j )  = cui( 1:3, ix(i)+j )  + w * u(1:3,i)
              dcui( 1:2, ix(i)+j ) = dcui( 1:2, ix(i)+j ) + w * du(1:2)
              amui( 1:3, ix(i)+j ) = amui( 1:3, ix(i)+j ) + w * u2(1:3)
            enddo
    
          enddo
    
          pp = pp + 1
        enddo
    
      enddo
    
      if ( id_proc_loc() == 0 ) then
! Correct the used values on axis by adding the value on the guard cell 0
#ifdef QUADRATIC
    cu0(3,1) = cu0(3,1) + cu0(3,0)
#endif
        ! guard cells on the axis are useless
        cu0(1:3,0)  = 0.0
        dcu0(1:2,0) = 0.0
        amu0(1:3,0) = 0.0
    
        cu0(1:2,1)  = 0.0; cu0(3,1) = 8.0 * cu0(3,1)
        dcu0(1:2,1) = 0.0
        amu0(1:3,1) = 0.0
    
        do mode = 1, max_mode
    
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          ! guard cells on the axis are useless
          cur(1:3,0)  = 0.0; cui(1:3,0)  = 0.0
          dcur(1:2,0) = 0.0; dcui(1:2,0) = 0.0
          amur(1:3,0) = 0.0; amui(1:3,0) = 0.0
            
          if ( mode == 1 ) then
#ifdef QUADRATIC
            cur(1:2,1) = cur(1:2,1) + cur(1:2,0)
            dcur(1:2,1) = dcur(1:2,1) + dcur(1:2,0)
            cui(1:2,1) = cui(1:2,1) + cui(1:2,0)
            dcui(1:2,1) = dcui(1:2,1) + dcui(1:2,0)
#endif
            cur(1:2,1)  = 8.0 * cur(1:2,1); cur(3,1) = 0.0
            dcur(1:2,1) = 8.0 * dcur(1:2,1)
            amur(1:3,1) = 0.0
            cui(1:2,1)  = 8.0 * cui(1:2,1); cui(3,1) = 0.0
            dcui(1:2,1) = 8.0 * dcui(1:2,1)
            amui(1:3,1) = 0.0
          elseif ( mode == 2 ) then
#ifdef QUADRATIC
            amur(1:3,1) = amur(1:3,1) + amur(1:3,0)
            amui(1:3,1) = amui(1:3,1) + amui(1:3,0)
#endif
            cur(1:3,1)  = 0.0
            dcur(1:2,1) = 0.0
            amur(1:3,1) = 8.0 * amur(1:3,1)
            cui(1:3,1)  = 0.0
            dcui(1:2,1) = 0.0
            amui(1:3,1) = 8.0 * amui(1:3,1)
          else
            cur(1:3,1)  = 0.0
            dcur(1:2,1) = 0.0
            amur(1:3,1) = 0.0
            cui(1:3,1)  = 0.0
            dcui(1:2,1) = 0.0
            amui(1:3,1) = 0.0
            endif
        enddo
    
        do j = 2, nrp + 1
          ir = 1.0 / ( j + noff - 1 )
          cu0(1:3,j)  = cu0(1:3,j)  * ir
          dcu0(1:2,j) = dcu0(1:2,j) * ir
          amu0(1:3,j) = amu0(1:3,j) * ir
        enddo
    
        do mode = 1, max_mode
    
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          do j = 2, nrp + 1
            ir = 1.0 / ( j + noff - 1 )
            cur(1:3,j)  = cur(1:3,j)  * ir; cui(1:3,j)  = cui(1:3,j)  * ir
            dcur(1:2,j) = dcur(1:2,j) * ir; dcui(1:2,j) = dcui(1:2,j) * ir
            amur(1:3,j) = amur(1:3,j) * ir; amui(1:3,j) = amui(1:3,j) * ir
          enddo
        enddo
    
      else
! #ifdef QUADRATIC
!         if (id_proc_loc() == num_procs_loc() - 1) then
!         ! Correct the values on the boundary by adding the value on the guard cell nrp+2
!             cu0(1:3, nrp) = cu0(1:3, nrp) + cu0(1:3, nrp+2)
!             dcu0(1:2, nrp) = dcu0(1:2, nrp) + dcu0(1:2, nrp+2)
!             amu0(1:3, nrp) = amu0(1:3, nrp) + amu0(1:3, nrp+2)
            
!             do mode = 1, max_mode
!                 cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
!                 dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
!                 amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()

!                 cur(1:3, nrp) = cur(1:3, nrp) + cur(1:3, nrp+2);  cui(1:3, nrp) = cui(1:3, nrp) + cui(1:3, nrp+2)
!                 dcur(1:2, nrp) = dcur(1:2, nrp) + dcur(1:2, nrp+2);  dcui(1:2, nrp) = dcui(1:2, nrp) + dcui(1:2, nrp+2)
!                 amur(1:3, nrp) = amur(1:3, nrp) + amur(1:3, nrp+2);  amui(1:3, nrp) = amui(1:3, nrp) + amui(1:3, nrp+2)
!             enddo
!         endif
! #endif
        do j = 0, nrp + 1
            ir = 1.0 / ( j + noff - 1 )
            cu0(1:3,j)  = cu0(1:3,j)  * ir
            dcu0(1:2,j) = dcu0(1:2,j) * ir
            amu0(1:3,j) = amu0(1:3,j) * ir
        enddo
    
        do mode = 1, max_mode
          cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
          dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
          amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()
    
          do j = 0, nrp+1
            ir = 1.0 / ( j + noff - 1 )
            cur(1:3,j)  = cur(1:3,j)  * ir; cui(1:3,j)  = cui(1:3,j)  * ir
            dcur(1:2,j) = dcur(1:2,j) * ir; dcui(1:2,j) = dcui(1:2,j) * ir
            amur(1:3,j) = amur(1:3,j) * ir; amui(1:3,j) = amui(1:3,j) * ir
          enddo
        enddo
    
      endif
      
      call stop_tprof( 'deposit 2D particles' )
      call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine AMJ_ROB_PGC


subroutine PUSH_STD(this, ef, bf, dt_)

    implicit none
  
    class(part2d), intent(inout) :: this
    class(field), intent(in) :: ef, bf
    real, intent(in), optional :: dt_
    ! local data
    character(len=18), save :: sname = 'push_u_std_part2d'
    type(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im
    integer :: i, np, max_mode, noff
    real :: qtmh, qtmh1, qtmh2, gam, dtc, ostq, dt
    real, dimension(:, :), allocatable :: weight
    real, dimension(p_cache_size) :: pcos, psin
    integer, dimension(p_cache_size) :: idx
    real, dimension(p_p_dim, p_cache_size) :: bp, ep, utmp
    integer(kind=LG) :: ptrcur, pp
  
    call write_dbg(cls_name, sname, cls_level, 'starts')
    call start_tprof( 'push 2D particles' )
    allocate(weight(LP:UP, p_cache_size))
  
    dt = this%dt
    if (present(dt_)) dt = dt_
    qtmh = this%qbm * dt * 0.5
    max_mode = ef%get_max_mode()
  
    ef_re => ef%get_rf_re()
    ef_im => ef%get_rf_im()
    bf_re => bf%get_rf_re()
    bf_im => bf%get_rf_im()
    noff = ef_re(0)%get_noff(1)
  
    do ptrcur = 1, this%npp, p_cache_size
  
      ! check if last copy of table and set np
      if( ptrcur + p_cache_size > this%npp ) then
        np = this%npp - ptrcur + 1
      else
        np = p_cache_size
      endif
  
      ! interpolate fields to particles
      call gen_interp_info(this%x, this%dr, noff, np, ptrcur, weight, idx, pcos, psin)
      call interp_field(ef_re, ef_im, max_mode, weight, idx, pcos, psin, np, ep)
      call interp_field(bf_re, bf_im, max_mode, weight, idx, pcos, psin, np, bp)
      call transform_to_cartesian(ep, np, pcos, psin)
      call transform_to_cartesian(bp, np, pcos, psin)
  
      pp = ptrcur
      do i = 1, np
  
        ! normalize E and B fields
        qtmh1 = qtmh / (1.0 - this%qbm * this%psi(pp))
        qtmh2 = qtmh1 * this%gamma(pp)
        ep(:, i) = ep(:, i) * qtmh2
        bp(:, i) = bp(:, i) * qtmh1
  
        ! first half of electric field acceleration
        utmp(:, i) = this%p(:, pp) + ep(:, i)
  
        ! rotation about magnetic field
        this%p(1, pp) = utmp(1, i) + utmp(2, i) * bp(3, i) - utmp(3, i) * bp(2, i)
        this%p(2, pp) = utmp(2, i) + utmp(3, i) * bp(1, i) - utmp(1, i) * bp(3, i)
        this%p(3, pp) = utmp(3, i) + utmp(1, i) * bp(2, i) - utmp(2, i) * bp(1, i)
  
        ostq = 2.0 / ( 1.0 + bp(1, i)**2 + bp(2, i)**2 + bp(3, i)**2 )
        bp(1, i) = bp(1, i) * ostq
        bp(2, i) = bp(2, i) * ostq
        bp(3, i) = bp(3, i) * ostq
  
        utmp(1, i) = utmp(1, i) + this%p(2, pp) * bp(3, i) - this%p(3, pp) * bp(2, i)
        utmp(2, i) = utmp(2, i) + this%p(3, pp) * bp(1, i) - this%p(1, pp) * bp(3, i)
        utmp(3, i) = utmp(3, i) + this%p(1, pp) * bp(2, i) - this%p(2, pp) * bp(1, i)
  
        ! second half of electric field acc.
        this%p(:, pp) = utmp(:, i) + ep(:, i)
  
        ! calculate time-centered gamma for push of position
        this%gamma(pp) = sqrt(1.0 + this%p(1, pp)**2 + this%p(2, pp)**2 + this%p(3, pp)**2)
  
        pp = pp + 1
      enddo
    enddo
    
    call stop_tprof( 'push 2D particles' )
    call write_dbg(cls_name, sname, cls_level, 'ends')
  
end subroutine PUSH_STD

subroutine PUSH_ROB(this, ef, bf, dt_)

    implicit none
  
    class(part2d), intent(inout) :: this
    class(field), intent(in) :: ef, bf
    real, intent(in), optional :: dt_
    ! local data
    character(len=18), save :: sname = 'push_u_robust_part2d'
    type(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im
  
    integer :: i, np, max_mode, noff
    real :: qtmh, qtmh1, qtmh2, gam, dtc, ostq, dt
    real, dimension(p_p_dim, p_cache_size) :: bp, ep, utmp
    real, dimension(:, :), allocatable :: weight
    real, dimension(p_cache_size) :: pcos, psin
    integer, dimension(p_cache_size) :: idx
    integer(kind=LG) :: ptrcur, pp
  
    call write_dbg(cls_name, sname, cls_level, 'starts')
    call start_tprof( 'push 2D particles' )
    allocate(weight(LP:UP, p_cache_size))
  
    dt = this%dt
    if (present(dt_)) dt = dt_
    qtmh = this%qbm * dt * 0.5
    max_mode = ef%get_max_mode()
  
    ef_re => ef%get_rf_re()
    ef_im => ef%get_rf_im()
    bf_re => bf%get_rf_re()
    bf_im => bf%get_rf_im()
    noff = ef_re(0)%get_noff(1)
  
    do ptrcur = 1, this%npp, p_cache_size
  
      ! check if last copy of table and set np
      if( ptrcur + p_cache_size > this%npp ) then
        np = this%npp - ptrcur + 1
      else
        np = p_cache_size
      endif
  
      ! interpolate fields to particles
      call gen_interp_info(this%x, this%dr, noff, np, ptrcur, weight, idx, pcos, psin)
      call interp_field(ef_re, ef_im, max_mode, weight, idx, pcos, psin, np, ep)
      call interp_field(bf_re, bf_im, max_mode, weight, idx, pcos, psin, np, bp)
      call transform_to_cartesian(ep, np, pcos, psin)
      call transform_to_cartesian(bp, np, pcos, psin)
  
      pp = ptrcur
      do i = 1, np
        gam = sqrt( 1.0 + this%p(1, pp)**2 + this%p(2, pp)**2 + this%p(3, pp)**2 )
        qtmh1 = qtmh / (gam - this%p(3, pp))
        qtmh2 = qtmh1 * gam
        ep(:, i) = ep(:, i) * qtmh2
        bp(:, i) = bp(:, i) * qtmh1
  
        ! first half of electric field acceleration
        utmp(:, i) = this%p(:, pp) + ep(:, i)
  
        ! rotation about magnetic field
        this%p(1, pp) = utmp(1, i) + utmp(2, i) * bp(3, i) - utmp(3, i) * bp(2, i)
        this%p(2, pp) = utmp(2, i) + utmp(3, i) * bp(1, i) - utmp(1, i) * bp(3, i)
        this%p(3, pp) = utmp(3, i) + utmp(1, i) * bp(2, i) - utmp(2, i) * bp(1, i)
  
        ostq = 2.0 / ( 1.0 + bp(1, i)**2 + bp(2, i)**2 + bp(3, i)**2 )
        bp(:, i) = bp(:, i) * ostq
  
        utmp(1, i) = utmp(1, i) + this%p(2, pp) * bp(3, i) - this%p(3, pp) * bp(2, i)
        utmp(2, i) = utmp(2, i) + this%p(3, pp) * bp(1, i) - this%p(1, pp) * bp(3, i)
        utmp(3, i) = utmp(3, i) + this%p(1, pp) * bp(2, i) - this%p(2, pp) * bp(1, i)
  
        ! second half of electric field acc.
        this%p(:, pp) = utmp(:, i) + ep(:, i)
  
        ! calculate time-centered gamma for push of position
        this%gamma(pp) = sqrt(1.0 + this%p(1, pp)**2 + this%p(2, pp)**2 + this%p(3, pp)**2)
  
        pp = pp + 1
      enddo
  
    enddo
    
    call stop_tprof( 'push 2D particles' )
    call write_dbg(cls_name, sname, cls_level, 'ends')
  
end subroutine PUSH_ROB
  
subroutine PUSH_ROB_PGC(this, ef, bf, af, dt_)

    implicit none
  
    class(part2d), intent(inout) :: this
    class(field), intent(in) :: ef, bf
    class(field_laser), intent(in) :: af
    real, intent(in), optional :: dt_
    ! local data
    character(len=18), save :: sname = 'push_u_robust_pgc_part2d'
    type(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im
    type(ufield), dimension(:), pointer :: ar_re, ar_im, ai_re, ai_im
    type(ufield), dimension(:), pointer :: ar_grad_re, ar_grad_im, ai_grad_re, ai_grad_im
  
    integer :: i, np, max_mode, noff
    real :: qtmh, qtmh_e, qtmh_b, gam, dtc, ostq, gam_corr, tmp, qbm2_hf, dt
    real, dimension(:, :), allocatable :: weight
    real, dimension(p_cache_size) :: pcos, psin
    integer, dimension(p_cache_size) :: idx
    real, dimension(p_p_dim, p_cache_size) :: bp, ep
    real, dimension(p_cache_size) :: apr, api
    real, dimension(3,p_cache_size) :: apr_grad, api_grad
    real, dimension(p_p_dim) :: utmp
    integer(kind=LG) :: ptrcur, pp
  
    call write_dbg(cls_name, sname, cls_level, 'starts')
    call start_tprof( 'push 2D particles' )
    allocate(weight(LP:UP, p_cache_size))
  
    dt = this%dt
    if (present(dt_)) dt = dt_
  
    qtmh = this%qbm * dt * 0.5
    qbm2_hf = this%qbm**2 * 0.5
    max_mode = ef%get_max_mode()
  
    ef_re => ef%get_rf_re()
    ef_im => ef%get_rf_im()
    bf_re => bf%get_rf_re()
    bf_im => bf%get_rf_im()
  
    ar_re => af%get_cfr_re()
    ar_im => af%get_cfr_im()
    ai_re => af%get_cfi_re()
    ai_im => af%get_cfi_im()
  
    ar_grad_re => af%ar_grad_re
    ar_grad_im => af%ar_grad_im
    ai_grad_re => af%ai_grad_re
    ai_grad_im => af%ai_grad_im
  
    noff = ef_re(0)%get_noff(1)
    do ptrcur = 1, this%npp, p_cache_size
  
      ! check if last copy of table and set np
      if( ptrcur + p_cache_size > this%npp ) then
        np = this%npp - ptrcur + 1
      else
        np = p_cache_size
      endif
  
      ! interpolate fields to particles
      call gen_interp_info(this%x, this%dr, noff, np, ptrcur, weight, idx, pcos, psin)
      call interp_field(ef_re, ef_im, max_mode, weight, idx, pcos, psin, np, ep)
      call interp_field(bf_re, bf_im, max_mode, weight, idx, pcos, psin, np, bp)
      call interp_field(ar_re, ar_im, max_mode, weight, idx, pcos, psin, np, apr)
      call interp_field(ai_re, ai_im, max_mode, weight, idx, pcos, psin, np, api)
      call interp_field(ar_grad_re, ar_grad_im, max_mode, weight, idx, pcos, psin, np, apr_grad)
      call interp_field(ai_grad_re, ai_grad_im, max_mode, weight, idx, pcos, psin, np, api_grad)
      call transform_to_cartesian(ep, np, pcos, psin)
      call transform_to_cartesian(bp, np, pcos, psin)
      call transform_to_cartesian(apr_grad, np, pcos, psin)
      call transform_to_cartesian(api_grad, np, pcos, psin)
  
      ! note that the gamma and momenta calculated in this subroutine are the ones after averaging.
      pp = ptrcur
      do i = 1, np
  
        ! calculate the correction factor for gamma due to ponderomotive force
        gam_corr = qbm2_hf * ( apr(i)**2 + api(i)**2 )
  
        ! calculate the effective electric fields
        ! note that this%gamma is already time-centered
        tmp = 0.5 * this%qbm / this%gamma(pp)
        ep(1,i) = ep(1,i) - tmp * ( apr(i) * apr_grad(1,i) + api(i) * api_grad(1,i) )
        ep(2,i) = ep(2,i) - tmp * ( apr(i) * apr_grad(2,i) + api(i) * api_grad(2,i) )
        ep(3,i) = ep(3,i) + tmp * ( apr(i) * apr_grad(3,i) + api(i) * api_grad(3,i) )
  
        ! half acceleration due to effective electric field
        ! note that this%psi (1 - (q/m)*psi) is already time-centered
        qtmh_b = qtmh / (1.0 - this%qbm * this%psi(pp))
        qtmh_e = qtmh_b * this%gamma(pp)
        ep(:,i) = ep(:,i) * qtmh_e
        utmp(:) = this%p(:,pp) + ep(:,i)
  
        ! rotation about magnetic field
        bp(:,i) = bp(:,i) * qtmh_b
        this%p(1,pp) = utmp(1) + utmp(2) * bp(3,i) - utmp(3) * bp(2,i)
        this%p(2,pp) = utmp(2) + utmp(3) * bp(1,i) - utmp(1) * bp(3,i)
        this%p(3,pp) = utmp(3) + utmp(1) * bp(2,i) - utmp(2) * bp(1,i)
  
        ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
        bp(1,i) = bp(1,i) * ostq
        bp(2,i) = bp(2,i) * ostq
        bp(3,i) = bp(3,i) * ostq
  
        utmp(1) = utmp(1) + this%p(2,pp) * bp(3,i) - this%p(3,pp) * bp(2,i)
        utmp(2) = utmp(2) + this%p(3,pp) * bp(1,i) - this%p(1,pp) * bp(3,i)
        utmp(3) = utmp(3) + this%p(1,pp) * bp(2,i) - this%p(2,pp) * bp(1,i)
  
        ! half acceleration due to effective electric field
        this%p(:,pp) = utmp(:) + ep(:,i)
  
        ! calculate time-centered gamma for push of position
        tmp = apr_grad(3, i) * dt; gam_corr = gam_corr + qbm2_hf * (apr(i) + 0.25 * tmp) * tmp
        tmp = api_grad(3, i) * dt; gam_corr = gam_corr + qbm2_hf * (api(i) + 0.25 * tmp) * tmp
        this%gamma(pp) = sqrt(1.0 + this%p(1, pp)**2 + this%p(2, pp)**2 + this%p(3, pp)**2 + gam_corr)
  
        pp = pp + 1
      enddo
  
    enddo
    
    call stop_tprof( 'push 2D particles' )
    call write_dbg(cls_name, sname, cls_level, 'ends')
  
  end subroutine PUSH_ROB_PGC

  subroutine PUSH_STD_PGC(this, ef, bf, af, dt_)

    implicit none
  
    class(part2d), intent(inout) :: this
    class(field), intent(in) :: ef, bf
    class(field_laser), intent(in) :: af
    real, intent(in), optional :: dt_
    ! local data
    character(len=18), save :: sname = 'push_u_std_pgc_part2d'
    type(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im
    type(ufield), dimension(:), pointer :: ar_re, ar_im, ai_re, ai_im
    type(ufield), dimension(:), pointer :: ar_grad_re, ar_grad_im, ai_grad_re, ai_grad_im
    integer :: i, np, max_mode, noff
    real :: qtmh, qtmh_e, qtmh_b, gam, dtc, ostq, gam_corr, tmp, dt, qbm2_hf
    real, dimension(:, :), allocatable :: weight
    real, dimension(p_cache_size) :: pcos, psin
    integer, dimension(p_cache_size) :: idx
    real, dimension(p_p_dim, p_cache_size) :: bp, ep
    real, dimension(p_cache_size) :: apr, api
    real, dimension(3,p_cache_size) :: apr_grad, api_grad
    real, dimension(p_p_dim) :: utmp
    integer(kind=LG) :: ptrcur, pp
  
    call write_dbg(cls_name, sname, cls_level, 'starts')
    call start_tprof( 'push 2D particles' )
    allocate(weight(LP:UP, p_cache_size))
  
    dt = this%dt
    if (present(dt_)) dt = dt_
  
    qtmh = this%qbm * this%dt * 0.5
    qbm2_hf = this%qbm**2 * 0.5
    max_mode = ef%get_max_mode()
  
    ef_re => ef%get_rf_re()
    ef_im => ef%get_rf_im()
    bf_re => bf%get_rf_re()
    bf_im => bf%get_rf_im()
  
    ar_re => af%get_cfr_re()
    ar_im => af%get_cfr_im()
    ai_re => af%get_cfi_re()
    ai_im => af%get_cfi_im()
  
    ar_grad_re => af%ar_grad_re
    ar_grad_im => af%ar_grad_im
    ai_grad_re => af%ai_grad_re
    ai_grad_im => af%ai_grad_im
  
    noff = ef_re(0)%get_noff(1)
    do ptrcur = 1, this%npp, p_cache_size
  
      ! check if last copy of table and set np
      if( ptrcur + p_cache_size > this%npp ) then
        np = this%npp - ptrcur + 1
      else
        np = p_cache_size
      endif
  
      ! interpolate fields to particles
      call gen_interp_info(this%x, this%dr, noff, np, ptrcur, weight, idx, pcos, psin)
      call interp_field(ef_re, ef_im, max_mode, weight, idx, pcos, psin, np, ep)
      call interp_field(bf_re, bf_im, max_mode, weight, idx, pcos, psin, np, bp)
      call interp_field(ar_re, ar_im, max_mode, weight, idx, pcos, psin, np, apr)
      call interp_field(ai_re, ai_im, max_mode, weight, idx, pcos, psin, np, api)
      call interp_field(ar_grad_re, ar_grad_im, max_mode, weight, idx, pcos, psin, np, apr_grad)
      call interp_field(ai_grad_re, ai_grad_im, max_mode, weight, idx, pcos, psin, np, api_grad)
      call transform_to_cartesian(ep, np, pcos, psin)
      call transform_to_cartesian(bp, np, pcos, psin)
      call transform_to_cartesian(apr_grad, np, pcos, psin)
      call transform_to_cartesian(api_grad, np, pcos, psin)
  
      ! note that the gamma and momenta calculated in this subroutine are the ones after averaging.
      pp = ptrcur
      do i = 1, np
  
        ! calculate the correction factor for gamma due to ponderomotive force
        gam_corr = qbm2_hf * (apr(i)**2 + api(i)**2)
  
        ! calculate the effective electric fields
        ! note that this%gamma is already time-centered
        tmp = 0.5 * this%qbm / this%gamma(pp)
        ep(1, i) = ep(1, i) - tmp * ( apr(i) * apr_grad(1, i) + api(i) * api_grad(1, i) )
        ep(2, i) = ep(2, i) - tmp * ( apr(i) * apr_grad(2, i) + api(i) * api_grad(2, i) )
        ep(3, i) = ep(3, i) + tmp * ( apr(i) * apr_grad(3, i) + api(i) * api_grad(3, i) )
  
        ! normalize the fields
        qtmh_b = qtmh / (1.0 - this%qbm * this%psi(pp))
        qtmh_e = qtmh_b * this%gamma(pp)
        ep(:, i) = ep(:, i) * qtmh_e
        bp(:, i) = bp(:, i) * qtmh_b
  
        ! half acceleration due to effective electric field
        utmp(:) = this%p(:, pp) + ep(:, i)
  
        ! rotation about magnetic field
        this%p(1, pp) = utmp(1) + utmp(2) * bp(3, i) - utmp(3) * bp(2, i)
        this%p(2, pp) = utmp(2) + utmp(3) * bp(1, i) - utmp(1) * bp(3, i)
        this%p(3, pp) = utmp(3) + utmp(1) * bp(2, i) - utmp(2) * bp(1, i)
  
        ostq = 2.0 / (1.0 + bp(1, i)**2 + bp(2, i)**2 + bp(3, i)**2)
        bp(1, i) = bp(1, i) * ostq
        bp(2, i) = bp(2, i) * ostq
        bp(3, i) = bp(3, i) * ostq
  
        utmp(1) = utmp(1) + this%p(2, pp) * bp(3, i) - this%p(3, pp) * bp(2, i)
        utmp(2) = utmp(2) + this%p(3, pp) * bp(1, i) - this%p(1, pp) * bp(3, i)
        utmp(3) = utmp(3) + this%p(1, pp) * bp(2, i) - this%p(2, pp) * bp(1, i)
  
        ! half acceleration due to effective electric field
        this%p(:, pp) = utmp(:) + ep(:, i)
  
        ! calculate time-centered gamma for push of position
        tmp = apr_grad(3, i) * dt; gam_corr = gam_corr + qbm2_hf * (apr(i) + 0.25 * tmp) * tmp
        tmp = api_grad(3, i) * dt; gam_corr = gam_corr + qbm2_hf * (api(i) + 0.25 * tmp) * tmp
        this%gamma(pp) = sqrt(1.0 + this%p(1, pp)**2 + this%p(2, pp)**2 + this%p(3, pp)**2 + gam_corr)
  
        pp = pp + 1
      enddo
  
    enddo
    
    call stop_tprof( 'push 2D particles' )
    call write_dbg(cls_name, sname, cls_level, 'ends')
  
end subroutine PUSH_STD_PGC
  
#undef SPLINE
#undef CHIDEPOSIT
#undef QDEPOSIT
! #undef AX_COOR
#undef LINEAR
#undef QUADRATIC
! #undef CORRNUM
#undef LP
#undef UP
#undef AMJ_STD
#undef AMJ_ROB
#undef AMJ_STD_PGC
#undef AMJ_ROB_PGC
#undef PUSH_STD
#undef PUSH_ROB
#undef PUSH_ROB_PGC
#undef PUSH_STD_PGC

#endif