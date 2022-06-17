module part2d_class

use param
use sysutil_module
use parallel_module
use sort_module
use options_class
use field_class
use ufield_class
use fdist2d_class
use hdf5io_class
use mpi
use interpolation
use debug_tool

implicit none

private

public :: part2d

type part2d

   ! qbm = particle charge/mass ratio
   ! dt = time interval between successive calculations
   ! dr = radial cell size
   real :: qbm, dt, dr

   ! maximum effective time step (only used for sub-cycling)
   real :: dt_eff_max

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
   real, dimension(:), allocatable :: w
   real, dimension(:), allocatable :: w0
   ! particle upper boundaries
   real :: edge
   ! particle buffer
   real, dimension(:), allocatable :: pbuf
   ! clamped value of gamma/ (1 + psi)
   real :: fac_clamp

   ! temporary arrays used for buffer reallocation
   real, private, dimension(:), allocatable :: tmp1
   real, private, dimension(:,:), allocatable :: tmp2

   contains

   procedure :: new                      => init_part2d
   procedure :: renew                    => renew_part2d
   procedure :: del                      => end_part2d
   procedure :: qdeposit                 => qdeposit_part2d
   procedure :: amjdeposit_robust        => amjdeposit_robust_part2d
   procedure :: amjdeposit_clamp         => amjdeposit_clamp_part2d
   procedure :: amjdeposit_robust_subcyc => amjdeposit_robust_subcyc_part2d
   procedure :: ionize                   => ionize_part2d
   procedure :: add_particles            => add_particles_part2d
   procedure :: push_robust              => push_robust_part2d
   procedure :: push_clamp               => push_clamp_part2d
   procedure :: push_robust_subcyc       => push_robust_subcyc_part2d
   procedure :: update_bound             => update_bound_part2d
   ! procedure :: extract_psi            => extract_psi_part2d
   procedure :: pipesend                 => pipesend_part2d
   procedure :: piperecv                 => piperecv_part2d
   procedure :: wr                       => writehdf5_part2d
   procedure :: realloc                  => realloc_part2d
   procedure :: sort                     => sort_part2d

end type part2d

save

character(len=20), parameter :: cls_name = "part2d"
integer, parameter :: cls_level = 2

real, dimension(:), allocatable :: recv_buf
integer :: recv_buf_size = 0

integer, dimension(:), allocatable :: sort_idx

integer, parameter :: p_max_subcyc = 1024
real, parameter :: p_buf_incr = 1.5

contains

subroutine init_part2d( this, opts, pf, qbm, dt, s, if_empty, ionization )

   implicit none

   class(part2d), intent(inout) :: this
   type(options), intent(in) :: opts
   class(fdist2d), intent(inout) :: pf
   real, intent(in) :: qbm, dt, s
   logical, intent(in), optional :: if_empty
   logical, intent(in), optional :: ionization

   ! local data
   character(len=18), save :: sname = 'init_part2d'
   integer :: npmax
   logical :: empty = .false.
   logical :: ionize = .false.

   call write_dbg( cls_name, sname, cls_level, 'starts' )

   this%qbm = qbm
   this%dt  = dt
   this%dt_eff_max = pf%dt_eff_max
   this%fac_clamp = pf%fac_clamp
   this%part_dim = 2 + p_p_dim + 5

   npmax      = pf%npmax
   this%npmax = npmax
   this%npp   = 0

   this%dr   = opts%get_dr()
   this%edge = opts%get_nd(1) * this%dr
   
   allocate( this%x( 2, npmax ) )
   allocate( this%p( p_p_dim, npmax ) )
   allocate( this%gamma( npmax ), this%q( npmax ), this%psi( npmax ), this%w( npmax ), this%w0( npmax ) )
   allocate( this%pbuf( this%part_dim * npmax ) )

   recv_buf_size = max( recv_buf_size, npmax )

   if ( present( if_empty ) ) empty = if_empty
   if ( present( ionization ) ) ionize = ionization

   ! initialize particle coordinates according to specified profile
   if ( .not. empty ) call pf%inject( this%x, this%p, this%gamma, this%psi, this%q, this%w, this%w0, this%npp, s, ionize )


   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_part2d

subroutine end_part2d(this)

   implicit none

   class(part2d), intent(inout) :: this
   character(len=18), save :: sname = 'end_part2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   deallocate( this%x, this%p, this%gamma, this%q, this%psi, this%w )

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
  this%tmp1( 1:this%npp ) = this%w( 1:this%npp )
  call move_alloc( this%tmp1, this%w )

  allocate( this%tmp1( npmax ) )
  this%tmp1 = 0.0
  this%tmp1( 1:this%npp ) = this%w0( 1:this%npp )
  call move_alloc( this%tmp1, this%w0 )

  allocate( this%tmp1( npmax ) )
  this%tmp1 = 0.0
  this%tmp1( 1:this%npp ) = this%q( 1:this%npp )
  call move_alloc( this%tmp1, this%q )

  deallocate( this%pbuf )
  allocate( this%pbuf( this%part_dim * npmax ) )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine realloc_part2d

subroutine renew_part2d( this, pf, s, if_empty, ionization )

   implicit none

   class(part2d), intent(inout) :: this
   class(fdist2d), intent(inout) :: pf
   real, intent(in) :: s
   logical, intent(in), optional :: if_empty
   logical, intent(in), optional :: ionization

   ! local data
   logical :: empty = .false.
   logical :: ionize = .false.
   character(len=18), save :: sname = 'renew_part2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%npp = 0

   if ( present( if_empty ) ) empty = if_empty
   if ( present( ionization ) ) ionize = ionization
   if ( .not. empty ) call pf%inject( this%x, this%p, this%gamma, this%psi, this%q, this%w, this%w0, this%npp, s, ionize )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine renew_part2d

subroutine qdeposit_part2d( this, q )
! deposit the charge density (rho - Jz)

  implicit none

  class(part2d), intent(in) :: this
  class(field), intent(in) :: q
  ! local data
  type(ufield), dimension(:), pointer :: q_re => null(), q_im => null()
  complex(kind=DB), dimension(p_cache_size) :: phase0
  complex(kind=DB) :: phase
  real, dimension(:,:), pointer :: q0 => null(), qr => null(), qi => null()
  real, dimension(0:1) :: wt ! interpolation weight
  real, dimension(p_cache_size) :: pos ! normalized position
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
       pos(i) = sqrt( this%x(1, pp)**2 + this%x(2, pp)**2 ) * idr
       phase0(i) = cmplx( this%x(1, pp), -this%x(2, pp) ) / pos(i) * idr
       pp = pp + 1
    enddo

    pp = ptrcur
    do i = 1, np
      nn = floor( pos(i) )
      ! in-cell position
      pos(i) = pos(i) - real(nn)
      nn = nn - noff + 1

      call spline_linear( pos(i), wt )

      phase = cmplx( 1.0, 0.0 ) * this%q(pp)
      ! deposit m=0 mode
      do j = 0, 1
        q0( 1, nn+j ) = q0( 1, nn+j ) + wt(j) * real(phase)
      enddo

      ! deposit m>0 mode
      do mode = 1, max_mode
        qr => q_re(mode)%get_f1()
        qi => q_im(mode)%get_f1()
        phase = phase * phase0(i)

        do j = 0, 1
          qr( 1, nn+j ) = qr( 1, nn+j ) + wt(j) * real(phase)
          qi( 1, nn+j ) = qi( 1, nn+j ) + wt(j) * aimag(phase)
        enddo
      enddo
       pp = pp + 1
    enddo

  enddo

  if (noff == 0) then

    q0(1,0) = 0.0 ! guard cell is useless on axis
    q0(1,1) = 8.0 * q0(1,1)
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

end subroutine qdeposit_part2d

subroutine amjdeposit_robust_part2d( this, ef, bf, cu, amu, dcu )
! deposit the current, acceleration and momentum flux

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: cu, amu, dcu
  class(field), intent(in) :: ef, bf
  ! local data
  character(len=32), save :: sname = 'amjdeposit_robust_part2d'
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
  real, dimension(p_p_dim, p_cache_size) :: bp, ep, wp, u0, u, utmp
  real, dimension(0:1, p_cache_size) :: wt
  real, dimension(p_cache_size) :: cc, ss
  real, dimension(p_p_dim) :: du, u2
  real :: qtmh, qtmh1, qtmh2, idt, gam, ostq, ipsi, dpsi, w, ir
  complex(kind=DB) :: phase, phase0

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'deposit 2D particles' )

  ef_re  => ef%get_rf_re();  ef_im  => ef%get_rf_im()
  bf_re  => bf%get_rf_re();  bf_im  => bf%get_rf_im()
  cu_re  => cu%get_rf_re();  cu_im  => cu%get_rf_im()
  dcu_re => dcu%get_rf_re(); dcu_im => dcu%get_rf_im()
  amu_re => amu%get_rf_re(); amu_im => amu%get_rf_im()

  idt = 1.0 / this%dt
  qtmh = 0.5 * this%qbm * this%dt
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
    call interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, this%x, this%dr, &
      bp, ep, np, ptrcur, p_cylindrical, weight = wt, ix = ix, pcos = cc, psin = ss )

    ! calculate wake field
    do i = 1, np
      wp(1,i) = ep(1,i) - bp(2,i)
      wp(2,i) = ep(2,i) + bp(1,i)
      wp(3,i) = ep(3,i)
    enddo

    ! transform momentum from Cartesian to cylindrical coordinates
    pp = ptrcur
    do i = 1, np
      u0(1,i) = this%p(1,pp) * cc(i) + this%p(2,pp) * ss(i)
      u0(2,i) = this%p(2,pp) * cc(i) - this%p(1,pp) * ss(i)
      u0(3,i) = this%p(3,pp)
      pp = pp + 1
    enddo

    ! half electric acceleration
    do i = 1, np
      gam = sqrt( 1.0 + u0(1,i)**2 + u0(2,i)**2 + u0(3,i)**2 )
      qtmh1 = qtmh * gam / ( gam - u0(3,i) )
      ep(:,i) = ep(:,i) * qtmh1
      utmp(:,i) = u0(:,i) + ep(:,i)
    enddo

    ! scale magnetic field
    do i = 1, np
      gam = sqrt( 1.0 + utmp(1,i)**2 + utmp(2,i)**2 + utmp(3,i)**2 )
      qtmh2 = qtmh / ( gam - utmp(3,i) )
      bp(:,i) = bp(:,i) * qtmh2
    enddo

    ! magnetic rotation
    do i = 1, np
      u(1,i) = utmp(1,i) + utmp(2,i) * bp(3,i) - utmp(3,i) * bp(2,i)
      u(2,i) = utmp(2,i) + utmp(3,i) * bp(1,i) - utmp(1,i) * bp(3,i)
      u(3,i) = utmp(3,i) + utmp(1,i) * bp(2,i) - utmp(2,i) * bp(1,i)

      ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
      bp(:,i) = bp(:,i) * ostq

      utmp(1,i) = utmp(1,i) + u(2,i) * bp(3,i) - u(3,i) * bp(2,i)
      utmp(2,i) = utmp(2,i) + u(3,i) * bp(1,i) - u(1,i) * bp(3,i)
      utmp(3,i) = utmp(3,i) + u(1,i) * bp(2,i) - u(2,i) * bp(1,i)
    enddo

    ! half electric acceleration
    do i = 1, np
      u(:,i) = utmp(:,i) + ep(:,i)
    enddo

    ! calculate and store time-centered values
    ! deposit momentum flux, acceleration density, and current density
    pp = ptrcur
    do i = 1, np

      du(1) = idt * ( u(1,i) - u0(1,i) )
      du(2) = idt * ( u(2,i) - u0(2,i) )
      
      u(:,i)  = 0.5 * ( u(:,i) + u0(:,i) )
      this%gamma(pp) = sqrt( 1.0 + u(1,i)**2 + u(2,i)**2 + u(3,i)**2 )
      this%psi(pp)   = this%gamma(pp) - u(3,i)

      ipsi = 1.0 / this%psi(pp)
      dpsi = this%qbm * ( wp(3,i) - ( wp(1,i) * u(1,i) + wp(2,i) * u(2,i) ) * ipsi )

      du(1) = du(1) + u(1,i) * dpsi * ipsi
      du(2) = du(2) + u(2,i) * dpsi * ipsi

      u2(1) = u(1,i) * u(1,i) * ipsi
      u2(2) = u(1,i) * u(2,i) * ipsi
      u2(3) = u(2,i) * u(2,i) * ipsi

      phase0 = cmplx( cc(i), -ss(i) )
      phase  = cmplx( 1.0, 0.0 ) * this%q(pp) * ipsi

      ! deposit m = 0 mode
      do j = 0, 1
        w = wt(j,i) * real(phase)
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

        do j = 0, 1
          w = wt(j,i) * real(phase)
          cur( 1:3, ix(i)+j )  = cur( 1:3, ix(i)+j )  + w * u(1:3,i)
          dcur( 1:2, ix(i)+j ) = dcur( 1:2, ix(i)+j ) + w * du(1:2)
          amur( 1:3, ix(i)+j ) = amur( 1:3, ix(i)+j ) + w * u2(1:3)

          w = wt(j,i) * aimag(phase)
          cui( 1:3, ix(i)+j )  = cui( 1:3, ix(i)+j )  + w * u(1:3,i)
          dcui( 1:2, ix(i)+j ) = dcui( 1:2, ix(i)+j ) + w * du(1:2)
          amui( 1:3, ix(i)+j ) = amui( 1:3, ix(i)+j ) + w * u2(1:3)
        enddo

      enddo

      pp = pp + 1
    enddo

  enddo

  if ( noff == 0 ) then

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
        cur(1:2,1)  = 8.0 * cur(1:2,1); cur(3,1) = 0.0
        dcur(1:2,1) = 8.0 * dcur(1:2,1)
        amur(1:3,1) = 0.0
        cui(1:2,1)  = 8.0 * cui(1:2,1); cui(3,1) = 0.0
        dcui(1:2,1) = 8.0 * dcui(1:2,1)
        amui(1:3,1) = 0.0
      elseif ( mode == 2 ) then
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

end subroutine amjdeposit_robust_part2d

subroutine amjdeposit_robust_subcyc_part2d( this, ef, bf, cu, amu, dcu )
! deposit the current, acceleration and momentum flux

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: cu, amu, dcu
  class(field), intent(in) :: ef, bf
  ! local data
  character(len=18), save :: sname = 'amjdeposit_robust_subcyc_part2d'
  type(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
  type(ufield), dimension(:), pointer :: bf_re => null(), bf_im => null()
  type(ufield), dimension(:), pointer :: cu_re => null(), cu_im => null()
  type(ufield), dimension(:), pointer :: dcu_re => null(), dcu_im => null()
  type(ufield), dimension(:), pointer :: amu_re => null(), amu_im => null()

  real, dimension(:,:), pointer :: cu0 => null(), dcu0 => null(), amu0 => null()
  real, dimension(:,:), pointer :: cur => null(), dcur => null(), amur => null()
  real, dimension(:,:), pointer :: cui => null(), dcui => null(), amui => null()

  integer(kind=LG) :: ptrcur, pp
  integer :: i, j, noff, nrp, np, np_subcyc, mode, max_mode, n_subcyc, exit_cnt, n_subcyc_max, np_fail
  integer, dimension(p_cache_size) :: ix, ix_subcyc, ndt_rem
  real, dimension(p_p_dim, p_cache_size) :: bp, ep, wp, p_old, p_new, p_subcyc
  real, dimension(2, p_cache_size) :: x_subcyc
  real, dimension(0:1, p_cache_size) :: wt
  real, dimension(p_cache_size) :: cc, ss, gam_subcyc, dt_subcyc
  real, dimension(p_p_dim) :: du, u2, utmp
  real :: qtmh1, qtmh2, idt, ostq, ipsi, dpsi, w, ir, dtc, pcos, psin
  complex(kind=DB) :: phase, phase0

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'deposit 2D particles' )

  ef_re  => ef%get_rf_re();  ef_im  => ef%get_rf_im()
  bf_re  => bf%get_rf_re();  bf_im  => bf%get_rf_im()
  cu_re  => cu%get_rf_re();  cu_im  => cu%get_rf_im()
  dcu_re => dcu%get_rf_re(); dcu_im => dcu%get_rf_im()
  amu_re => amu%get_rf_re(); amu_im => amu%get_rf_im()

  idt = 1.0 / this%dt
  max_mode = ef%get_max_mode()

  noff = cu_re(0)%get_noff(1)
  nrp  = cu_re(0)%get_ndp(1)

  cu0  => cu_re(0)%get_f1()
  dcu0 => dcu_re(0)%get_f1()
  amu0 => amu_re(0)%get_f1()

  n_subcyc_max = 0
  np_fail = 0

  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! get the time-centered fields by interpolation
    call interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, this%x, this%dr, &
      bp, ep, np, ptrcur, p_cylindrical, weight = wt, ix = ix, pcos = cc, psin = ss )

    ! calculate wake field and transform momentum from Cartesian to cylindrical coordinates
    pp = ptrcur
    do i = 1, np
      wp(1,i) = ep(1,i) - bp(2,i)
      wp(2,i) = ep(2,i) + bp(1,i)
      wp(3,i) = ep(3,i)
      pp = pp + 1
    enddo

    ! initialize sub-cycling
    pp = ptrcur
    do i = 1, np

      ! momentum for sub-cycling is in Cartesian coordinates
      p_subcyc(:,i) = this%p(:,pp)

      ! synchronize the particle position with the momentum for sub-cycling
      gam_subcyc(i) = sqrt( 1.0 + p_subcyc(1,i)**2 + p_subcyc(2,i)**2 + p_subcyc(3,i)**2 )
      dtc = 0.5 * this%dt / ( gam_subcyc(i) - p_subcyc(3,i) )
      x_subcyc(1,i) = this%x(1,pp) - dtc * p_subcyc(1,i)
      x_subcyc(2,i) = this%x(2,pp) - dtc * p_subcyc(2,i)

      ! get the old momentum in cylindrical coordinates
      ir = 1.0 / sqrt( x_subcyc(1,i)**2 + x_subcyc(2,i)**2 )
      pcos = x_subcyc(1,i) * ir
      psin = x_subcyc(2,i) * ir
      p_old(1,i) = p_subcyc(1,i) * pcos + p_subcyc(2,i) * psin
      p_old(2,i) = p_subcyc(2,i) * pcos - p_subcyc(1,i) * psin
      p_old(3,i) = p_subcyc(3,i)
      
      ! particle index for sub-cycling
      ix_subcyc(i) = pp      
      
      ! calculate the sub-cycling time step
      ndt_rem(i) = ceiling( this%dt / min( this%dt, this%dt_eff_max * (1.0 - p_subcyc(3,i) / gam_subcyc(i)) ) )
      dt_subcyc(i) = this%dt / ndt_rem(i)

      ! advance position by a half sub-cycling time step
      dtc = 0.5 * dt_subcyc(i) / ( gam_subcyc(i) - p_subcyc(3,i) )
      x_subcyc(1,i) = x_subcyc(1,i) + p_subcyc(1,i) * dtc
      x_subcyc(2,i) = x_subcyc(2,i) + p_subcyc(2,i) * dtc

      pp = pp + 1
    enddo

    n_subcyc = 0
    np_subcyc = np
    ! begin sub-cycling
    do while ( np_subcyc > 0 )

      ! interpolate fields to particles
      call interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, x_subcyc, this%dr, &
        bp, ep, np_subcyc, int(1, kind=LG), p_cartesian )

      do i = 1, np_subcyc

        qtmh1 = 0.5 * dt_subcyc(i) * this%qbm / ( gam_subcyc(i) - p_subcyc(3,i) )
        qtmh2 = qtmh1 * gam_subcyc(i)

        ! scale the fields
        ep(:,i) = ep(:,i) * qtmh2
        bp(:,i) = bp(:,i) * qtmh1

        ! first half of electric field acceleration
        utmp = p_subcyc(:,i) + ep(:,i)

        ! rotation about magnetic field
        p_subcyc(1,i) = utmp(1) + utmp(2) * bp(3,i) - utmp(3) * bp(2,i)
        p_subcyc(2,i) = utmp(2) + utmp(3) * bp(1,i) - utmp(1) * bp(3,i)
        p_subcyc(3,i) = utmp(3) + utmp(1) * bp(2,i) - utmp(2) * bp(1,i)

        ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
        bp(:,i) = bp(:,i) * ostq

        utmp(1) = utmp(1) + p_subcyc(2,i) * bp(3,i) - p_subcyc(3,i) * bp(2,i)
        utmp(2) = utmp(2) + p_subcyc(3,i) * bp(1,i) - p_subcyc(1,i) * bp(3,i)
        utmp(3) = utmp(3) + p_subcyc(1,i) * bp(2,i) - p_subcyc(2,i) * bp(1,i)

        ! second half of electric field acceleration
        p_subcyc(:,i) = utmp + ep(:,i)

        ! update the remaining time step numbers
        ndt_rem(i) = ndt_rem(i) - 1

        ! advance position
        gam_subcyc(i) = sqrt( 1.0 + p_subcyc(1,i)**2 + p_subcyc(2,i)**2 + p_subcyc(3,i)**2 )
        dtc = dt_subcyc(i) / ( gam_subcyc(i) - p_subcyc(3,i) )
        x_subcyc(1,i) = x_subcyc(1,i) + p_subcyc(1,i) * dtc
        x_subcyc(2,i) = x_subcyc(2,i) + p_subcyc(2,i) * dtc

      enddo

      ! store the particles that finished sub-cycling and rearrange sub-cycling array
      exit_cnt = 0; i = 1
      do while ( i <= np_subcyc - exit_cnt )

        ! the particles that have used up the global 2D time step or went out of the boundary
        ! will exit the sub-cycling.
        ir = sqrt( x_subcyc(1,i)**2 + x_subcyc(2,i)**2 )
        if ( ndt_rem(i) == 0 .or. ir >= this%edge ) then

          ! store the new momentum in cylindrical coordinates
          pp = ix_subcyc(i) - ptrcur + 1
          pcos = x_subcyc(1,i) / ir
          psin = x_subcyc(2,i) / ir
          p_new(1,pp) = p_subcyc(1,i) * pcos + p_subcyc(2,i) * psin
          p_new(2,pp) = p_subcyc(2,i) * pcos - p_subcyc(1,i) * psin
          p_new(3,pp) = p_subcyc(3,i)

          ! move the last element of sub-cycling array to current position
          x_subcyc(:,i) = x_subcyc( :, np_subcyc - exit_cnt )
          p_subcyc(:,i) = p_subcyc( :, np_subcyc - exit_cnt )
          ix_subcyc(i)  = ix_subcyc( np_subcyc - exit_cnt )
          dt_subcyc(i)  = dt_subcyc( np_subcyc - exit_cnt )
          gam_subcyc(i) = gam_subcyc( np_subcyc - exit_cnt )
          ndt_rem(i)    = ndt_rem( np_subcyc - exit_cnt )

          exit_cnt = exit_cnt + 1

        else
          i = i + 1
        endif

      enddo
      np_subcyc = np_subcyc - exit_cnt
      n_subcyc = n_subcyc + 1

      ! when reaching the max number of sub-cycling, store the unfinished particles
      ! and pop out a warning.
      if ( n_subcyc == p_max_subcyc ) then
        np_fail = np_fail + np_subcyc
        do i = 1, np_subcyc
          pp = ix_subcyc(i) - ptrcur + 1
          ir = 1.0 / sqrt( x_subcyc(1,i)**2 + x_subcyc(2,i)**2 )
          pcos = x_subcyc(1,i) * ir
          psin = x_subcyc(2,i) * ir
          p_new(1,pp) = p_subcyc(1,i) * pcos + p_subcyc(2,i) * psin
          p_new(2,pp) = p_subcyc(2,i) * pcos - p_subcyc(1,i) * psin
          p_new(3,pp) = p_subcyc(3,i)
        enddo
        exit
      endif

    enddo ! sub-cycling

    n_subcyc_max = max(n_subcyc_max, n_subcyc)

    ! calculate and store time-centered values
    ! deposit momentum flux, acceleration density, and current density
    pp = ptrcur
    do i = 1, np

      du(1) = idt * ( p_new(1,i) - p_old(1,i) )
      du(2) = idt * ( p_new(2,i) - p_old(2,i) )
      
      ! store time-centered values
      utmp  = 0.5 * ( p_new(:,i) + p_old(:,i) )
      this%gamma(pp) = sqrt( 1.0 + utmp(1)**2 + utmp(2)**2 + utmp(3)**2 )
      this%psi(pp)   = this%gamma(pp) - utmp(3)

      ipsi = 1.0 / this%psi(pp)
      dpsi = this%qbm * ( wp(3,i) - ( wp(1,i) * utmp(1) + wp(2,i) * utmp(2) ) * ipsi )

      du(1) = du(1) + utmp(1) * dpsi * ipsi
      du(2) = du(2) + utmp(2) * dpsi * ipsi

      u2(1) = utmp(1) * utmp(1) * ipsi
      u2(2) = utmp(1) * utmp(2) * ipsi
      u2(3) = utmp(2) * utmp(2) * ipsi

      phase0 = cmplx( cc(i), -ss(i) )
      phase  = cmplx( 1.0, 0.0 ) * this%q(pp) * ipsi

      ! deposit m = 0 mode
      do j = 0, 1
        w = wt(j,i) * real(phase)
        cu0( 1:3, ix(i)+j )  = cu0( 1:3, ix(i)+j )  + w * utmp(1:3)
        dcu0( 1:2, ix(i)+j ) = dcu0( 1:2, ix(i)+j ) + w * du(1:2)
        amu0( 1:3, ix(i)+j ) = amu0( 1:3, ix(i)+j ) + w * u2(1:3)
      enddo

      ! deposit m > 0 mode
      do mode = 1, max_mode

        cur  => cu_re(mode)%get_f1();  cui  => cu_im(mode)%get_f1()
        dcur => dcu_re(mode)%get_f1(); dcui => dcu_im(mode)%get_f1()
        amur => amu_re(mode)%get_f1(); amui => amu_im(mode)%get_f1()

        phase = phase * phase0

        do j = 0, 1
          w = wt(j,i) * real(phase)
          cur( 1:3, ix(i)+j )  = cur( 1:3, ix(i)+j )  + w * utmp(1:3)
          dcur( 1:2, ix(i)+j ) = dcur( 1:2, ix(i)+j ) + w * du(1:2)
          amur( 1:3, ix(i)+j ) = amur( 1:3, ix(i)+j ) + w * u2(1:3)

          w = wt(j,i) * aimag(phase)
          cui( 1:3, ix(i)+j )  = cui( 1:3, ix(i)+j )  + w * utmp(1:3)
          dcui( 1:2, ix(i)+j ) = dcui( 1:2, ix(i)+j ) + w * du(1:2)
          amui( 1:3, ix(i)+j ) = amui( 1:3, ix(i)+j ) + w * u2(1:3)
        enddo

      enddo

      pp = pp + 1
    enddo

  enddo ! chunk loop

  ! call write_stdout( "[amj] max sub-cycling times = " // num2str(n_subcyc_max) )
  if ( np_fail > 0 ) then
    call write_stdout( '[amj] Max number of sub-cycling reached. ' // num2str(np_fail) // &
            ' particles have not yet finished sub-cycling.' )
  endif

  if ( noff == 0 ) then

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
        cur(1:2,1)  = 8.0 * cur(1:2,1); cur(3,1) = 0.0
        dcur(1:2,1) = 8.0 * dcur(1:2,1)
        amur(1:3,1) = 0.0
        cui(1:2,1)  = 8.0 * cui(1:2,1); cui(3,1) = 0.0
        dcui(1:2,1) = 8.0 * dcui(1:2,1)
        amui(1:3,1) = 0.0
      elseif ( mode == 2 ) then
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

end subroutine amjdeposit_robust_subcyc_part2d

subroutine amjdeposit_clamp_part2d( this, ef, bf, cu, amu, dcu )
! deposit the current, acceleration and momentum flux

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: cu, amu, dcu
  class(field), intent(in) :: ef, bf
  ! local data
  character(len=18), save :: sname = 'amjdeposit_clamp_part2d'
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
  real, dimension(p_p_dim, p_cache_size) :: bp, ep, wp, u0, u, utmp
  real, dimension(0:1, p_cache_size) :: wt
  real, dimension(p_cache_size) :: cc, ss
  real, dimension(p_p_dim) :: du, u2
  real :: qtmh, qtmh1, qtmh2, idt, gam, ostq, ipsi, dpsi, w, ir, psi_plus1
  complex(kind=DB) :: phase, phase0

  integer :: stat

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'deposit 2D particles' )

  ef_re  => ef%get_rf_re();  ef_im  => ef%get_rf_im()
  bf_re  => bf%get_rf_re();  bf_im  => bf%get_rf_im()
  cu_re  => cu%get_rf_re();  cu_im  => cu%get_rf_im()
  dcu_re => dcu%get_rf_re(); dcu_im => dcu%get_rf_im()
  amu_re => amu%get_rf_re(); amu_im => amu%get_rf_im()

  idt = 1.0 / this%dt
  qtmh = 0.5 * this%qbm * this%dt
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
    call interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, this%x, this%dr, &
      bp, ep, np, ptrcur, p_cylindrical, weight = wt, ix = ix, pcos = cc, psin = ss )

    ! calculate wake field
    do i = 1, np
      wp(1,i) = ep(1,i) - bp(2,i)
      wp(2,i) = ep(2,i) + bp(1,i)
      wp(3,i) = ep(3,i)
    enddo

    ! transform momentum from Cartesian to cylindrical coordinates
    pp = ptrcur
    do i = 1, np
      u0(1,i) = this%p(1,pp) * cc(i) + this%p(2,pp) * ss(i)
      u0(2,i) = this%p(2,pp) * cc(i) - this%p(1,pp) * ss(i)
      u0(3,i) = this%p(3,pp)
      pp = pp + 1
    enddo

    ! half electric acceleration
    do i = 1, np
      gam = sqrt( 1.0 + u0(1,i)**2 + u0(2,i)**2 + u0(3,i)**2 )

      ! clamp the value of gamma/(1 + psi)
      psi_plus1 = gam - u0(3,i)
      call clamp_momentum( this%fac_clamp, u0(:,i), gam, psi_plus1 )

      qtmh1 = qtmh * gam / psi_plus1
      ep(:,i) = ep(:,i) * qtmh1
      utmp(:,i) = u0(:,i) + ep(:,i)
    enddo

    ! scale magnetic field
    do i = 1, np
      gam = sqrt( 1.0 + utmp(1,i)**2 + utmp(2,i)**2 + utmp(3,i)**2 )

      psi_plus1 = gam - utmp(3,i)
      qtmh2 = qtmh / psi_plus1
      bp(:,i) = bp(:,i) * qtmh2
    enddo

    ! magnetic rotation
    do i = 1, np
      u(1,i) = utmp(1,i) + utmp(2,i) * bp(3,i) - utmp(3,i) * bp(2,i)
      u(2,i) = utmp(2,i) + utmp(3,i) * bp(1,i) - utmp(1,i) * bp(3,i)
      u(3,i) = utmp(3,i) + utmp(1,i) * bp(2,i) - utmp(2,i) * bp(1,i)

      ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
      bp(:,i) = bp(:,i) * ostq

      utmp(1,i) = utmp(1,i) + u(2,i) * bp(3,i) - u(3,i) * bp(2,i)
      utmp(2,i) = utmp(2,i) + u(3,i) * bp(1,i) - u(1,i) * bp(3,i)
      utmp(3,i) = utmp(3,i) + u(1,i) * bp(2,i) - u(2,i) * bp(1,i)
    enddo

    ! half electric acceleration
    do i = 1, np
      u(:,i) = utmp(:,i) + ep(:,i)
    enddo

    ! calculate and store time-centered values
    ! deposit momentum flux, acceleration density, and current density
    pp = ptrcur
    do i = 1, np

      du(1) = idt * ( u(1,i) - u0(1,i) )
      du(2) = idt * ( u(2,i) - u0(2,i) )
      
      u(:,i)  = 0.5 * ( u(:,i) + u0(:,i) )
      this%gamma(pp) = sqrt( 1.0 + u(1,i)**2 + u(2,i)**2 + u(3,i)**2 )

      this%psi(pp)   = this%gamma(pp) - u(3,i)

      ipsi = 1.0 / this%psi(pp)
      dpsi = this%qbm * ( wp(3,i) - ( wp(1,i) * u(1,i) + wp(2,i) * u(2,i) ) * ipsi )

      du(1) = du(1) + u(1,i) * dpsi * ipsi
      du(2) = du(2) + u(2,i) * dpsi * ipsi

      u2(1) = u(1,i) * u(1,i) * ipsi
      u2(2) = u(1,i) * u(2,i) * ipsi
      u2(3) = u(2,i) * u(2,i) * ipsi

      phase0 = cmplx( cc(i), -ss(i) )
      phase  = cmplx( 1.0, 0.0 ) * this%q(pp) * ipsi

      ! deposit m = 0 mode
      do j = 0, 1
        w = wt(j,i) * real(phase)
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

        do j = 0, 1
          w = wt(j,i) * real(phase)
          cur( 1:3, ix(i)+j )  = cur( 1:3, ix(i)+j )  + w * u(1:3,i)
          dcur( 1:2, ix(i)+j ) = dcur( 1:2, ix(i)+j ) + w * du(1:2)
          amur( 1:3, ix(i)+j ) = amur( 1:3, ix(i)+j ) + w * u2(1:3)

          w = wt(j,i) * aimag(phase)
          cui( 1:3, ix(i)+j )  = cui( 1:3, ix(i)+j )  + w * u(1:3,i)
          dcui( 1:2, ix(i)+j ) = dcui( 1:2, ix(i)+j ) + w * du(1:2)
          amui( 1:3, ix(i)+j ) = amui( 1:3, ix(i)+j ) + w * u2(1:3)
        enddo

      enddo

      pp = pp + 1
    enddo

  enddo

  if ( noff == 0 ) then

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
        cur(1:2,1)  = 8.0 * cur(1:2,1); cur(3,1) = 0.0
        dcur(1:2,1) = 8.0 * dcur(1:2,1)
        amur(1:3,1) = 0.0
        cui(1:2,1)  = 8.0 * cui(1:2,1); cui(3,1) = 0.0
        dcui(1:2,1) = 8.0 * dcui(1:2,1)
        amui(1:3,1) = 0.0
      elseif ( mode == 2 ) then
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

end subroutine amjdeposit_clamp_part2d

subroutine interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, x, dr, bp, ep, np, ptrcur, &
  geom, weight, ix, pcos, psin )

   implicit none

   type(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im, bf_re, bf_im
   integer, intent(in) :: max_mode, np, geom
   real, intent(in) :: dr
   real, dimension(:,:), intent(in) :: x
   real, dimension(:,:), intent(inout) :: bp, ep
   integer(kind=LG), intent(in) :: ptrcur
   real, dimension(:,:), intent(out), optional :: weight
   integer, dimension(:), intent(out), optional :: ix
   real, dimension(:), intent(out), optional :: pcos, psin

   real, dimension(:,:), pointer :: e0, b0, er, ei, br, bi
   integer :: noff, i, j, nn, mode
   integer(kind=LG) :: pp
   real :: pos, idr, ph_r, ph_i, cc, ss
   real, dimension(0:1) :: wt
   complex(kind=DB) :: phase0, phase

   idr = 1.0 / dr
   noff = ef_re(0)%get_noff(1)

   e0 => ef_re(0)%get_f1()
   b0 => bf_re(0)%get_f1()

   ep = 0.0
   bp = 0.0

   pp = ptrcur
   do i = 1, np
      pos = sqrt( x(1,pp)**2 + x(2,pp)**2 ) * idr
      ! cosine and sine
      cc = x(1,pp) / pos * idr
      ss = x(2,pp) / pos * idr
      phase0 = cmplx( cc, ss )

      ! in-cell position
      nn  = int( pos )
      pos = pos - real(nn)

      ! cell index
      nn = nn - noff + 1

      ! get interpolation weight factor
      call spline_linear( pos, wt )

      if ( present(weight) ) weight(:,i) = wt
      if ( present(ix) ) ix(i) = nn
      if ( present(pcos) ) then
        pcos(i) = cc
        psin(i) = ss
      endif

      ! interpolate m=0 mode
      do j = 0, 1
        ep(:,i) = ep(:,i) + e0(:,nn+j) * wt(j)
        bp(:,i) = bp(:,i) + b0(:,nn+j) * wt(j)
      enddo

      ! interpolate m>0 modes
      phase = cmplx( 1.0, 0.0 )
      do mode = 1, max_mode
         phase = phase * phase0
         ph_r = 2.0 * real(phase)
         ph_i = 2.0 * aimag(phase)

         er => ef_re(mode)%get_f1()
         ei => ef_im(mode)%get_f1()
         br => bf_re(mode)%get_f1()
         bi => bf_im(mode)%get_f1()

         do j = 0, 1
            ep(:,i) = ep(:,i) + ( er(:,nn+j) * ph_r - ei(:,nn+j) * ph_i ) * wt(j)
            bp(:,i) = bp(:,i) + ( br(:,nn+j) * ph_r - bi(:,nn+j) * ph_i ) * wt(j)
         enddo
      enddo

      ! transform from cylindrical geometry to Cartesian geometry
      if ( geom == p_cartesian ) then
        ! ph_r, ph_i are temporary variables here
        ph_r = ep(1,i) * cc - ep(2,i) * ss
        ph_i = ep(1,i) * ss + ep(2,i) * cc
        ep(1,i) = ph_r
        ep(2,i) = ph_i

        ph_r = bp(1,i) * cc - bp(2,i) * ss
        ph_i = bp(1,i) * ss + bp(2,i) * cc
        bp(1,i) = ph_r
        bp(2,i) = ph_i
      endif

      pp = pp + 1
    enddo

end subroutine interp_emf_part2d

subroutine interp_ef_part2d( ef_re, ef_im, max_mode, x, dr, ep, np, ptrcur, &
  geom, weight, ix, pcos, psin )

  implicit none

   type(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im
   integer, intent(in) :: max_mode, np, geom
   real, intent(in) :: dr
   real, dimension(:,:), intent(in) :: x
   real, dimension(:,:), intent(inout) :: ep
   integer(kind=LG), intent(in) :: ptrcur
   real, dimension(:,:), intent(out), optional :: weight
   integer, dimension(:), intent(out), optional :: ix
   real, dimension(:), intent(out), optional :: pcos, psin

   real, dimension(:,:), pointer :: e0, er, ei
   integer :: noff, i, j, nn, mode
   integer(kind=LG) :: pp
   real :: pos, idr, ph_r, ph_i, cc, ss
   real, dimension(0:1) :: wt
   complex(kind=DB) :: phase0, phase

   idr = 1.0 / dr
   noff = ef_re(0)%get_noff(1)

   e0 => ef_re(0)%get_f1()
  
   ep = 0.0

   pp = ptrcur
   do i = 1, np
      pos = sqrt( x(1,pp)**2 + x(2,pp)**2 ) * idr
      ! cosine and sine
      cc = x(1,pp) / pos * idr
      ss = x(2,pp) / pos * idr
      phase0 = cmplx( cc, ss )

      ! in-cell position
      nn  = int( pos )
      pos = pos - real(nn)

      ! cell index
      nn = nn - noff + 1

      ! get interpolation weight factor
      call spline_linear( pos, wt )

      if ( present(weight) ) weight(:,i) = wt
      if ( present(ix) ) ix(i) = nn
      if ( present(pcos) ) then
        pcos(i) = cc
        psin(i) = ss
      endif

      ! interpolate m=0 mode
      do j = 0, 1
        ep(:,i) = ep(:,i) + e0(:,nn+j) * wt(j)
      enddo

      ! interpolate m>0 modes
      phase = cmplx( 1.0, 0.0 )
      do mode = 1, max_mode
         phase = phase * phase0
         ph_r = 2.0 * real(phase)
         ph_i = 2.0 * aimag(phase)

         er => ef_re(mode)%get_f1()
         ei => ef_im(mode)%get_f1()

         do j = 0, 1
            ep(:,i) = ep(:,i) + ( er(:,nn+j) * ph_r - ei(:,nn+j) * ph_i ) * wt(j)
         enddo
      enddo

      ! transform from cylindrical geometry to Cartesian geometry
      if ( geom == p_cartesian ) then
        ! ph_r, ph_i are temporary variables here
        ph_r = ep(1,i) * cc - ep(2,i) * ss
        ph_i = ep(1,i) * ss + ep(2,i) * cc
        ep(1,i) = ph_r
        ep(2,i) = ph_i

      endif

      pp = pp + 1
    enddo

end subroutine interp_ef_part2d

subroutine ionize_part2d( this, prof, ef, wp, dt, adk_coef )

   implicit none
   
   class(part2d), intent(inout) :: this
   class(fdist2d), intent(inout) :: prof
   double precision, dimension(:,:), intent(in) :: adk_coef
   class(field), intent(in) :: ef
   real, intent(in) :: wp, dt
  ! local data
   character(len=18), save :: sname = 'ionize_neutral'
   type(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
   real, dimension(:,:), pointer :: ef0 => null()
   real :: rn, dr, w_ion, cons, eff,esum,r
   integer :: noff, i, j, k, l, nn, mode, max_mode, np, nrp
   real, dimension(p_cache_size) :: cc, ss
   real, dimension(p_cache_size) :: pcos, psin
   integer(kind=LG) :: ptrcur, pp
   real :: idr
   integer, dimension(p_cache_size) :: ix
   real, dimension(0:1, p_cache_size) :: wt
   real, dimension(p_p_dim, p_cache_size) :: ep
   complex(kind=DB) :: phase0, phase
   
  idr = 1.0 / this%dr

  ef_re => ef%get_rf_re()
  ef_im => ef%get_rf_im()

  max_mode = ef%get_max_mode()

  noff = ef_re(0)%get_noff(1)
  nrp  = ef_re(0)%get_ndp(1)
  ef0  => ef_re(0)%get_f1()
  esum = 0.0
  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    call interp_ef_part2d( ef_re, ef_im, max_mode, this%x, this%dr, ep, np, &
      ptrcur, p_cylindrical)
!     write(2,*) wp, "ionize_part2d"
    pp = ptrcur
      do i = 1, np
        eff = sqrt(ep(1,i)**2 + ep(2,i)**2 + ep(3,i)**2)
        eff = eff * wp * 1.708e-12
!         write(2,*) eff, "ionize_part2d"
        if ((eff .gt. 1.0e-6) .and. (this%q(pp) .gt. 1.0e-6 )) then 
          w_ion = adk_coef(1,1)*eff**(-adk_coef(3,1))*exp(-adk_coef(2,1)/eff&
            &)/wp 

                ! w_ion is in normalized unit now
                !2nd rk (dn/dt==n*w_ion)
                !>>>>>> w_ion* dt/(1.0 + 0.5 * w_ion * dt) could be larger than 1. 
                !>>>>>> should add:
          cons = w_ion* dt
          cons = cons/(1.0 + 0.5 * cons)
          cons = cons*this%gamma(pp)/(this%gamma(pp) - this%p(3,pp))
          this%w(pp) = (1.0-this%w(pp))*cons+this%w(pp)
          ! write(2,*) "ionize_part2d",this%w(pp),dt
          if (this%w(pp) .gt. 1.0) this%w(pp) = 1.0
        endif
!         r = sqrt(this%x(1,pp)**2 + this%x(2,pp)**2)
!         if ((r .gt. 4.5).and.(r .lt. 5)) then
! !           write(2,*) w_ion, "w_ion"
!           write(2,*) eff, "eff"
!           write(2,*) w_ion, "w_ion"
!           write(2,*) r, "particle_r"
!         endif
        pp = pp + 1
        esum = esum + eff
      enddo
  enddo
  write(2,*) esum, "ionize_part2d_esum"
end subroutine ionize_part2d

subroutine add_particles_part2d( this, prof, ppart1, ppart2, multi_max, m, sec, s)

    implicit none
    class(part2d), intent(inout) :: this
    class(fdist2d), intent(inout) :: prof
    class(part2d), intent(inout) :: ppart1
    class(part2d), intent(inout) :: ppart2
    real, intent(in) :: s, sec
    integer, intent(in) :: multi_max, m

    ! local
    character(len=18), save :: sname = 'add_particles'
    integer :: noff, i, j, k, pp1, pp2, m1, m2, np, l
    integer(kind=LG) :: ptrcur, pp
    real :: dxp, ww, a, b

    call write_dbg( cls_name, sname, cls_level, 'starts' )

      if ( m .le. 2 ) then
        m1 = 1
        m2 = m
      else
        m1 = m - 1
        m2 = m
      end if
      write(2,*) m1,m2,"m1,m2"
      pp1 = ppart1%npp
      pp2 = ppart2%npp
!       write(2,*) pp2, "particles_e"
     do ptrcur = 1, this%npp, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > this%npp ) then
          np = this%npp - ptrcur + 1
        else
          np = p_cache_size
        endif

        pp = ptrcur
        do i = 1, np
          !Adding particles in segments according to the gain of ionization rate
          if (sec.ge.2) then
            ww = this%w(pp) - this%w0(pp)
            if ( this%w(pp) .ge. 0.95) then

                  dxp = this%q(pp)/m1
                  pp1 = pp1 + 1
                  pp2 = pp2 + 1
                  ppart1%x(1,pp1)   = this%x(1,pp)
                  ppart1%x(2,pp1)   = this%x(2,pp)
                  ppart1%q(pp1)     = dxp*m2
                  ppart1%p(1,pp1)   = this%p(1,pp)
                  ppart1%p(2,pp1)   = this%p(2,pp)
                  ppart1%p(3,pp1)   = this%p(3,pp)
                  ppart1%gamma(pp1) = 1.0
                  ppart1%psi(pp1) = 1.0
                  ppart1%w(pp1) = 0.0
                  ppart1%w0(pp1) = 0.0

                  ppart2%x(1,pp2)   = this%x(1,pp)
                  ppart2%x(2,pp2)   = this%x(2,pp)
                  ppart2%q(pp2)     = -dxp
                  ppart2%p(1,pp2)   = this%p(1,pp)
                  ppart2%p(2,pp2)   = this%p(2,pp)
                  ppart2%p(3,pp2)   = this%p(3,pp)
                  ppart2%gamma(pp2) = 1.0
                  ppart2%psi(pp2) = 1.0
                  ppart2%w(pp2) = 0.0
                  ppart2%w0(pp2) = 0.0

                  this%q(pp) = 0.0
                  this%w0(pp) = 1.0
                  this%w(pp) = 0.0

            else

                do l = 1, sec - 1
                  b = l/sec
                  a = 1 - b
                  if( ww .ge. a) then
                    dxp = this%q(pp)/m1
                    pp1 = pp1 + 1
                    pp2 = pp2 + 1

                    ppart1%x(1,pp1)   = this%x(1,pp)
                    ppart1%x(2,pp1)   = this%x(2,pp)
                    ppart1%q(pp1)     = dxp*m2*a/(1-this%w0(pp))
  !                   ppart1%q(pp1)     = dxp*m2*0.75
                    ppart1%p(1,pp1)   = this%p(1,pp)
                    ppart1%p(2,pp1)   = this%p(2,pp)
                    ppart1%p(3,pp1)   = this%p(3,pp)
                    ppart1%gamma(pp1) = 1.0
                    ppart1%psi(pp1) = 1.0
                    ppart1%w(pp1) = 0.0
                    ppart1%w0(pp1) = 0.0

                    ppart2%x(1,pp2)   = this%x(1,pp)
                    ppart2%x(2,pp2)   = this%x(2,pp)
                    ppart2%q(pp2)     = -dxp*a/(1-this%w0(pp))
                    ppart2%p(1,pp2)   = this%p(1,pp)
                    ppart2%p(2,pp2)   = this%p(2,pp)
                    ppart2%p(3,pp2)   = this%p(3,pp)
                    ppart2%gamma(pp2) = 1.0
                    ppart2%psi(pp2) = 1.0
                    ppart2%w(pp2) = 0.0
                    ppart2%w0(pp2) = 0.0

                    this%q(pp) = this%q(pp)*((1-this%w0(pp)-a)/(1-this%w0(pp)))
                    this%w0(pp) = a + this%w0(pp)
                    exit
                  endif
                enddo
            endif
          elseif (sec.ge.1) then
            if ( this%w(pp) .ge. 0.95) then

                  dxp = this%q(pp)/m1
                  pp1 = pp1 + 1
                  pp2 = pp2 + 1
                  ppart1%x(1,pp1)   = this%x(1,pp)
                  ppart1%x(2,pp1)   = this%x(2,pp)
                  ppart1%q(pp1)     = dxp*m2
                  ppart1%p(1,pp1)   = this%p(1,pp)
                  ppart1%p(2,pp1)   = this%p(2,pp)
                  ppart1%p(3,pp1)   = this%p(3,pp)
                  ppart1%gamma(pp1) = 1.0
                  ppart1%psi(pp1) = 1.0
                  ppart1%w(pp1) = 0.0
                  ppart1%w0(pp1) = 0.0

                  ppart2%x(1,pp2)   = this%x(1,pp)
                  ppart2%x(2,pp2)   = this%x(2,pp)
                  ppart2%q(pp2)     = -dxp
                  ppart2%p(1,pp2)   = this%p(1,pp)
                  ppart2%p(2,pp2)   = this%p(2,pp)
                  ppart2%p(3,pp2)   = this%p(3,pp)
                  ppart2%gamma(pp2) = 1.0
                  ppart2%psi(pp2) = 1.0
                  ppart2%w(pp2) = 0.0
                  ppart2%w0(pp2) = 0.0

                  this%q(pp) = 0.0
                  this%w0(pp) = 1.0
                  this%w(pp) = 0.0
            endif
          endif

          pp = pp + 1
        enddo
     enddo
     ppart1%npp = pp1
     ppart2%npp = pp2
!     write(2,*) pp1, "add_particles"
    write(2,*) pp1, "add_particles_ion"
    write(2,*) pp2, "add_particles_e"
    call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_particles_part2d


subroutine push_robust_part2d( this, ef, bf )

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: ef, bf
  ! local data
  character(len=18), save :: sname = 'push_robust_part2d'
  type(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im

  integer :: i, np, max_mode
  real :: qtmh, qtmh1, qtmh2, gam, dtc, ostq
  real, dimension(p_p_dim, p_cache_size) :: bp, ep, utmp
  integer(kind=LG) :: ptrcur, pp

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'push 2D particles' )

  qtmh = this%qbm * this%dt * 0.5
  max_mode = ef%get_max_mode()

  ef_re => ef%get_rf_re()
  ef_im => ef%get_rf_im()
  bf_re => bf%get_rf_re()
  bf_im => bf%get_rf_im()

  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! interpolate fields to particles
    call interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, this%x, this%dr, &
      bp, ep, np, ptrcur, p_cartesian )

    pp = ptrcur
    do i = 1, np
      gam = sqrt( 1.0 + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
      ! qtmh1 = qtmh / this%psi(pp)
      ! qtmh2 = qtmh1 * this%gamma(pp)
      qtmh1 = qtmh / ( gam - this%p(3,pp) )
      qtmh2 = qtmh1 * gam
      ep(:,i) = ep(:,i) * qtmh2
      bp(:,i) = bp(:,i) * qtmh1
      pp = pp + 1
    enddo

    ! first half of electric field acceleration
    pp = ptrcur
    do i = 1, np
      utmp(:,i) = this%p(:,pp) + ep(:,i)
      pp = pp + 1
    enddo

    ! rotation about magnetic field
    pp = ptrcur
    do i = 1, np
      this%p(1,pp) = utmp(1,i) + utmp(2,i) * bp(3,i) - utmp(3,i) * bp(2,i)
      this%p(2,pp) = utmp(2,i) + utmp(3,i) * bp(1,i) - utmp(1,i) * bp(3,i)
      this%p(3,pp) = utmp(3,i) + utmp(1,i) * bp(2,i) - utmp(2,i) * bp(1,i)
      pp = pp + 1
    enddo

    do i = 1, np
      ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
      bp(1,i) = bp(1,i) * ostq
      bp(2,i) = bp(2,i) * ostq
      bp(3,i) = bp(3,i) * ostq
    enddo

    pp = ptrcur
    do i = 1, np
      utmp(1,i) = utmp(1,i) + this%p(2,pp) * bp(3,i) - this%p(3,pp) * bp(2,i)
      utmp(2,i) = utmp(2,i) + this%p(3,pp) * bp(1,i) - this%p(1,pp) * bp(3,i)
      utmp(3,i) = utmp(3,i) + this%p(1,pp) * bp(2,i) - this%p(2,pp) * bp(1,i)
      pp = pp + 1
    enddo

    ! second half of electric field acc.
    pp = ptrcur
    do i = 1, np
      this%p(:,pp) = utmp(:,i) + ep(:,i)
      pp = pp + 1
    enddo

    ! advance particle position
    pp = ptrcur
    do i = 1, np
      gam = sqrt( 1.0 + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
      dtc = this%dt / ( gam - this%p(3,pp) )
      this%x(1,pp) = this%x(1,pp) + this%p(1,pp) * dtc
      this%x(2,pp) = this%x(2,pp) + this%p(2,pp) * dtc
      pp = pp + 1
    enddo

  enddo

  call stop_tprof( 'push 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_robust_part2d

subroutine push_robust_subcyc_part2d( this, ef, bf )

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: ef, bf
  ! local data
  character(len=18), save :: sname = 'push_robust_subcyc_part2d'
  type(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im

  integer :: i, np, max_mode, exit_cnt
  integer :: n_subcyc, n_subcyc_max, np_fail
  real :: qtmh1, qtmh2, dtc, ostq
  real, dimension(p_p_dim, p_cache_size) :: bp, ep, p_subcyc
  real, dimension(2, p_cache_size) :: x_subcyc
  real, dimension(p_cache_size) :: dt_subcyc, gam_subcyc
  real, dimension(p_p_dim) :: utmp
  integer, dimension(p_cache_size) :: ndt_rem
  integer(kind=LG), dimension(p_cache_size) :: ix_subcyc
  integer(kind=LG) :: ptrcur, pp

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'push 2D particles' )

  n_subcyc_max = 0
  np_fail = 0
  max_mode = ef%get_max_mode()

  ef_re => ef%get_rf_re()
  ef_im => ef%get_rf_im()
  bf_re => bf%get_rf_re()
  bf_im => bf%get_rf_im()

  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! initialize sub-cycling
    pp = ptrcur
    do i = 1, np
      p_subcyc(:,i) = this%p(:,pp)
      
      ! synchronize the particle position with the momentum for sub-cycling
      gam_subcyc(i) = sqrt( 1.0 + p_subcyc(1,i)**2 + p_subcyc(2,i)**2 + p_subcyc(3,i)**2 )
      dtc = 0.5 * this%dt / ( gam_subcyc(i) - p_subcyc(3,i) )
      x_subcyc(1,i) = this%x(1,pp) - dtc * p_subcyc(1,i)
      x_subcyc(2,i) = this%x(2,pp) - dtc * p_subcyc(2,i)

      ! particle index for sub-cycling
      ix_subcyc(i) = pp

      ! calculate the current sub-cycling time step
      ndt_rem(i) = ceiling( this%dt / min( this%dt, this%dt_eff_max * (1.0 - p_subcyc(3,i) / gam_subcyc(i)) ) )
      dt_subcyc(i) = this%dt / ndt_rem(i)

      ! advance position by half sub-cycling time step
      dtc = 0.5 * dt_subcyc(i) / ( gam_subcyc(i) - p_subcyc(3,i) )
      x_subcyc(1,i) = x_subcyc(1,i) + p_subcyc(1,i) * dtc
      x_subcyc(2,i) = x_subcyc(2,i) + p_subcyc(2,i) * dtc

      pp = pp + 1
    enddo

    n_subcyc = 0
    ! begin sub-cycling
    do while( np > 0 )

      ! interpolate fields to particles for sub-cycling
      call interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, x_subcyc, this%dr, &
        bp, ep, np, int(1, kind=LG), p_cartesian )
      
      do i = 1, np

        qtmh1 = 0.5 * dt_subcyc(i) * this%qbm / ( gam_subcyc(i) - p_subcyc(3,i) )
        qtmh2 = qtmh1 * gam_subcyc(i)

        ep(:,i) = ep(:,i) * qtmh2
        bp(:,i) = bp(:,i) * qtmh1

        ! first half of electric field acceleration
        utmp = p_subcyc(:,i) + ep(:,i)

        ! rotation about magnetic field
        p_subcyc(1,i) = utmp(1) + utmp(2) * bp(3,i) - utmp(3) * bp(2,i)
        p_subcyc(2,i) = utmp(2) + utmp(3) * bp(1,i) - utmp(1) * bp(3,i)
        p_subcyc(3,i) = utmp(3) + utmp(1) * bp(2,i) - utmp(2) * bp(1,i)

        ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
        bp(1,i) = bp(1,i) * ostq
        bp(2,i) = bp(2,i) * ostq
        bp(3,i) = bp(3,i) * ostq

        utmp(1) = utmp(1) + p_subcyc(2,i) * bp(3,i) - p_subcyc(3,i) * bp(2,i)
        utmp(2) = utmp(2) + p_subcyc(3,i) * bp(1,i) - p_subcyc(1,i) * bp(3,i)
        utmp(3) = utmp(3) + p_subcyc(1,i) * bp(2,i) - p_subcyc(2,i) * bp(1,i)

        ! second half of electric field acceleration
        p_subcyc(:,i) = utmp + ep(:,i)

        ! update the remaining time step numbers
        ndt_rem(i) = ndt_rem(i) - 1

        ! advance position
        gam_subcyc(i) = sqrt( 1.0 + p_subcyc(1,i)**2 + p_subcyc(2,i)**2 + p_subcyc(3,i)**2 )
        ! dt_subcyc(i) = min( dt_rem(i), this%dt_eff_max * (1.0 - p_subcyc(3,i) / gam_subcyc(i)) )
        dtc = dt_subcyc(i) / ( gam_subcyc(i) - p_subcyc(3,i) )
        x_subcyc(1,i) = x_subcyc(1,i) + p_subcyc(1,i) * dtc
        x_subcyc(2,i) = x_subcyc(2,i) + p_subcyc(2,i) * dtc

      enddo

      ! store the particles that finished sub-cycling and rearrange sub-cycling array
      exit_cnt = 0; i = 1
      do while ( i <= np - exit_cnt )

        ! the particles that have used up the global 2D time step or went out of the boundary
        ! will exit the sub-cycling.
        if ( ndt_rem(i) == 0 .or. &
          x_subcyc(1,i)**2 + x_subcyc(2,i)**2 >= this%edge**2 ) then

          ! store the momentum          
          this%p( :, ix_subcyc(i) ) = p_subcyc(:,i)

          ! advance the position by a full 2D time step and store it
          dtc = this%dt / ( sqrt( 1.0 + p_subcyc(1,i)**2 + p_subcyc(2,i)**2 + p_subcyc(3,i)**2 ) - p_subcyc(3,i) )
          this%x( 1, ix_subcyc(i) ) = this%x( 1, ix_subcyc(i) ) + p_subcyc(1,i) * dtc
          this%x( 2, ix_subcyc(i) ) = this%x( 2, ix_subcyc(i) ) + p_subcyc(2,i) * dtc

          ! move the last element of sub-cycling array to current position
          x_subcyc(:,i) = x_subcyc( :, np - exit_cnt )
          p_subcyc(:,i) = p_subcyc( :, np - exit_cnt )
          ix_subcyc(i)  = ix_subcyc( np - exit_cnt )
          dt_subcyc(i)  = dt_subcyc( np - exit_cnt )
          gam_subcyc(i) = gam_subcyc( np - exit_cnt )
          ndt_rem(i)    = ndt_rem( np - exit_cnt )

          exit_cnt = exit_cnt + 1

        else
          i = i + 1
        endif

      enddo
      np = np - exit_cnt
      n_subcyc = n_subcyc + 1

      ! when reaching the max number of sub-cycling, store the unfinished particles
      ! and pop out a warning.
      if ( n_subcyc == p_max_subcyc ) then
        np_fail = np_fail + np
        do i = 1, np
          this%p(:, ix_subcyc(i) ) = p_subcyc(:,i)
          dtc = this%dt / ( sqrt( 1.0 + p_subcyc(1,i)**2 + p_subcyc(2,i)**2 + p_subcyc(3,i)**2 ) - p_subcyc(3,i) )
          this%x( 1, ix_subcyc(i) ) = this%x( 1, ix_subcyc(i) ) + p_subcyc(1,i) * dtc
          this%x( 2, ix_subcyc(i) ) = this%x( 2, ix_subcyc(i) ) + p_subcyc(2,i) * dtc
        enddo
        exit
      endif

    enddo ! sub-cycling

    ! DEBUG
    n_subcyc_max = max(n_subcyc_max, n_subcyc)

  enddo ! chunk loop

  ! DEBUG
  ! call write_stdout( "[push] max sub-cycling times = " // num2str(n_subcyc_max) )
  if ( np_fail > 0 ) then
    call write_stdout( '[push] Max number of sub-cycling reached. ' // num2str(np_fail) // &
          ' particles have not yet finished sub-cycling.' )
  endif

  call stop_tprof( 'push 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_robust_subcyc_part2d

subroutine push_clamp_part2d( this, ef, bf )

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: ef, bf
  ! local data
  character(len=18), save :: sname = 'push_clamp_part2d'
  type(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im

  integer :: i, np, max_mode
  real :: qtmh, qtmh1, qtmh2, gam, dtc, ostq, psi_plus1
  real, dimension(p_p_dim, p_cache_size) :: bp, ep, utmp
  integer(kind=LG) :: ptrcur, pp

  integer :: stat

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'push 2D particles' )

  qtmh = this%qbm * this%dt * 0.5
  max_mode = ef%get_max_mode()

  ef_re => ef%get_rf_re()
  ef_im => ef%get_rf_im()
  bf_re => bf%get_rf_re()
  bf_im => bf%get_rf_im()

  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! interpolate fields to particles
    call interp_emf_part2d( ef_re, ef_im, bf_re, bf_im, max_mode, this%x, this%dr, &
      bp, ep, np, ptrcur, p_cartesian )

    pp = ptrcur
    do i = 1, np
      ! Clamp the value of gamma/(1 + psi). It should be already checked in amjdeposit.
      ! Here check it again for safety
      call clamp_momentum( this%fac_clamp, this%p(:,pp), this%gamma(pp), this%psi(pp) )

      qtmh1 = qtmh / this%psi(pp)
      qtmh2 = qtmh1 * this%gamma(pp)
      ep(:,i) = ep(:,i) * qtmh2
      bp(:,i) = bp(:,i) * qtmh1
      pp = pp + 1
    enddo

    ! first half of electric field acceleration
    pp = ptrcur
    do i = 1, np
      utmp(:,i) = this%p(:,pp) + ep(:,i)
      pp = pp + 1
    enddo

    ! rotation about magnetic field
    pp = ptrcur
    do i = 1, np
      this%p(1,pp) = utmp(1,i) + utmp(2,i) * bp(3,i) - utmp(3,i) * bp(2,i)
      this%p(2,pp) = utmp(2,i) + utmp(3,i) * bp(1,i) - utmp(1,i) * bp(3,i)
      this%p(3,pp) = utmp(3,i) + utmp(1,i) * bp(2,i) - utmp(2,i) * bp(1,i)
      pp = pp + 1
    enddo

    do i = 1, np
      ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
      bp(1,i) = bp(1,i) * ostq
      bp(2,i) = bp(2,i) * ostq
      bp(3,i) = bp(3,i) * ostq
    enddo

    pp = ptrcur
    do i = 1, np
      utmp(1,i) = utmp(1,i) + this%p(2,pp) * bp(3,i) - this%p(3,pp) * bp(2,i)
      utmp(2,i) = utmp(2,i) + this%p(3,pp) * bp(1,i) - this%p(1,pp) * bp(3,i)
      utmp(3,i) = utmp(3,i) + this%p(1,pp) * bp(2,i) - this%p(2,pp) * bp(1,i)
      pp = pp + 1
    enddo

    ! second half of electric field acc.
    pp = ptrcur
    do i = 1, np
      this%p(:,pp) = utmp(:,i) + ep(:,i)
      pp = pp + 1
    enddo

    ! advance particle position
    pp = ptrcur
    do i = 1, np
      gam = sqrt( 1.0 + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
      psi_plus1 = gam - this%p(3,pp)

      ! clamp the updated gamma / (1 + psi)
      call clamp_momentum( this%fac_clamp, this%p(:,pp), gam, psi_plus1 )

      dtc = this%dt / psi_plus1
      this%x(1,pp) = this%x(1,pp) + this%p(1,pp) * dtc
      this%x(2,pp) = this%x(2,pp) + this%p(2,pp) * dtc
      pp = pp + 1
    enddo

  enddo

  call stop_tprof( 'push 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_clamp_part2d

subroutine clamp_momentum( fac_clamp, p, gam, psi_plus1 )

  implicit none

  real, intent(in) :: fac_clamp
  real, intent(inout), dimension(3) :: p
  real, intent(inout) :: gam, psi_plus1

  real :: scale

  if ( gam / psi_plus1 > fac_clamp ) then

    scale = ( fac_clamp - 1.0 )**2 * ( 1.0 + p(1)**2 + p(2)**2 )
    scale = scale / ( (2.0 * fac_clamp - 1.0) * p(3)**2 )
    scale = sqrt(scale)

    p(3) = p(3) * scale
    gam = sqrt( 1.0 + p(1)**2 + p(2)**2 + p(3)**2 )
    psi_plus1 = gam - p(3)

  endif

end subroutine clamp_momentum

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
         this%w(i)     = this%w(this%npp)
         this%w0(i)    = this%w0(this%npp)
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

! subroutine extract_psi_part2d(this,psi)

!    implicit none

!    class(part2d), intent(inout) :: this
!    class(field), intent(in) :: psi
! ! local data
!    character(len=18), save :: sname = 'extractpsi'
!    class(ufield), dimension(:), pointer :: psi_re => null(), psi_im => null()

!    call write_dbg(cls_name, sname, cls_level, 'starts')

!    psi_re => psi%get_rf_re()
!    psi_im => psi%get_rf_im()
!    ! call part2d_extractpsi(this%part,this%npp,this%dr,this%qbm,psi_re,psi_im,psi%get_max_mode())
   
!    call part2d_extractpsi(this%x, this%p, this%gamma, this%q, this%psi,&
!     this%npp,this%dr,this%qbm,psi_re,psi_im,psi%get_max_mode())

!    call write_dbg(cls_name, sname, cls_level, 'ends')

! end subroutine extract_psi_part2d

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

  ! rearrange the particle 1 + psi
  this%pbuf(1:this%npp) = this%psi( 1:this%npp )
  do i = 1, this%npp
    this%psi( sort_idx(i) ) = this%pbuf(i)
  enddo

  deallocate( ix )

  call stop_tprof( 'sort 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')
    
end subroutine sort_part2d

end module part2d_class