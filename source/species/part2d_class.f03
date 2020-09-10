! part2d class for QPAD

module part2d_class

use param
use sysutil
use parallel_pipe_class
use grid_class
use field_class
use ufield_class
use fdist2d_class
use hdf5io_class
use part2d_lib
use mpi
use interpolation
use debug_tool

implicit none

private

public :: part2d

type part2d

   class(parallel_pipe), pointer :: pp => null()

   ! qbm = particle charge/mass ratio
   ! dt = time interval between successive calculations
   ! dr = radial cell size
   real :: qbm, dt, dr

   ! nbmax = size of buffer for passing particles between processors
   ! npp = number of particles in current partition
   ! npmax = maximum number of particles in each partition
   integer(kind=LG) :: npmax, nbmax, npp

   ! dimension of particle coordinates
   integer :: part_dim

   ! array for particle position
   real, dimension(:,:), pointer :: x => null()
   ! array for particle momenta
   real, dimension(:,:), pointer :: p => null()
   ! array for time-centered gamma
   real, dimension(:), pointer :: gamma => null()
   ! array for particle charge
   real, dimension(:), pointer :: q => null()
   ! array for psi
   real, dimension(:), pointer :: psi => null()

   ! particle upper boundaries
   real :: edge

   contains

   procedure :: new          => init_part2d
   procedure :: renew        => renew_part2d
   procedure :: del          => end_part2d
   procedure :: qdeposit     => qdeposit_part2d
   procedure :: amjdeposit   => amjdeposit_part2d
   procedure :: push         => push_part2d
   procedure :: update_bound => update_bound_part2d
   ! procedure :: extract_psi  => extract_psi_part2d
   procedure :: pipesend     => pipesend_part2d
   procedure :: piperecv     => piperecv_part2d

   ! TODO: particle manager to be rewritten
   ! procedure :: pmv => pmove

   procedure :: wr => writehdf5_part2d

end type

save

character(len=20), parameter :: cls_name = "part2d"
integer, parameter :: cls_level = 2

! TODO: data communication, to be deleted
real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
integer(kind=LG), dimension(:), allocatable :: ihole

contains

subroutine init_part2d( this, pp, gp, pf, qbm, dt, s )

   implicit none

   class(part2d), intent(inout) :: this
   class(parallel_pipe), intent(in), pointer :: pp
   class(grid), intent(in), pointer :: gp
   class(fdist2d), intent(inout) :: pf
   real, intent(in) :: qbm, dt, s

   ! local data
   character(len=18), save :: sname = 'init_part2d'
   integer :: npmax, nbmax

   call write_dbg( cls_name, sname, cls_level, 'starts' )

   this%pp  => pp
   this%qbm = qbm
   this%dt  = dt
   this%part_dim = 2 + p_p_dim + 3

   npmax      = pf%getnpmax()
   nbmax      = max(int(0.01*npmax),100)
   this%npmax = npmax
   this%nbmax = nbmax
   this%npp   = 0

   this%dr   = pf%getdex()
   this%edge = gp%get_nd(1) * this%dr
   
   allocate( this%x( 2, npmax ) )
   allocate( this%p( p_p_dim, npmax ) )
   allocate( this%gamma( npmax ), this%q( npmax ), this%psi( npmax ) )

   ! initialize particle coordinates according to specified profile
   call pf%dist( this%x, this%p, this%gamma, this%q, this%psi, this%npp, s )
   if (.not. allocated(sbufl)) then
      allocate( sbufl( this%part_dim, nbmax ), sbufr( this%part_dim, nbmax ) )
      allocate( rbufl( this%part_dim, nbmax ), rbufr( this%part_dim, nbmax ) )
      allocate( ihole( nbmax*2 ) )
   end if

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

subroutine renew_part2d( this, pf, s )

   implicit none

   class(part2d), intent(inout) :: this
   class(fdist2d), intent(inout) :: pf
   real, intent(in) :: s
   ! local data
   character(len=18), save :: sname = 'renew_part2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   call pf%dist( this%x, this%p, this%gamma, this%q, this%psi, this%npp, s )

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

  max_mode = q%get_num_modes()

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

subroutine amjdeposit_part2d( this, ef, bf, cu, amu, dcu )
! deposit the current, acceleration and momentum flux

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: cu, amu, dcu
  class(field), intent(in) :: ef, bf
  ! local data
  character(len=18), save :: sname = 'amjdeposit_part2d'
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
  max_mode = ef%get_num_modes()

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

end subroutine amjdeposit_part2d

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

subroutine push_part2d( this, ef, bf )

  implicit none

  class(part2d), intent(inout) :: this
  class(field), intent(in) :: ef, bf
  ! local data
  character(len=18), save :: sname = 'push_part2d'
  type(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im

  integer :: i, np, max_mode
  real :: qtmh, qtmh1, qtmh2, gam, dtc, ostq
  real, dimension(p_p_dim, p_cache_size) :: bp, ep, utmp
  integer(kind=LG) :: ptrcur, pp

  integer :: stat

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'push 2D particles' )

  qtmh = this%qbm * this%dt * 0.5
  max_mode = ef%get_num_modes()

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
      dtc = this%dt / ( gam - this%p(3,pp) )
      this%x(1,pp) = this%x(1,pp) + this%p(1,pp) * dtc
      this%x(2,pp) = this%x(2,pp) + this%p(2,pp) * dtc
      pp = pp + 1
    enddo

  enddo

  call stop_tprof( 'push 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_part2d

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

! subroutine pmove(this,fd)

!    implicit none

!    class(part2d), intent(inout) :: this
!    class(field), intent(in) :: fd
! ! local data
!    character(len=18), save :: sname = 'pmove'
!    class(ufield), pointer :: ud

!    integer :: stat

!    call write_dbg(cls_name, sname, cls_level, 'starts')

!    ud => fd%get_rf_re(0)
!    ! call part2d_pmove(this%part,this%pp,this%npp,this%dr,this%xdim,this%npmax,&
!    ! &this%nbmax,ud,sbufl,sbufr,rbufl,rbufr,ihole)

!    call part2d_pmove(this%x, this%p, this%gamma, this%q, this%psi,&
!    &this%pp,this%npp,this%dr,this%part_dim,this%npmax,&
!    &this%nbmax,ud,sbufl,sbufr,rbufl,rbufr,ihole)

!    call write_dbg(cls_name, sname, cls_level, 'ends')

! end subroutine pmove

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
!    ! call part2d_extractpsi(this%part,this%npp,this%dr,this%qbm,psi_re,psi_im,psi%get_num_modes())
   
!    call part2d_extractpsi(this%x, this%p, this%gamma, this%q, this%psi,&
!     this%npp,this%dr,this%qbm,psi_re,psi_im,psi%get_num_modes())

!    call write_dbg(cls_name, sname, cls_level, 'ends')

! end subroutine extract_psi_part2d

subroutine pipesend_part2d(this, tag, id)

  implicit none

  class(part2d), intent(inout) :: this
  integer, intent(in) :: tag
  integer, intent(inout) :: id

  ! local data
  character(len=18), save :: sname = 'pipesend_part2d'
  integer :: des, ierr, i
  real, dimension(:,:), allocatable, save :: sbuf

  call write_dbg(cls_name, sname, cls_level, 'starts')

  if ( .not. allocated(sbuf) ) then
    allocate( sbuf( this%part_dim, this%npmax ) )
  endif

  des = this%pp%getidproc() + this%pp%getlnvp()

  if (des >= this%pp%getnvp()) then
    id = MPI_REQUEST_NULL
    call write_dbg(cls_name, sname, cls_level, 'ends')
    return
  endif

  ! to be implemented if using tile
  ! call this%pcb()

  do i = 1, this%npp
    sbuf(1:2,i) = this%x(:,i)
    sbuf(3:5,i) = this%p(:,i)
    sbuf(6,i)   = this%gamma(i)
    sbuf(7,i)   = this%psi(i)
    sbuf(8,i)   = this%q(i)
  enddo

  ! NOTE: npp*xdim might be larger than MAX_INT32
  call MPI_ISEND(sbuf, int(this%npp*this%part_dim), this%pp%getmreal(), &
    des, tag, this%pp%getlworld(), id, ierr)

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
  integer, dimension(MPI_STATUS_SIZE) :: istat
  integer :: nps, des, ierr, i
  real, dimension(:,:), allocatable, save :: rbuf

  call write_dbg(cls_name, sname, cls_level, 'starts')

  if ( .not. allocated(rbuf) ) then
    allocate( rbuf( this%part_dim, this%npmax ) )
  endif

  des = this%pp%getidproc() - this%pp%getlnvp()

  if (des < 0) then
    call write_dbg(cls_name, sname, cls_level, 'ends')
    return
  endif

  ! NOTE: npp*xdim might be larger than MAX_INT32
  call MPI_RECV(rbuf, int(this%npmax*this%part_dim), this%pp%getmreal(), &
    des, tag, this%pp%getlworld(), istat, ierr)

  call MPI_GET_COUNT(istat, this%pp%getmreal(), nps, ierr)

  this%npp = nps/this%part_dim

  do i = 1, this%npp
    this%x(:,i)   = rbuf(1:2,i)
    this%p(:,i)   = rbuf(3:5,i)
    this%gamma(i) = rbuf(6,i)
    this%psi(i)   = rbuf(7,i)
    this%q(i)     = rbuf(8,i)
  enddo

  ! to be implemented if using tile
  ! call this%pcp(fd)

  ! check for errors
  if (ierr /= 0) then
    call write_err('MPI failed')
  endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine piperecv_part2d
!
subroutine writehdf5_part2d(this,file)

   implicit none

   class(part2d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'writehdf5_part2d'
   integer :: ierr

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ! call pwpart(this%pp,file,this%part,this%npp,1,ierr)
   call pwpart(this%pp,file,this%x, this%p, this%q, this%npp,1,ierr)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_part2d
!
end module part2d_class