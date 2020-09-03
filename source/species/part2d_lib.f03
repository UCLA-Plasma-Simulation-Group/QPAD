! part2d_lib module for QPAD
! update: 04/18/2016

module part2d_lib

use parallel_pipe_class
use ufield_class
use param
use sysutil
use mpi

implicit none

public
!
character(len=20), parameter :: cls_name = "part2d_lib"
integer, parameter :: cls_level = 2
character(len=128) :: erstr

private :: cls_name, cls_level, erstr

contains
!
! subroutine part2d_qdeposit(part,npp,dr,q_re,q_im,num_modes)

!    implicit none

!    real, dimension(:,:), pointer, intent(in) :: part
!    integer(kind=LG), intent(in) :: npp
!    real, intent(in) :: dr
!    class(ufield), dimension(:), pointer, intent(in) :: q_re, q_im
!    integer, intent(in) :: num_modes
! ! local data
!    character(len=20), save :: sname = "part2d_qdeposit"
!    integer :: i, nn, noff, n1p, j
!    integer(kind=LG) :: ii
!    real, dimension(:,:), pointer :: q0, qr, qi
!    real :: r0, r, qc, lq, rq, rcr, rci
!    complex(kind=DB) :: rc, rc0

!    call write_dbg(cls_name, sname, cls_level, 'starts')
!    call start_tprof( 'deposit 2D particles' )

!    noff = q_re(0)%get_noff(1)
!    n1p = q_re(0)%get_ndp(1)
!    q0 => q_re(0)%get_f1()
!    qr => null(); qi => null()

!    do ii = 1, npp
!       r0 = sqrt( part(1,ii)**2 + part(2,ii)**2 )
!       qc = part(8,ii)
!       rc0 = cmplx( part(1,ii), -part(2,ii), kind=DB) / r0
!       r = r0/dr + 1.0
!       nn = r
!       rq = qc * (r - real(nn))
!       lq = qc - rq
!       nn = nn - noff
!       q0(1,nn)   = q0(1,nn)   + lq
!       q0(1,nn+1) = q0(1,nn+1) + rq
!       rc = rc0
!       do i = 1, num_modes
!          rcr = real(rc)
!          rci = aimag(rc)
!          qr => q_re(i)%get_f1()
!          qi => q_im(i)%get_f1()
!          qr(1,nn)   = qr(1,nn)   + lq * rcr
!          qr(1,nn+1) = qr(1,nn+1) + rq * rcr
!          qi(1,nn)   = qi(1,nn)   + lq * rci
!          qi(1,nn+1) = qi(1,nn+1) + rq * rci
!          rc = rc * rc0
!       enddo
!    enddo

!    if (noff == 0) then

!       q0(1,0) = 0.0 ! guard cell is useless on axis
!       q0(1,1) = 8.0 * q0(1,1)
!       do j = 2, n1p+1
!          r = j + noff - 1
!          q0(1,j) = q0(1,j) / r
!       enddo

!       do i = 1, num_modes
!          qr => q_re(i)%get_f1()
!          qi => q_im(i)%get_f1()
!          qr(1,0) = 0.0
!          qi(1,0) = 0.0
!          qr(1,1) = 0.0
!          qi(1,1) = 0.0
!          do j = 2, n1p+1
!             r = j + noff - 1
!             qr(1,j) = qr(1,j) / r
!             qi(1,j) = qi(1,j) / r
!          enddo
!       enddo

!    else

!       do j = 0, n1p+1
!          r = j + noff -1
!          q0(1,j) = q0(1,j) / r
!       enddo

!       do i = 1, num_modes
!          qr => q_re(i)%get_f1()
!          qi => q_im(i)%get_f1()
!          do j = 0, n1p+1
!             r = j + noff -1
!             qr(1,j) = qr(1,j) / r
!             qi(1,j) = qi(1,j) / r
!          enddo
!       enddo

!    endif

!    call stop_tprof( 'deposit 2D particles' )
!    call write_dbg(cls_name, sname, cls_level, 'ends')

! end subroutine part2d_qdeposit
!
subroutine part2d_amjdeposit(x,p,gamma,q,psi,npp,dr,dt,qbm,ef_re,ef_im,bf_re,bf_im,&
&cu_re,cu_im,dcu_re,dcu_im,amu_re,amu_im,num_modes)

   implicit none

   ! real, dimension(:,:), pointer, intent(in) :: part
   real, dimension(:,:), pointer, intent(in) :: x, p
   real, dimension(:), pointer, intent(in) :: gamma, q, psi
   integer(kind=LG), intent(in) :: npp
   real, intent(in) :: dr, dt, qbm
   class(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im, &
   &bf_re, bf_im, cu_re, cu_im, dcu_re, dcu_im, amu_re, amu_im
   integer, intent(in) :: num_modes
! local data
   character(len=20), save :: sname = "part2d_amjdeposit"
   integer :: i, j
   integer(kind=LG) :: ii
   real, dimension(:,:), pointer :: e0,b0,cu0,dcu0,amu0
   real, dimension(:,:), pointer :: er,br,cur,dcur,amur
   real, dimension(:,:), pointer :: ei,bi,cui,dcui,amui
   integer :: nn, noff, n1p
   real :: qtmh, qtmh1, qtmh2, dti, p6, p7, ip7
   real :: r0, r, qc, qc1, dd, ad, rcr, rci
   real, dimension(3) :: dx, dxx, ox, oxx
   real :: ddx, ddy, vx, vy, acx, acy, omzt, omt, anorm, rot1, rot2
   real :: v1, v2, v3, pcos, psin
   complex(kind=DB) :: rc, rc0

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'deposit 2D particles' )

   noff = ef_re(0)%get_noff(1)
   n1p = ef_re(0)%get_ndp(1)
   e0 => ef_re(0)%get_f1()
   b0 => bf_re(0)%get_f1()
   cu0 => cu_re(0)%get_f1()
   dcu0 => dcu_re(0)%get_f1()
   amu0 => amu_re(0)%get_f1()
   er => null(); ei => null(); br => null(); bi => null()
   cur => null(); cui => null(); dcur => null(); dcui => null()
   amur => null(); amui => null()

   qtmh = 0.5*qbm*dt
   dti = 1.0/dt
   do ii = 1, npp
      ! r0 = sqrt( part(1,ii)**2 + part(2,ii)**2 )
      ! pcos = part(1,ii) / r0
      ! psin = part(2,ii) / r0
      ! p6 = part(6,ii)
      ! p7 = part(7,ii)
      ! qc = part(8,ii)
      r0 = sqrt( x(1,ii)**2 + x(2,ii)**2 )
      pcos = x(1,ii) / r0
      psin = x(2,ii) / r0
      p6 = gamma(ii)
      p7 = psi(ii)
      qc = q(ii)
      rc0 = cmplx( pcos, psin, kind=DB )
      r = r0/dr + 1.0
      nn = r
      dd = r - real(nn)
      ad = 1.0 - dd
      nn = nn - noff
      ip7 = 1.0 / p7
      qtmh1 = qtmh * ip7
      qtmh2 = qtmh1 * p6
! find E and B field
      dx(1:3) = ad*e0(1:3,nn)
      dx(1:3) = dd*e0(1:3,nn+1) + dx(1:3)
      ox(1:3) = ad*b0(1:3,nn)
      ox(1:3) = dd*b0(1:3,nn+1) + ox(1:3)
      rc = rc0
      do i = 1, num_modes
         rcr = 2.0*real(rc)
         rci = 2.0*aimag(rc)
         er => ef_re(i)%get_f1()
         ei => ef_im(i)%get_f1()
         dxx(1:3) = (ad*er(1:3,nn)+dd*er(1:3,nn+1))*rcr
         dxx(1:3) = dxx(1:3) - (ad*ei(1:3,nn)+dd*ei(1:3,nn+1))*rci
         dx(1:3) = dx(1:3) + dxx(1:3)
         br => bf_re(i)%get_f1()
         bi => bf_im(i)%get_f1()
         oxx(1:3) = (ad*br(1:3,nn)+dd*br(1:3,nn+1))*rcr
         oxx(1:3) = oxx(1:3) - (ad*bi(1:3,nn)+dd*bi(1:3,nn+1))*rci
         ox(1:3) = ox(1:3) + oxx(1:3)
         rc = rc*rc0
      end do
      dx(1) = -dx(1) + ox(2)
      dx(2) = -dx(2) - ox(1)
! calculate half impulse
      ddx = (-1.0)*qtmh2*dx(1) + qtmh*ox(2)
      ddy = (-1.0)*qtmh2*dx(2) - qtmh*ox(1)
! half acceleration
      ! vx = part(3,ii) * pcos + part(4,ii) * psin
      ! vy = part(4,ii) * pcos - part(3,ii) * psin
      vx = p(1,ii) * pcos + p(2,ii) * psin
      vy = p(2,ii) * pcos - p(1,ii) * psin
      acx = vx + ddx
      acy = vy + ddy
! find inverse gamma
! renormalize magnetic field
! calculate cyclotron frequency
! calculate rotation matrix
      omzt = qtmh1*ox(3)
      omt = omzt*omzt
      anorm = 2.0/(1.0 + omt)
      rot1 = 0.5*(1.0 - omt)
      rot2 = omzt
! new momentum
      v1 = (rot1*acx + rot2*acy)*anorm + ddx
      v2 = (rot1*acy - rot2*acx)*anorm + ddy
! deposit momentum flux, acceleration density, and current density
      qc1 = qc*ip7
      ad = qc1*ad
      dd = qc1*dd

      ox(1) = 0.5*(v1 + vx)
      ox(2) = 0.5*(v2 + vy)
      ox(3) = 0.5*(1.0+ox(1)*ox(1)+ox(2)*ox(2))*ip7-0.5*p7

      ! part(6,ii) = p7 + ox(3)
      gamma(ii) = p7 + ox(3)
      vx = (v1 - vx)*dti
      vy = (v2 - vy)*dti
      dx(3) = qbm*(dx(3)+(dx(1)*ox(1)+dx(2)*ox(2))*ip7)

      v1 = ox(1)*ox(1)*ip7
      v2 = ox(2)*ox(1)*ip7
      v3 = ox(2)*ox(2)*ip7

      vx = vx+ox(1)*dx(3)*ip7
      vy = vy+ox(2)*dx(3)*ip7

      ! dx(1) = ox(1)*ox(2)/r0*ip7
      ! dx(2) = ox(2)*ox(2)/r0*ip7

      dxx(1) = v1 * dd
      v1 = v1 * ad
      dxx(2) = v2 * dd
      v2 = v2 * ad
      rot1 = v3 * dd
      v3 = v3 * ad

      oxx(1) = vx * dd
      vx = vx * ad
      oxx(2) = vy * dd
      vy = vy * ad

      ! rot1 = dx(1) * dd
      ! dx(1) = dx(1) * ad
      ! rot2 = dx(2) * dd
      ! dx(2) = dx(2) * ad

      omzt = ox(1) * dd
      ox(1) = ox(1) * ad
      omt = ox(2) * dd
      ox(2) = ox(2) * ad
      anorm = ox(3) * dd
      ox(3) = ox(3) * ad

      amu0(1,nn) = amu0(1,nn) + v1
      amu0(2,nn) = amu0(2,nn) + v2
      amu0(3,nn) = amu0(3,nn) + v3
      amu0(1,nn+1) = amu0(1,nn+1) + dxx(1)
      amu0(2,nn+1) = amu0(2,nn+1) + dxx(2)
      amu0(3,nn+1) = amu0(3,nn+1) + rot1

      dcu0(1,nn) = dcu0(1,nn) + vx
      dcu0(2,nn) = dcu0(2,nn) + vy
      dcu0(1,nn+1) = dcu0(1,nn+1) + oxx(1)
      dcu0(2,nn+1) = dcu0(2,nn+1) + oxx(2)

      cu0(1,nn) = cu0(1,nn) + ox(1)
      cu0(2,nn) = cu0(2,nn) + ox(2)
      cu0(3,nn) = cu0(3,nn) + ox(3)
      cu0(1,nn+1) = cu0(1,nn+1) + omzt
      cu0(2,nn+1) = cu0(2,nn+1) + omt
      cu0(3,nn+1) = cu0(3,nn+1) + anorm

      rc0 = cmplx( pcos, -psin, kind=DB )
      rc = rc0
      do i = 1, num_modes
         rcr = real(rc)
         rci = aimag(rc)
         amur => amu_re(i)%get_f1()
         amui => amu_im(i)%get_f1()
         amur(1,nn) = amur(1,nn) + v1*rcr
         amur(2,nn) = amur(2,nn) + v2*rcr
         amur(3,nn) = amur(3,nn) + v3*rcr
         amur(1,nn+1) = amur(1,nn+1) + dxx(1)*rcr
         amur(2,nn+1) = amur(2,nn+1) + dxx(2)*rcr
         amur(3,nn+1) = amur(3,nn+1) + rot1*rcr
         amui(1,nn) = amui(1,nn) + v1*rci
         amui(2,nn) = amui(2,nn) + v2*rci
         amui(3,nn) = amui(3,nn) + v3*rci
         amui(1,nn+1) = amui(1,nn+1) + dxx(1)*rci
         amui(2,nn+1) = amui(2,nn+1) + dxx(2)*rci
         amui(3,nn+1) = amui(3,nn+1) + rot1*rci

         dcur => dcu_re(i)%get_f1()
         dcui => dcu_im(i)%get_f1()
         dcur(1,nn) = dcur(1,nn) + vx*rcr
         dcur(2,nn) = dcur(2,nn) + vy*rcr
         dcur(1,nn+1) = dcur(1,nn+1) + oxx(1)*rcr
         dcur(2,nn+1) = dcur(2,nn+1) + oxx(2)*rcr
         dcui(1,nn) = dcui(1,nn) + vx*rci
         dcui(2,nn) = dcui(2,nn) + vy*rci
         dcui(1,nn+1) = dcui(1,nn+1) + oxx(1)*rci
         dcui(2,nn+1) = dcui(2,nn+1) + oxx(2)*rci

         cur => cu_re(i)%get_f1()
         cui => cu_im(i)%get_f1()
         cur(1,nn) = cur(1,nn) + ox(1)*rcr
         cur(2,nn) = cur(2,nn) + ox(2)*rcr
         cur(3,nn) = cur(3,nn) + ox(3)*rcr
         cur(1,nn+1) = cur(1,nn+1) + omzt*rcr
         cur(2,nn+1) = cur(2,nn+1) + omt*rcr
         cur(3,nn+1) = cur(3,nn+1) + anorm*rcr
         cui(1,nn) = cui(1,nn) + ox(1)*rci
         cui(2,nn) = cui(2,nn) + ox(2)*rci
         cui(3,nn) = cui(3,nn) + ox(3)*rci
         cui(1,nn+1) = cui(1,nn+1) + omzt*rci
         cui(2,nn+1) = cui(2,nn+1) + omt*rci
         cui(3,nn+1) = cui(3,nn+1) + anorm*rci

         rc = rc*rc0
      end do
   end do

   if (noff == 0) then

      ! guard cells on the axis are useless
      cu0(1:3,0)  = 0.0
      dcu0(1:2,0) = 0.0
      amu0(1:3,0) = 0.0
      !
      cu0(1:2,1)  = 0.0
      cu0(3,1)    = 8.0 * cu0(3,1)
      dcu0(1:2,1) = 0.0
      amu0(1:3,1) = 0.0
      do i = 1, num_modes
         amur => amu_re(i)%get_f1()
         amui => amu_im(i)%get_f1()
         dcur => dcu_re(i)%get_f1()
         dcui => dcu_im(i)%get_f1()
         cur => cu_re(i)%get_f1()
         cui => cu_im(i)%get_f1()
         ! guard cells on the axis are useless
         cur(1:3,0)  = 0.0
         dcur(1:2,0) = 0.0
         amur(1:3,0) = 0.0
         cui(1:3,0)  = 0.0
         dcui(1:2,0) = 0.0
         amui(1:3,0) = 0.0

         if ( i == 1 ) then
            cur(1:2,1)  = 8.0 * cur(1:2,1)
            cur(3,1)    = 0.0
            dcur(1:2,1) = 8.0 * dcur(1:2,1)
            amur(1:3,1) = 0.0
            cui(1:2,1)  = 8.0 * cui(1:2,1)
            cui(3,1)    = 0.0
            dcui(1:2,1) = 8.0 * dcui(1:2,1)
            amui(1:3,1) = 0.0
         elseif ( i == 2 ) then
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
      end do

      do j = 2, n1p+1
         r = j + noff - 1
         cu0(1:3,j) = cu0(1:3,j) / r
         dcu0(1:2,j) = dcu0(1:2,j) / r
         amu0(1:3,j) = amu0(1:3,j) / r
      end do

      do i = 1, num_modes
         amur => amu_re(i)%get_f1()
         amui => amu_im(i)%get_f1()
         dcur => dcu_re(i)%get_f1()
         dcui => dcu_im(i)%get_f1()
         cur => cu_re(i)%get_f1()
         cui => cu_im(i)%get_f1()
         do j = 2, n1p+1
            r = j + noff - 1
            cur(1:3,j)  = cur(1:3,j) / r
            dcur(1:2,j) = dcur(1:2,j) / r
            amur(1:3,j) = amur(1:3,j) / r
            cui(1:3,j)  = cui(1:3,j) / r
            dcui(1:2,j) = dcui(1:2,j) / r
            amui(1:3,j) = amui(1:3,j) / r
         end do
      end do

   else

      do j = 0, n1p+1
         r = j + noff - 1
         cu0(1:3,j) = cu0(1:3,j) / r
         dcu0(1:2,j) = dcu0(1:2,j) / r
         amu0(1:3,j) = amu0(1:3,j) / r
      end do

      do i = 1, num_modes
         amur => amu_re(i)%get_f1()
         amui => amu_im(i)%get_f1()
         dcur => dcu_re(i)%get_f1()
         dcui => dcu_im(i)%get_f1()
         cur => cu_re(i)%get_f1()
         cui => cu_im(i)%get_f1()
         do j = 0, n1p+1
            r = j + noff - 1
            cur(1:3,j) = cur(1:3,j) / r
            dcur(1:2,j) = dcur(1:2,j) / r
            amur(1:3,j) = amur(1:3,j) / r
            cui(1:3,j) = cui(1:3,j) / r
            dcui(1:2,j) = dcui(1:2,j) / r
            amui(1:3,j) = amui(1:3,j) / r
         end do
      end do

   end if

   call stop_tprof( 'deposit 2D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine part2d_amjdeposit
!
subroutine part2d_push(x,p,gamma,q,psi,npp,dr,xdim,dt,qbm,ef_re,ef_im,&
&bf_re,bf_im,num_modes)

   implicit none

   ! real, dimension(:,:), pointer, intent(inout) :: part
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:), pointer, intent(inout) :: gamma, psi, q
   integer(kind=LG), intent(inout) :: npp
   real, intent(in) :: dr, dt, qbm
   class(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im, &
   &bf_re, bf_im
   integer, intent(in) :: xdim, num_modes
! local data
   character(len=20), save :: sname = "part2d_push"
   integer :: i
   integer(kind=LG) :: ii
   real, dimension(:,:), pointer :: e0,b0,er,ei,br,bi
   integer :: n1, nn, noff, n1p
   real :: edge, qtmh, qtmh1, qtmh2, p6, p7, ip7
   real :: r0, r, rn, qc, dd, ad, rcr, rci
   real, dimension(3) :: dx, dxx, ox, oxx, om
   real :: acx, acy, acz, omt, anorm
   real :: rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
   real :: dtc1, pos_x, pos_y, pcos, psin
   complex(kind=DB) :: rc, rc0

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'push 2D particles' )

   noff = ef_re(0)%get_noff(1)
   n1p = ef_re(0)%get_ndp(1)
   n1 = ef_re(0)%get_nd(1)
   e0 => ef_re(0)%get_f1()
   b0 => bf_re(0)%get_f1()
   er => null(); ei => null(); br => null(); bi => null()

   qtmh = 0.5*qbm*dt
   edge = real(n1) * dr

   ii = 1
   do
      if (ii > npp) exit
      ! r0 = sqrt( part(1,ii)**2 + part(2,ii)**2 )
      ! ! th = part(2,ii)
      ! pcos = part(1,ii) / r0
      ! psin = part(2,ii) / r0
      ! p6 = part(6,ii)
      ! p7 = part(7,ii)
      ! qc = part(8,ii)
      r0 = sqrt( x(1,ii)**2 + x(2,ii)**2 )
      pcos = x(1,ii) / r0
      psin = x(2,ii) / r0
      p6 = gamma(ii)
      p7 = psi(ii)
      qc = q(ii)
      rc0 = cmplx(pcos, psin, kind=DB)
      r = r0/dr + 1.0
      nn = r
      dd = r - real(nn)
      ad = 1.0 - dd
      nn = nn - noff
      ip7 = 1.0 / p7
      qtmh1 = qtmh * ip7
      qtmh2 = qtmh1 * p6

      dx(1:3) = ad*e0(1:3,nn)
      dx(1:3) = dd*e0(1:3,nn+1) + dx(1:3)
      ox(1:3) = ad*b0(1:3,nn)
      ox(1:3) = dd*b0(1:3,nn+1) + ox(1:3)
      rc = rc0
      do i = 1, num_modes
         rcr = 2.0*real(rc)
         rci = 2.0*aimag(rc)
         er => ef_re(i)%get_f1()
         ei => ef_im(i)%get_f1()
         dxx(1:3) = (ad*er(1:3,nn)+dd*er(1:3,nn+1))*rcr
         dxx(1:3) = dxx(1:3) - (ad*ei(1:3,nn)+dd*ei(1:3,nn+1))*rci
         dx(1:3) = dx(1:3) + dxx(1:3)
         br => bf_re(i)%get_f1()
         bi => bf_im(i)%get_f1()
         oxx(1:3) = (ad*br(1:3,nn)+dd*br(1:3,nn+1))*rcr
         oxx(1:3) = oxx(1:3) - (ad*bi(1:3,nn)+dd*bi(1:3,nn+1))*rci
         ox(1:3) = ox(1:3) + oxx(1:3)
         rc = rc*rc0
      end do

      ! transform (Er, Eth) to (Ex, Ey)
      ! dxx, oxx are for temporary use here
      dxx(1) = dx(1) * pcos - dx(2) * psin
      dxx(2) = dx(1) * psin + dx(2) * pcos
      dx(1:2) = dxx(1:2)
      ! transform (Br, Bth) to (Bx, By)
      oxx(1) = ox(1) * pcos - ox(2) * psin
      oxx(2) = ox(1) * psin + ox(2) * pcos
      ox(1:2) = oxx(1:2)

! calculate half impulse
      dx(1:3) = qtmh2*dx(1:3)
! half acceleration
      acx = p(1,ii) + dx(1)
      acy = p(2,ii) + dx(2)
      acz = p(3,ii) + dx(3)

      om(1:3) = qtmh1*ox(1:3)
      omt = om(1)*om(1) + om(2)*om(2) + om(3)*om(3)
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = om(1)*om(2)
      rot7 = om(1)*om(3)
      rot8 = om(2)*om(3)
      rot1 = omt + om(1)*om(1)
      rot5 = omt + om(2)*om(2)
      rot9 = omt + om(3)*om(3)
      rot2 = om(3) + rot4
      rot4 = -om(3) + rot4
      rot3 = -om(2) + rot7
      rot7 = om(2) + rot7
      rot6 = om(1) + rot8
      rot8 = -om(1) + rot8
! new momentum
      dx(1) = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx(1)
      dx(2) = (rot4*acx + rot5*acy + rot6*acz)*anorm + dx(2)
      dx(3) = (rot7*acx + rot8*acy + rot9*acz)*anorm + dx(3)

      dtc1 = dt/(sqrt(1+dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))-dx(3))

      ! new position
      pos_x = x(1,ii) + dx(1) * dtc1
      pos_y = x(2,ii) + dx(2) * dtc1
      rn = sqrt(pos_x**2 + pos_y**2)

      if (rn >= edge) then
         if (ii == npp) then
            npp = npp -1
            exit
         else
            do i = 1, xdim
               ! part(i,ii) = part(i,npp)
               x(:,ii)   = x(:,npp)
               p(:,ii)   = p(:,npp)
               gamma(ii) = gamma(npp)
               q(ii)     = q(npp)
               psi(ii)   = psi(npp)
            end do
            npp = npp - 1
            cycle
         end if
      else
         x(1,ii) = pos_x
         x(2,ii) = pos_y
         p(1,ii) = dx(1)
         p(2,ii) = dx(2)
         p(3,ii) = dx(3)
         ii = ii + 1
      end if
   end do

   call stop_tprof( 'push 2D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine part2d_push
!
subroutine part2d_pmove(x,p,gamma,q,psi,pp,npp,dr,xdim,npmax,nbmax,ud,sbufl,sbufr,rbufl,rbufr,ihole)

   implicit none

   ! real, dimension(:,:), pointer, intent(inout) :: part
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:), pointer, intent(inout) :: gamma, q, psi
   real, dimension(:,:), intent(inout) :: sbufl,sbufr,rbufl,rbufr
   integer(kind=LG), intent(inout) :: npp
   real, intent(in) :: dr
   integer(kind=LG), intent(in) :: npmax, nbmax
   integer(kind=LG), dimension(:), intent(inout) :: ihole
   integer, intent(in) :: xdim
   class(parallel_pipe), pointer, intent(in) :: pp
   class(ufield), intent(in) :: ud
! local data
   character(len=20), save :: sname = "part2d_pmove"
   integer :: i, j1, j2
   integer(kind=LG) :: j, npt
   integer, dimension(2) :: jsr, jsl, jss
   integer, dimension(4) :: ibflg, msid, iwork
   integer, dimension(MPI_STATUS_SIZE) :: istatus
   integer, dimension(10) :: info
   integer :: itermax, nvp, n1, n1p, noff, iter, nter, nps
   integer :: id, idl, idr, mter, nbsize
   integer :: mreal, mint, lworld, lgrp, ierr
   real :: edgel,edger,an,xt

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'move 2D particles' )

   nbsize = nbmax*xdim
   itermax = 20000
   ierr = 0
! buffer outgoing particles, first in r direction
   nvp = pp%getlnvp()
   mreal = pp%getmreal()
   mint = pp%getmint()
   lworld = pp%getlgrp()
   lgrp = pp%getlgrp()
   id = pp%getlidproc()
   idl = id - 1
   if (idl < 0) idl = idl + nvp
   idr = id + 1
   if (idr >= nvp) idr = idr - nvp
   do j = 1, 10
      info(j) = 0
   end do
   n1 = ud%get_nd(1)
   noff = ud%get_noff(1)
   n1p = ud%get_ndp(1)
   an = real(n1) - 0.5
   if (id == 0) then
      edgel = 0.0
      edger = real(n1p) * dr
   else
      edgel = real(noff) * dr
      edger = edgel + real(n1p)*dr
   end if
   iter = 2
   nter = 0
   do
      mter = 0
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do j = 1, npp
         ! xt = sqrt( part(1,j)**2 + part(2,j)**2 )
         xt = sqrt( x(1,j)**2 + x(2,j)**2 )
! particles going down
         if (xt < edgel) then
            if (jsl(1) < nbmax) then
               jsl(1) = jsl(1) + 1
               ! do i = 1, xdim
               !    sbufl(i,jsl(1)) = part(i,j)
               ! end do
               sbufl(1:2,jsl(1)) = x(:,j)
               sbufl(3:5,jsl(1)) = p(:,j)
               sbufl(6,jsl(1)) = gamma(j)
               sbufl(7,jsl(1)) = psi(j)
               sbufl(8,jsl(1)) = q(j)

               ihole(jsl(1)+jsr(1)) = j
            else
               jss(2) = 1
               exit
            end if
! particles going up
         else if (xt >= edger) then
            if (jsr(1) < nbmax) then
               jsr(1) = jsr(1) + 1
               ! do i = 1, xdim
               !    sbufr(i,jsr(1)) = part(i,j)
               ! end do
               sbufr(1:2,jsr(1)) = x(:,j)
               sbufr(3:5,jsr(1)) = p(:,j)
               sbufr(6,jsr(1)) = gamma(j)
               sbufr(7,jsr(1)) = psi(j)
               sbufr(8,jsr(1)) = q(j)

               ihole(jsl(1)+jsr(1)) = j
            else
               jss(2) = 1
               exit
            end if
         end if
      end do
      jss(1) = jsl(1) + jsr(1)
! check for full buffer condition
      nps = 0
      nps = max0(nps,jss(2))
      ibflg(3) = nps
! copy particle buffers
      do
         iter = iter + 2
         mter = mter + 1
! post receive
         call MPI_IRECV(rbufl,nbsize,mreal,idl,iter-1,lgrp,msid(1),ierr)
         call MPI_IRECV(rbufr,nbsize,mreal,idr,iter,lgrp,msid(2),ierr)
! send particles
         call MPI_ISEND(sbufr,xdim*jsr(1),mreal,idr,iter-1,lgrp,msid(3),ierr)
         call MPI_ISEND(sbufl,xdim*jsl(1),mreal,idl,iter,lgrp,msid(4),ierr)
! wait for particles to arrive
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsl(2) = nps/xdim
         call MPI_WAIT(msid(2),istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsr(2) = nps/xdim
! check if particles must be passed further
         nps = 0
! check if any particles coming from above belong here
         jsl(1) = 0
         jsr(1) = 0
         jss(2) = 0
         do j = 1, jsr(2)
            xt = sqrt(rbufr(1,j)**2+rbufr(2,j)**2)
            if (xt < edgel) jsl(1) = jsl(1) + 1
            if (xt >= edger) jsr(1) = jsr(1) + 1
         end do
         if (jsr(1) /= 0) then
            write (erstr,*) 'Info:',jsr(1),' particles returning above'
            call write_dbg(cls_name, sname, cls_level, erstr)
         end if
! check if any particles coming from below belong here
         do j = 1, jsl(2)
            xt = sqrt(rbufl(1,j)**2+rbufl(2,j)**2)
            if (xt >= edger) jsr(1) = jsr(1) + 1
            if (xt < edgel) jss(2) = jss(2) + 1
         end do
         if (jss(2) /= 0) then
            write (erstr,*) 'Info:',jss(2),' particles returning below'
            call write_dbg(cls_name, sname, cls_level, erstr)
         end if
         jsl(1) = jsl(1) + jss(2)
         nps = max0(nps,jsl(1)+jsr(1))

         ibflg(2) = nps
! make sure sbufr and sbufl have been sent
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
         if (nps /= 0) then
! remove particles which do not belong here
! first check particles coming from above
            jsl(1) = 0
            jsr(1) = 0
            jss(2) = 0
            do j = 1, jsr(2)
               xt = sqrt(rbufr(1,j)**2+rbufr(2,j)**2)
! particles going down
               if (xt < edgel) then
                  jsl(1) = jsl(1) + 1
                  do i = 1, xdim
                     sbufl(i,jsl(1)) = rbufr(i,j)
                  end do
! particles going up, should not happen
               else if (xt >= edger) then
                  jsr(1) = jsr(1) + 1
                  do i = 1, xdim
                     sbufr(i,jsr(1)) = rbufr(i,j)
                  end do
! particles staying here
               else
                  jss(2) = jss(2) + 1
                  do i = 1, xdim
                     rbufr(i,jss(2)) = rbufr(i,j)
                  end do
               end if
            end do
            jsr(2) = jss(2)
! next check particles coming from below
            jss(2) = 0
            do j = 1, jsl(2)
               xt = sqrt(rbufl(1,j)**2+rbufl(2,j)**2)
! particles going up
               if (xt >= edger) then
                  if (jsr(1) < nbmax) then
                     jsr(1) = jsr(1) + 1
                     do i = 1, xdim
                        sbufr(i,jsr(1)) = rbufl(i,j)
                     end do
                  else
                     jss(2) = 2*npmax
                     exit
                  end if
! particles going down back, should not happen
               else if (xt < edgel) then
                  if (jsl(1) < nbmax) then
                     jsl(1) = jsl(1) + 1
                     do i = 1, xdim
                        sbufl(i,jsl(1)) = rbufl(i,j)
                     end do
                  else
                     jss(2) = 2*npmax
                     exit
                  end if
! particles staying here
               else
                  jss(2) = jss(2) + 1
                  do i = 1, xdim
                     rbufl(i,jss(2)) = rbufl(i,j)
                  end do
               end if
            end do
            jsl(2) = jss(2)
         end if
! check if move would overflow particle array
         nps = 0
         npt = npmax
         jss(2) = npp + jsl(2) + jsr(2) - jss(1)
         nps = max0(nps,jss(2))
         npt = min0(npt,jss(2))
         ibflg(1) = nps
         ibflg(4) = -npt
         iwork = ibflg
         call MPI_ALLREDUCE(iwork,ibflg,4,mint,MPI_MAX,lgrp,ierr)
         info(2) = ibflg(1)
         info(3) = -ibflg(4)
         ierr = ibflg(1) - npmax
         if (ierr > 0) then
            write (erstr,*) 'particle overflow error, ierr = ', ierr
            call write_dbg(cls_name, sname, cls_level, erstr)
            info(1) = ierr
            return
         end if
! distribute incoming particles from buffers
! distribute particles coming from below into holes
         jss(2) = min0(jss(1),jsl(2))
         do j = 1, jss(2)
            ! do i = 1, xdim
            !    part(i,ihole(j)) = rbufl(i,j)
            ! end do
            x(:,ihole(j))   = rbufl(1:2,j)
            p(:,ihole(j))   = rbufl(3:5,j)
            gamma(ihole(j)) = rbufl(6,j)
            psi(ihole(j))   = rbufl(7,j)
            q(ihole(j))     = rbufl(8,j)
         end do
         if (jss(1) > jsl(2)) then
            jss(2) = min0(jss(1)-jsl(2),jsr(2))
         else
            jss(2) = jsl(2) - jss(1)
         end if
         do j = 1, jss(2)
! no more particles coming from below or back
! distribute particles coming from above or front into holes
            if (jss(1) > jsl(2)) then
               ! do i = 1, xdim
               !    part(i,ihole(j+jsl(2))) = rbufr(i,j)
               ! end do
               x(:,ihole(j+jsl(2)))   = rbufr(1:2,j)
               p(:,ihole(j+jsl(2)))   = rbufr(3:5,j)
               gamma(ihole(j+jsl(2))) = rbufr(6,j)
               psi(ihole(j+jsl(2)))   = rbufr(7,j)
               q(ihole(j+jsl(2)))     = rbufr(8,j)
            else
! no more holes
! distribute remaining particles from below or back into bottom
               ! do i = 1, xdim
                  ! part(i,j+npp) = rbufl(i,j+jss(1))
               ! end do
               x(:,j+npp)   = rbufl(1:2,j+jss(1))
               p(:,j+npp)   = rbufl(3:5,j+jss(1))
               gamma(j+npp) = rbufl(6,j+jss(1))
               psi(j+npp)   = rbufl(7,j+jss(1))
               q(j+npp)     = rbufl(8,j+jss(1))
            end if
            end do
            if (jss(1) <= jsl(2)) then
               npp = npp + (jsl(2) - jss(1))
               jss(1) = jsl(2)
            end if
            jss(2) = jss(1) - (jsl(2) + jsr(2))
            if (jss(2) > 0) then
               jss(1) = (jsl(2) + jsr(2))
               jsr(2) = jss(2)
            else
               jss(1) = jss(1) - jsl(2)
               jsr(2) = -jss(2)
            end if
            do j = 1, jsr(2)
! holes left over
! fill up remaining holes in particle array with particles from bottom
            if (jss(2) > 0) then
               j1 = npp - j + 1
               j2 = jss(1) + jss(2) - j + 1
               if (j1 > ihole(j2)) then
! move particle only if it is below current hole
                  ! do i = 1, xdim
                     ! part(i,ihole(j2)) = part(i,j1)
                  ! end do
                  x(:,ihole(j2))   = x(:,j1)
                  p(:,ihole(j2))   = p(:,j1)
                  gamma(ihole(j2)) = gamma(j1)
                  psi(ihole(j2))   = psi(j1)
                  q(ihole(j2))     = q(j1)
               end if
            else
! no more holes
! distribute remaining particles from above or front into bottom
               ! do i = 1, xdim
               !    part(i,j+npp) = rbufr(i,j+jss(1))
               ! end do
               x(:,j+npp)   = rbufr(1:2,j+jss(1))
               p(:,j+npp)   = rbufr(3:5,j+jss(1))
               gamma(j+npp) = rbufr(6,j+jss(1))
               psi(j+npp)   = rbufr(7,j+jss(1))
               q(j+npp)     = rbufr(8,j+jss(1))
            end if
         end do
         if (jss(2) > 0) then
            npp = npp - jsr(2)
         else
            npp = npp + jsr(2)
         end if
         jss(1) = 0
! check if any particles have to be passed further
         info(5) = max0(info(5),mter)
         if (ibflg(2) <= 0) exit
         write (erstr,*) 'Info: particles being passed further = ', ibflg(2)
         call write_dbg(cls_name, sname, cls_level, erstr)
         if (ibflg(3) > 0) ibflg(3) = 1
         if (iter >= itermax) then
            ierr = -((iter-2)/2)
            write (erstr,*) 'Iteration overflow, iter = ', ierr
            call write_err(erstr)
            info(1) = ierr
            if (nter > 0) then
               write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
            endif
            return
         end if
      end do
! check if buffer overflowed and more particles remain to be checked
      if (ibflg(3) <=  0) exit
      nter = nter + 1
      info(4) = nter
      write (erstr,*) "new loop, nter=", nter
      call write_dbg(cls_name, sname, cls_level, erstr)
   end do
   if (nter > 0) then
      write (erstr,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      call write_dbg(cls_name, sname, cls_level, erstr)
   end if

   call stop_tprof( 'move 2D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')
   return
end subroutine part2d_pmove
!
subroutine part2d_extractpsi(x,p,gamma,q,psi,npp,dr,qbm,psi_re,psi_im,num_modes)

   implicit none

   ! real, dimension(:,:), pointer, intent(inout) :: part
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:), pointer, intent(inout) :: gamma, q, psi
   integer(kind=LG), intent(in) :: npp
   integer, intent(in) :: num_modes
   real, intent(in) :: dr,qbm
   class(ufield), dimension(:), pointer, intent(in) :: psi_re, psi_im
! local data
   character(len=20), save :: sname = "part2d_extractpsi"
   real, dimension(:,:), pointer :: psi0,psir,psii
   real :: r0, r, vx, vy, dx, ad, dd, rci, rcr
   complex(kind=DB) :: rc, rc0
   integer(kind=LG) :: ii
   integer :: i, noff, n1p, nn

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'extract psi' )

   noff = psi_re(0)%get_noff(1)
   n1p = psi_re(0)%get_ndp(1)
   psi0 => psi_re(0)%get_f1()
   psir => null(); psii => null()

   do ii = 1, npp
      ! r0 = sqrt(part(1,ii)**2+part(2,ii)**2)
      ! ! th = part(2,ii)
      ! vx = part(3,ii)
      ! vy = part(4,ii)
      ! rc0 = cmplx(part(1,ii),part(2,ii),kind=DB) / r0
      r0 = sqrt(x(1,ii)**2+x(2,ii)**2)
      vx = p(1,ii)
      vy = p(2,ii)
      rc0 = cmplx(x(1,ii),x(2,ii),kind=DB) / r0
      r = r0/dr + 1.0
      nn = r
      dd = r - real(nn)
      ad = 1.0 - dd
      nn = nn - noff
      dx = ad*psi0(1,nn)
      dx = dd*psi0(1,nn+1) + dx
      rc = rc0
      do i = 1, num_modes
         rcr = 2.0*real(rc)
         rci = 2.0*aimag(rc)
         psir => psi_re(i)%get_f1()
         psii => psi_im(i)%get_f1()
         dx = dx + (ad*psir(1,nn)+dd*psir(1,nn+1))*rcr
         dx = dx - (ad*psii(1,nn)+dd*psii(1,nn+1))*rci
         rc = rc*rc0
      end do
      dx = - dx*qbm
      ! part(7,ii) = 1.0 + dx
      ! part(6,ii) = (vx**2+vy**2+1.0)/(2.0*(1.0+dx))+0.5*(1.0+dx)
      psi(ii) = 1.0 + dx
      gamma(ii) = (vx**2+vy**2+1.0)/(2.0*(1.0+dx))+0.5*(1.0+dx)
   end do

   call stop_tprof( 'extract psi' )
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine part2d_extractpsi
!
end module part2d_lib
