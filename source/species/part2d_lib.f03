! part2d_lib_h module for QuickPIC Open Source 1.0
! update: 04/18/2016

module part2d_lib

use parallel_pipe_class
use ufield_class
use param
use system
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
subroutine part2d_qdeposit(part,npp,q_re,q_im,num_modes)

   implicit none

   real, dimension(:,:), pointer, intent(in) :: part
   integer(kind=LG), intent(in) :: npp
   class(ufield), dimension(:), pointer, intent(in) :: q_re, q_im
   integer, intent(in) :: num_modes
! local data
   character(len=20), save :: sname = "part2d_qdeposit"
   integer :: i, nn, noff, n1p
   integer(kind=LG) :: ii
   real, dimension(:,:), pointer :: q0, qr, qi
   real :: r, qc, th, dd, ad, rcr, rci
   complex(kind=DB) :: rc

   call write_dbg(cls_name, sname, cls_level, 'starts')

   noff = q_re(0)%get_noff(1)
   n1p = q_re(0)%get_ndp(1)
   q0 => q_re(0)%get_f1()
   qr => null(); qi => null

   do ii = 1, npp
      r = part(1,ii)
      th = part(2,ii)
      qc = part(8,ii)
      rc = cmplx(r*cos(th),-r*sin(th),kind=DB)
      r = r + 0.5
      nn = r - noff
      dd = qc*(r - real(nn))
      ad = qc - dd
      r = q0(1,nn) + ad
      qc = q0(1,nn+1) + dd
      q0(1,nn) = r
      q0(1,nn+1) = qc
      do i = 1, num_modes
         rcr = real(rc)
         rci = aimag(rc)
         qr => q_re(i)%get_f1()
         qi => q_im(i)%get_f1()
         r = qr(1,nn) + ad*rcr
         qc = qr(1,nn+1) + dd*rcr
         qr(1,nn) = r
         qr(1,nn+1) = qc
         r = qi(1,nn) + ad*rci
         qc = qi(1,nn+1) + dd*rci
         qi(1,nn) = r
         qi(1,nn+1) = qc
         rc = rc*rc
      end do
   end do

   if (noff == 0) then
      q0(1,0) = q0(1,0)/0.5
      do i = 1, num_modes
         qr => q_re(i)%get_f1()
         qi => q_im(i)%get_f1()
         qr(1,0) = qr(1,0)/0.5
         qi(1,0) = qi(1,0)/0.5
      end do      
   else
      q0(1,0) = q0(1,0)/(0.5+noff-1)
      do i = 1, num_modes
         qr => q_re(i)%get_f1()
         qi => q_im(i)%get_f1()
         qr(1,0) = qr(1,0)/(0.5+noff-1)
         qi(1,0) = qi(1,0)/(0.5+noff-1)
      end do      
   end if

   do i = 1, n1p
      r = 0.5 + i + noff - 1
      q0(1,i) = q0(1,i)/r
   end do
      
   do i = 1, num_modes
      qr => q_re(i)%get_f1()
      qi => q_im(i)%get_f1()
      do j = 1, num_modes
         r = 0.5 + j + noff - 1
         qr(1,j) = qr(1,j)/r
         qi(1,j) = qi(1,j)/r
      end do      
   end do
   call write_dbg(cls_name, sname, cls_level, 'ends')
end subroutine part2d_qdeposit
!
subroutine part2d_amjdeposit(part,npp,dt,qbm,ef_re,ef_im,bf_re,bf_im,&
&cu_re,cu_im,dcu_re,dcu_im,amu_re,amu_im,num_modes)

   implicit none

   real, dimension(:,:), pointer, intent(in) :: part
   integer(kind=LG), intent(in) :: npp
   real, intent(in) :: dt, qbm
   class(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im, &
   &bf_re, bf_im, cu_re, cu_im, dcu_re, dcu_im, amu_re, amu_im
   integer, intent(in) :: num_modes
! local data
   character(len=20), save :: sname = "part2d_amjdeposit"
   integer :: i, j, k
   integer(kind=LG) :: ii, jj, kk
   real, dimension(:,:), pointer :: e0,b0,cu0,dcu0,amu0
   real, dimension(:,:), pointer :: er,br,cur,dcur,amur
   real, dimension(:,:), pointer :: ei,bi,cui,dcui,amui
   integer :: nn, noff, n1p
   real :: qtmh, qtmh1, qtmh2, dti, p6, p7, ip7
   real :: r0, r, qc, qc1, th, dd, ad, rcr, rci
   real, dimension(3) :: dx, dxx, ox, oxx
   real :: ddx, ddy, vx, vy, acx, acy, omzt, omt, anorm, rot1, rot2
   real :: v1, v2, v3, v1d, v2d
   complex(kind=DB) :: rc, rc0

   call write_dbg(cls_name, sname, cls_level, 'starts')

   noff = ef_re(0)%get_noff(1)
   n1p = ef_re(0)%get_ndp(1)
   e0 => ef_re(0)%get_f1()
   b0 => bf_re(0)%get_f1()
   cu0 => cu_re(0)%get_f1()
   dcu0 => dcu_re(0)%get_f1()
   amu0 => amu_re(0)%get_f1()
   er => null(); ei => null; br => null(); bi => null
   cur => null(); cui => null; dcur => null(); dcui => null
   amur => null(); amui => null

   qtmh = 0.5*qbm*dt
   dti = 1.0/dt
   do ii = 1, npp
      r0 = part(1,ii)
      th = part(2,ii)
      p6 = part(6,ii)
      p7 = part(7,ii)
      qc = part(8,ii)
      rc0 = cmplx(r0*cos(th),-r0*sin(th),kind=DB)
      r = r0 + 0.5      
      nn = r - noff
      dd = r - real(nn)
      ad = 1.0 - dd
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
         rcr = real(rc)
         rci = aimag(rc)
         er => ef_re(i)%get_f1()
         ei => ef_im(i)%get_f1()
         dxx(1:3) = (ad*er(1:3,nn)+dd*er(1:3,nn+1))*rcr
         dxx(1:3) = dxx(1:3) + (ad*ei(1:3,nn)+dd*ei(1:3,nn+1))*rci
         dx(1:3) = dx(1:3) + dxx(1:3)
         br => bf_re(i)%get_f1()
         bi => bf_im(i)%get_f1()
         oxx(1:3) = (ad*br(1:3,nn)+dd*br(1:3,nn+1))*rcr
         oxx(1:3) = oxx(1:3) + (ad*bi(1:3,nn)+dd*bi(1:3,nn+1))*rci
         ox(1:3) = ox(1:3) + oxx(1:3)
         rc = rc*rc
      end do
! calculate half impulse
      ddx = (-1.0)*qtmh2*dx(1) + qtmh*ox(2)
      ddy = (-1.0)*qtmh2*dx(2) - qtmh*ox(1)
! half acceleration
      vx = part(3,ii)
      vy = part(4,ii)
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
      ox(3) = 0.5*(1.0+ox(1)*ox(1)+oy(1)*oy(1))*ip7-0.5*p7

      part(6,ii) = p7 + ox(3)
      vx = (v1 - vx)*dti
      vy = (v2 - vy)*dti
      dx(3) = qbm*(dx(3)+(dx(1)*ox(1)+dx(2)*ox(2))*ip7)

      v1 = ox(1)*ox(1)*ip7
      v2 = ox(2)*ox(1)*ip7
      ! v3 = ox(1)*ox(2)*ip7

      vx = vx+ox(1)*dx(3)*ip7
      vy = vy+ox(2)*dx(3)*ip7

      dx(1) = ox(1)*ox(2)/r0*ip7
      dx(2) = ox(2)*ox(2)/r0*ip7

      dxx(1) = v1 * dd
      v1 = v1 * ad
      dxx(2) = v2 * dd
      v2 = v2 * ad

      oxx(1) = vx * dd
      vx = vx * ad
      oxx(2) = vy * dd
      vy = vy * ad

      rot1 = dx(1) * dd
      dx(1) = dx(1) * ad
      rot2 = dx(2) * dd
      dx(2) = dx(2) * ad

      omzt = ox(1) * dd
      ox(1) = ox(1) * ad
      omt = ox(2) * dd
      ox(2) = ox(2) * ad
      anorm = ox(3) * dd
      ox(3) = ox(3) * ad

      amu0(1,nn) = amu0(1,nn) + v1
      amu0(2,nn) = amu0(2,nn) + v2
      amu0(1,nn+1) = amu0(1,nn+1) + dxx(1)
      amu0(2,nn+1) = amu0(2,nn+1) + dxx(2)

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

      rc = rc0
      do i = 1, num_modes
         rcr = real(rc)
         rci = aimag(rc)
         amur => amu_re(i)%get_f1()
         amui => amu_im(i)%get_f1()
         amur(1,nn) = amur(1,nn) + v1*rcr
         amur(2,nn) = amur(2,nn) + v2*rcr
         amur(1,nn+1) = amur(1,nn+1) + dxx(1)*rcr
         amur(2,nn+1) = amur(2,nn+1) + dxx(2)*rcr
         amui(1,nn) = amui(1,nn) + v1*rci
         amui(2,nn) = amui(2,nn) + v2*rci
         amui(1,nn+1) = amui(1,nn+1) + dxx(1)*rci
         amui(2,nn+1) = amui(2,nn+1) + dxx(2)*rci

         dcur => dcu_re(i)%get_f1()
         dcui => dcu_im(i)%get_f1()
         dcur(1,nn) = dcur(1,nn) + vx*rcr + dx(1)*rci*real(i)
         dcur(2,nn) = dcur(2,nn) + vy*rcr + dx(2)*rci*real(i)
         dcur(1,nn+1) = dcur(1,nn+1) + oxx(1)*rcr + rot1*rci*real(i)
         dcur(2,nn+1) = dcur(2,nn+1) + oxx(2)*rcr + rot2*rci*real(i)
         dcui(1,nn) = dcui(1,nn) + vx*rci - dx(1)*rcr*real(i)
         dcui(2,nn) = dcui(2,nn) + vy*rci - dx(2)*rcr*real(i)
         dcui(1,nn+1) = dcui(1,nn+1) + oxx(1)*rci - rot1*rcr*real(i)
         dcui(2,nn+1) = dcui(2,nn+1) + oxx(2)*rci - rot2*rcr*real(i)

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

         rc = rc*rc
      end do
   end do

   if (noff == 0) then
      cu0(1:3,0) = cu0(1:3,0)/0.5
      dcu0(1:2,0) = dcu0(1:2,0)/0.5
      amu0(1:2,0) = amu0(1:2,0)/0.5
      do i = 1, num_modes
         amur => amu_re(i)%get_f1()
         amui => amu_im(i)%get_f1()
         dcur => dcu_re(i)%get_f1()
         dcui => dcu_im(i)%get_f1()
         cur => cu_re(i)%get_f1()
         cui => cu_im(i)%get_f1()

         cur(1:3,0) = cur(1:3,0)/0.5
         dcur(1:2,0) = dcur(1:2,0)/0.5
         amur(1:2,0) = amur(1:2,0)/0.5
         cui(1:3,0) = cui(1:3,0)/0.5
         dcui(1:2,0) = dcui(1:2,0)/0.5
         amui(1:2,0) = amui(1:2,0)/0.5
      end do      
   else
      cu0(1:3,0) = cu0(1:3,0)/(0.5+noff-1)
      dcu0(1:2,0) = dcu0(1:2,0)/(0.5+noff-1)
      amu0(1:2,0) = amu0(1:2,0)/(0.5+noff-1)
      do i = 1, num_modes
         amur => amu_re(i)%get_f1()
         amui => amu_im(i)%get_f1()
         dcur => dcu_re(i)%get_f1()
         dcui => dcu_im(i)%get_f1()
         cur => cu_re(i)%get_f1()
         cui => cu_im(i)%get_f1()

         cur(1:3,0) = cur(1:3,0)/(0.5+noff-1)
         dcur(1:2,0) = dcur(1:2,0)/(0.5+noff-1)
         amur(1:2,0) = amur(1:2,0)/(0.5+noff-1)
         cui(1:3,0) = cui(1:3,0)/(0.5+noff-1)
         dcui(1:2,0) = dcui(1:2,0)/(0.5+noff-1)
         amui(1:2,0) = amui(1:2,0)/(0.5+noff-1)
      end do      
   end if

   do j = 1, n1p
      r = 0.5 + j + noff - 1
      cu0(1:3,j) = cu0(1:3,j)/r
      dcu0(1:2,j) = dcu0(1:2,j)/r
      amu0(1:2,j) = amu0(1:2,j)/r
   end do

   do i = 1, num_modes
      amur => amu_re(i)%get_f1()
      amui => amu_im(i)%get_f1()
      dcur => dcu_re(i)%get_f1()
      dcui => dcu_im(i)%get_f1()
      cur => cu_re(i)%get_f1()
      cui => cu_im(i)%get_f1()
      do j = 1, n1p
         r = 0.5 + j + noff - 1
         cur(1:3,j) = cur(1:3,j)/r
         dcur(1:2,j) = dcur(1:2,j)/r
         amur(1:2,j) = amur(1:2,j)/r
         cui(1:3,j) = cui(1:3,j)/r
         dcui(1:2,j) = dcui(1:2,j)/r
         amui(1:2,j) = amui(1:2,j)/r
      end do
   end do      
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine part2d_amjdeposit
!
subroutine part2d_push(part,npp,dt,qbm,dex,ef_re,ef_im,&
&bf_re,bf_im,num_modes)

   implicit none

   real, dimension(:,:), pointer, intent(inout) :: part
   integer(kind=LG), intent(inout) :: npp
   real, intent(in) :: dt, qbm, dex
   class(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im, &
   &bf_re, bf_im
   integer, intent(in) :: num_modes
! local data
   character(len=20), save :: sname = "part2d_push"
   integer :: i
   integer(kind=LG) :: ii
   real, dimension(:,:), pointer :: e0,b0,er,ei,br,bi
   integer :: n1, nn, noff, n1p
   real :: idex, edge, qtmh, qtmh1, qtmh2, dti, p6, p7, ip7
   real :: r0, r, rn, qc, qc1, th, th1 dd, ad, rcr, rci
   real, dimension(3) :: dx, dxx, ox, oxx, om
   real :: ddx, ddy, vx, vy, acx, acy, acz, omt, anorm
   real :: rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8
   real :: dtc1, v1, v2
   complex(kind=DB) :: rc, rc0

   call write_dbg(cls_name, sname, cls_level, 'starts')

   noff = ef_re(0)%get_noff(1)
   n1p = ef_re(0)%get_ndp(1)
   n1 = ef_re(0)%get_nd(1)
   e0 => ef_re(0)%get_f1()
   b0 => bf_re(0)%get_f1()
   er => null(); ei => null; br => null(); bi => null

   idex = 1.0/dex
   qtmh = 0.5*qbm*dt
   edge = real(n1) - 0.5

   do ii = 1, npp
      r0 = part(1,ii)
      th = part(2,ii)
      p6 = part(6,ii)
      p7 = part(7,ii)
      qc = part(8,ii)
      rc0 = cmplx(r0*cos(th),-r0*sin(th),kind=DB)
      r = r0 + 0.5      
      nn = r - noff
      dd = r - real(nn)
      ad = 1.0 - dd
      ip7 = 1.0 / p7
      qtmh1 = qtmh * ip7
      qtmh2 = qtmh1 * p6

      dx(1:3) = ad*e0(1:3,nn)
      dx(1:3) = dd*e0(1:3,nn+1) + dx(1:3)
      ox(1:3) = ad*b0(1:3,nn)
      ox(1:3) = dd*b0(1:3,nn+1) + ox(1:3)
      rc = rc0
      do i = 1, num_modes
         rcr = real(rc)
         rci = aimag(rc)
         er => ef_re(i)%get_f1()
         ei => ef_im(i)%get_f1()
         dxx(1:3) = (ad*er(1:3,nn)+dd*er(1:3,nn+1))*rcr
         dxx(1:3) = dxx(1:3) + (ad*ei(1:3,nn)+dd*ei(1:3,nn+1))*rci
         dx(1:3) = dx(1:3) + dxx(1:3)
         br => bf_re(i)%get_f1()
         bi => bf_im(i)%get_f1()
         oxx(1:3) = (ad*br(1:3,nn)+dd*br(1:3,nn+1))*rcr
         oxx(1:3) = oxx(1:3) + (ad*bi(1:3,nn)+dd*bi(1:3,nn+1))*rci
         ox(1:3) = ox(1:3) + oxx(1:3)
         rc = rc*rc
      end do

! calculate half impulse
      dx(1:3) = qtmh2*dx(1:3)
! half acceleration
      acx = part(3,ii) + dx(1)
      acy = ppart(4,ii) + dx(2)
      acz = ppart(5,ii) + dx(3)

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
      dtc1 = dt/dex/(sqrt(1+dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))-dx(3))
      v1 = r0 + dx(1)*dtc1
      v2 = dx(2)*dtc1
      rn = sqrt(v1*v1+v2*v2)
      th1 = th + atan(v2/v1)
      part(1,ii) = rn
      part(2,ii) = th1
      part(3,ii) = (dx(1)*v1+dx(2)*v2)/rn
      part(4,ii) = r0*dx(2)/rn
      part(5,ii) = dx(3)
   end do

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine part2d_push
!
subroutine part2d_pmove(part,pp,npp,ud)

   implicit none

   real, dimension(:,:), pointer, intent(inout) :: part
   integer(kind=LG), intent(inout) :: npp
   class(parallel_pipe), pointer, intent(in) :: pp
   class(ufield), intent(in) :: ud
! local data
   character(len=20), save :: sname = "part2d_pmove"
   integer :: i
   integer(kind=LG) :: ii
   real, dimension(:,:), pointer :: e0,b0,er,ei,br,bi
   integer :: n1, nn, noff, n1p

   

end subroutine part2d_pmove
!
subroutine part2d_extractpsi(part,npp,qbm,dex,psi_re,psi_im)

   implicit none

   real, dimension(:,:), pointer, intent(inout) :: part
   integer(kind=LG), intent(in) :: npp
   class(ufield), dimension(:), pointer, intent(in) :: psi_re, psi_im
! local data
   character(len=20), save :: sname = "part2d_extractpsi"
   real, dimension(:,:), pointer :: psi0,psir,psii
   real :: dx2, r0, r, th, vx, vy, dx, dxx
   complex(kind=DB) :: rc, rc0
   integer(kind=LG) :: ii
   integer :: i, noff, n1p

   call write_dbg(cls_name, sname, cls_level, 'starts')

   noff = psi_re(0)%get_noff(1)
   n1p = psi_re(0)%get_ndp(1)
   psi0 => psi_re(0)%get_f1()
   psir => null(); psii => null

   dx2 = dex * dex
   do ii = 1, npp
      r0 = part(1,ii)
      th = part(2,ii)
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      rc0 = cmplx(r0*cos(th),-r0*sin(th),kind=DB)
      r = r0 + 0.5      
      nn = r - noff
      dd = r - real(nn)
      ad = 1.0 - dd
      dx = ad*psi0(1,nn)
      dx = dd*psi0(1,nn+1) + dx
      rc = rc0
      do i = 1, num_modes
         rcr = real(rc)
         rci = aimag(rc)
         psir => psi_re(i)%get_f1()
         psii => psi_im(i)%get_f1()
         dx = (ad*psir(1,nn)+dd*psir(1,nn+1))*rcr + dx
         dx = (ad*psii(1,nn)+dd*psii(1,nn+1))*rci + dx
         rc = rc*rc
      end do
      dx = - dx*qbm
      ppart(7,ii) = 1.0 + dx
      ppart(6,ii) = (vx**2+vy**2+1.0)/(2.0*(1.0+dx))+0.5
   end do

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine part2d_extractpsi
!
end module part2d_lib
