! part3d_lib module for QPAD

module part3d_lib

use parallel_pipe_class
use ufield_class
use param
use system
use mpi

implicit none

character(len=10), private, save :: cls_name = 'part3d_lib'
character(len=128), private, save :: erstr
integer, private, save :: cls_level = 2

contains
!
function ranorm()
! this program calculates a random number y from a gaussian distribution
! with zero mean and unit variance, according to the method of
! mueller and box:
!    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
!    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
! where x is a random number uniformly distributed on (0,1).
! written for the ibm by viktor k. decyk, ucla
   integer, save :: r1 = 885098780, r2 = 1824280461
   integer, save :: iflg = 0, r4 = 1396483093, r5 = 55318673
   real(kind=DB), save :: h1l = 65531.0d0, h1u = 32767.0d0
   real(kind=DB), save :: h2l = 65525.0d0,r0 = 0.0d0
   real(kind=DB) :: ranorm,r3,asc,bsc,temp
   integer :: isc, i1

   if (iflg == 0) then
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
   else
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   end if
end function ranorm
!
subroutine beam_dist000(part,qm,edges,npp,dr,dz,nps,vtx,vty,vtz,vdx,&
&vdy,vdz,npx,npy,npz,n1,n2,idimp,npmax,sigx,&
&sigy,sigz,x0,y0,z0,cx,cy,lquiet,ierr)

   implicit none

   integer, intent(in) :: npx,npy,npz,idimp,n1,n2
   integer, intent(inout) :: ierr
   integer(kind=LG), intent(inout) :: npp
   integer(kind=LG), intent(in) :: npmax,nps
   real, intent(in) :: dr,dz,qm,sigx,sigy,sigz,x0,y0,z0
   real, intent(in) :: vtx,vty,vtz,vdx,vdy,vdz
   real, dimension(:,:), intent(inout) :: part
   real, dimension(4), intent(in) :: edges
   real, dimension(3), intent(in) :: cx,cy
   logical, intent(in) :: lquiet
! local data
   character(len=20), save :: sname = "beam_dist000"
   integer(kind=LG) :: i, np, npt
   integer :: j,k,l
   real :: tempx,tempy,tempr,tempxx,tempyy,x2,y2,tempz,tvtx,tvty,tvtz
   real :: borderlz, borderz, r0, sigmar0

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ierr = 0

   npt = 1

   i = 1
   
   x2 = 2.0 * x0
   y2 = 2.0 * y0
   r0 = sqrt(x0**2+y0**2)
   sigmar0 = sqrt(sigx**2+sigy**2)
   borderlz = max((z0-5.0*sigz),1.0)
   borderz = min((z0+5.0*sigz),real(n2-1))
   np = npx*npy*npz
   
   do
      if (i > np) exit
      do
         tempz = z0+sigz*ranorm()
         if (tempz < borderz .and. tempz > borderlz) then
            exit          
         end if
      end do

      tempxx = -cx(1)*(tempz-z0)**2-cx(2)*(tempz-z0)-cx(3)                
      tempyy = -cy(1)*(tempz-z0)**2-cy(2)*(tempz-z0)-cy(3)
      do
         tempx = ranorm()
         tempy = ranorm()
         if ((tempx**2+tempy**2) > 25.0) then
            cycle          
         end if
         tempx = x0 + sigx*tempx + tempxx
         tempy = y0 + sigy*tempy + tempyy
         tempr = sqrt(tempx**2+tempy**2)
         exit
      end do
      tvtx = vtx*ranorm() + vdx
      tvty = vty*ranorm() + vdy
      tvtz = vtz*ranorm() + vdz
      tvtz = sqrt(tvtz*tvtz-1-tvtx*tvtx-tvty*tvty)
      if ((tempr >= edges(1)) .and. (tempr < edges(2)) .and.&
      &(tempz >= edges(3)) .and. (tempz < edges(4))) then
         if (npt < npmax) then
            part(3,npt) = tempz
            part(1,npt) = tempr
            if (tempx == 0.0) then
               if (tempy == 0.0) part(2,npt) = 0.0
               if (tempy > 0.0) part(2,npt) = pi/2.0
               if (tempy < 0.0) part(2,npt) = -pi/2.0
            else if (tempx > 0) then
               part(2,npt) = atan(tempx/tempy)
            else
               part(2,npt) = atan(tempx/tempy) + pi
            end if
            part(4,npt) = tvtx*cos(part(2,npt))+tvty*sin(part(2,npt))
            part(5,npt) = -tvtx*sin(part(2,npt))+tvty*cos(part(2,npt))
            part(6,npt) = tvtz 
            ! part(7,npt) = qm/tempr*dr
            part(7,npt) = qm
            npt = npt + 1
         else
            ierr = ierr + 1
         end if
      end if
      i = i + 1   
      if (lquiet) then
         if (npt < npmax) then
            tempx = x2 - tempx + 2.0*tempxx
            tempy = y2 - tempy + 2.0*tempyy
            tempr = sqrt(tempx**2+tempy**2)
            if ((tempr >= edges(1)) .and. (tempr < edges(2)) .and.&
            &(tempz >= edges(3)) .and. (tempz < edges(4))) then
               part(3,npt) = tempz
               part(1,npt) = tempr
               if (tempx == 0.0) then
                  if (tempy == 0.0) part(2,npt) = 0.0
                  if (tempy > 0.0) part(2,npt) = pi/2.0
                  if (tempy < 0.0) part(2,npt) = -pi/2.0
               else if (tempx > 0) then
                  part(2,npt) = atan(tempx/tempy)
               else
                  part(2,npt) = atan(tempx/tempy) + pi
               end if
               part(4,npt) = -tvtx*cos(part(2,npt))-tvty*sin(part(2,npt))
               part(5,npt) = tvtx*sin(part(2,npt))-tvty*cos(part(2,npt))
               part(6,npt) = tvtz 
               ! part(7,npt) = qm/tempr*dr
               part(7,npt) = qm
               npt = npt + 1
            end if
         else 
            ierr = ierr + 1
         end if
         i = i + 1
      end if
   enddo
   npp = npt - 1
   call write_dbg(cls_name, sname, cls_level, 'ends')
   return
end
!
subroutine part3d_qdeposit(part,npp,dr,dz,q_re,q_im,num_modes)
! For 3D particles

   implicit none

   real, dimension(:,:), pointer, intent(in) :: part
   integer(kind=LG), intent(in) :: npp
   real, intent(in) :: dr, dz
   class(ufield), dimension(:), pointer, intent(in) :: q_re, q_im
   integer, intent(in) :: num_modes
! local data
   character(len=20), save :: sname = "part3d_qdeposit"
   integer :: i, j, nn, mm, noff1, noff2, n1p, n2p
   integer(kind=LG) :: ii
   real, dimension(:,:,:), pointer :: q0, qr, qi
   real :: r, qc, th, zz, dd, ad, zd, za, rcr, rci
   complex(kind=DB) :: rc, rc0

   call write_dbg(cls_name, sname, cls_level, 'starts')

   noff1 = q_re(0)%get_noff(1)
   noff2 = q_re(0)%get_noff(2)
   n1p = q_re(0)%get_ndp(1)
   n2p = q_re(0)%get_ndp(2)
   q0 => q_re(0)%get_f2()
   qr => null(); qi => null()

   do ii = 1, npp
      r = part(1,ii)/dr
      th = part(2,ii)
      zz = part(3,ii)/dz
      qc = part(7,ii)
      rc0 = cmplx(cos(th),-sin(th),kind=DB)
      r = r + 0.5
      nn = r
      mm = zz
      dd = qc*(r - real(nn))
      ad = qc - dd
      zd = zz - real(mm)
      za = 1.0 - zd
      nn = nn - noff1
      mm = mm - noff2 + 1
      q0(1,nn,mm) = q0(1,nn,mm) + ad*za
      q0(1,nn+1,mm) = q0(1,nn+1,mm) + dd*za
      q0(1,nn,mm+1) = q0(1,nn,mm+1) + ad*zd
      q0(1,nn+1,mm+1) = q0(1,nn+1,mm+1) + dd*zd
      rc = rc0
      do i = 1, num_modes
         rcr = real(rc)
         rci = aimag(rc)
         qr => q_re(i)%get_f2()
         qi => q_im(i)%get_f2()
         qr(1,nn,mm) = qr(1,nn,mm) + ad*za*rcr
         qr(1,nn+1,mm) = qr(1,nn+1,mm) + dd*za*rcr
         qr(1,nn,mm+1) = qr(1,nn,mm+1) + ad*zd*rcr
         qr(1,nn+1,mm+1) = qr(1,nn+1,mm+1) + dd*zd*rcr
         qi(1,nn,mm) = qi(1,nn,mm) + ad*za*rci
         qi(1,nn+1,mm) = qi(1,nn+1,mm) + dd*za*rci
         qi(1,nn,mm+1) = qi(1,nn,mm+1) + ad*zd*rci
         qi(1,nn+1,mm+1) = qi(1,nn+1,mm+1) + dd*zd*rci
         rc = rc*rc0
      end do
   end do

   if (noff1 == 0) then
      q0(1,0,:) = q0(1,0,:)/0.5
      do i = 1, num_modes
         qr => q_re(i)%get_f2()
         qi => q_im(i)%get_f2()
         qr(1,0,:) = qr(1,0,:)/0.5
         qi(1,0,:) = qi(1,0,:)/0.5
      end do      
   else
      q0(1,0,:) = q0(1,0,:)/(0.5+noff1-1)
      do i = 1, num_modes
         qr => q_re(i)%get_f2()
         qi => q_im(i)%get_f2()
         qr(1,0,:) = qr(1,0,:)/(0.5+noff1-1)
         qi(1,0,:) = qi(1,0,:)/(0.5+noff1-1)
      end do      
   end if

   do i = 1, n1p+1
      r = 0.5 + i + noff1 - 1
      q0(1,i,:) = q0(1,i,:)/r
   end do
      
   do i = 1, num_modes
      qr => q_re(i)%get_f2()
      qi => q_im(i)%get_f2()
      do j = 1, n1p+1
         r = 0.5 + j + noff1 - 1
         qr(1,j,:) = qr(1,j,:)/r
         qi(1,j,:) = qi(1,j,:)/r
      end do      
   end do
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine part3d_qdeposit
!
subroutine part3d_push(part,npp,dr,dz,xdim,dt,qbm,ef_re,ef_im,&
&bf_re,bf_im,num_modes)

   implicit none

   real, dimension(:,:), pointer, intent(inout) :: part
   integer(kind=LG), intent(inout) :: npp
   real, intent(in) :: dr, dz, dt, qbm
   class(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im, &
   &bf_re, bf_im
   integer, intent(in) :: xdim, num_modes
! local data
   character(len=20), save :: sname = "part3d_push"
   integer :: i
   integer(kind=LG) :: ii
   real, dimension(:,:,:), pointer :: e0,b0,er,ei,br,bi
   integer :: n1, n2, nn, mm, noff1, n1p, noff2, n2p
   real :: edge1, edge2, qtmh, qtmh1, qtmh2, dti
   real :: r0, r, rn, qc, qc1, th, th1, dd, ad, za, zd, rcr, rci
   real :: zz, z
   real, dimension(3) :: dx, dxx, ox, oxx, tmp
   real :: ddx, ddy, vx, vy, acx, acy, acz, omt, anorm
   real :: rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8
   real :: dtc1, v1, v2, ngamma, p2, dtx1, dtz1
   complex(kind=DB) :: rc, rc0

   call write_dbg(cls_name, sname, cls_level, 'starts')

   noff1 = ef_re(0)%get_noff(1)
   n1p = ef_re(0)%get_ndp(1)
   n1 = ef_re(0)%get_nd(1)
   noff2 = ef_re(0)%get_noff(2)
   n2p = ef_re(0)%get_ndp(2)
   n2 = ef_re(0)%get_nd(2)
   e0 => ef_re(0)%get_f2()
   b0 => bf_re(0)%get_f2()
   er => null(); ei => null(); br => null(); bi => null()

   ! idex = 1.0/dx0
   ! dtx = dt/dx0
   ! dtz = dt/dz0
   qtmh = qbm*dt
   edge1 = (real(n1) - 0.5)*dr
   edge2 = (real(n2) - 1.0)*dz

   ii = 1
   do
      if (ii > npp) exit
      r0 = part(1,ii)
      th = part(2,ii)
      z = part(3,ii)
      qc = part(7,ii)
      rc0 = cmplx(cos(th),sin(th),kind=DB)
      r = r0/dr + 0.5  
      zz = z/dz    
      nn = r
      mm = zz
      dd = r - real(nn)
      ad = 1.0 - dd
      zd = zz - real(mm)
      za = 1.0 - zd
      nn = r - noff1
      mm = zz - noff2 + 1
      dx(1:3) = ad*e0(1:3,nn,mm)
      dx(1:3) = za*(dd*e0(1:3,nn+1,mm) + dx(1:3))
      tmp(1:3) = ad*e0(1:3,nn,mm+1)
      dx(1:3) = dx(1:3) + zd*(dd*e0(1:3,nn+1,mm+1) + tmp(1:3))

      ox(1:3) = ad*b0(1:3,nn,mm)
      ox(1:3) = za*(dd*b0(1:3,nn+1,mm) + ox(1:3))
      tmp(1:3) = ad*b0(1:3,nn,mm+1)
      ox(1:3) = ox(1:3) + zd*(dd*b0(1:3,nn+1,mm+1) + tmp(1:3))
      rc = rc0
      do i = 1, num_modes
         rcr = 2.0*real(rc)
         rci = 2.0*aimag(rc)
         er => ef_re(i)%get_f2()
         ei => ef_im(i)%get_f2()
         dxx(1:3) = ad*(er(1:3,nn,mm)*rcr - ei(1:3,nn,mm)*rci)
         dxx(1:3) = za*(dd*(er(1:3,nn+1,mm)*rcr-ei(1:3,nn+1,mm)*rci) + dxx(1:3))
         tmp(1:3) = ad*(er(1:3,nn,mm+1)*rcr - ei(1:3,nn,mm+1)*rci)
         dxx(1:3) = dxx(1:3) + zd*(dd*(er(1:3,nn+1,mm+1)*rcr-ei(1:3,nn+1,mm+1)*rci) + tmp(1:3))
         dx(1:3) = dx(1:3) + dxx(1:3)
         br => bf_re(i)%get_f2()
         bi => bf_im(i)%get_f2()
         oxx(1:3) = ad*(br(1:3,nn,mm)*rcr - bi(1:3,nn,mm)*rci)
         oxx(1:3) = za*(dd*(br(1:3,nn+1,mm)*rcr - bi(1:3,nn+1,mm)*rci) + oxx(1:3))
         tmp(1:3) = ad*(br(1:3,nn,mm+1)*rcr - bi(1:3,nn,mm+1)*rci)
         oxx(1:3) = oxx(1:3) + zd*(dd*(br(1:3,nn+1,mm+1)*rcr - bi(1:3,nn+1,mm+1)*rci) + tmp(1:3))
         ox(1:3) = ox(1:3) + oxx(1:3)
         rc = rc*rc0
      end do
      dx(1) = dx(1) - ox(2)
      dx(2) = dx(2) + ox(1)
! calculate half impulse
      dx(1:3) = qtmh*dx(1:3)
! half acceleration
      acx = part(4,ii) + dx(1)
      acy = part(5,ii) + dx(2)
      acz = part(6,ii) + dx(3)
      part(6,ii) = acz
      p2 = acx**2 + acy**2
      ngamma = sqrt(1.0 + p2 + acz**2)
      dtx1 = dt/ngamma
      dtz1 = dt*(1.0+p2)/(acz*(acz+ngamma))

      v1 = r0 + acx*dtx1
      v2 = acy*dtx1
      rn = sqrt(v1*v1+v2*v2)
      if (v1 > 0) then
         th1 = th + atan(v2/v1)
      else if (v1 < 0) then
         th1 = th + atan(v2/v1) + pi
      else
         if (v2 > 0) then
            th1 = th + pi/2.0
         else if (v2 < 0) then
            th1 = th - pi/2.0
         else
            th1 = th
         end if
      end if
      acz = z + dtz1
      if (rn > edge1) then 
         if (ii == npp) then
            npp = npp -1
            exit
         else
            do i = 1, xdim
               part(i,ii) = part(i,npp)
            end do
            npp = npp - 1
            cycle
         end if
      else if (acz > edge2) then
         if (ii == npp) then
            npp = npp -1
            exit
         else
            do i = 1, xdim
               part(i,ii) = part(i,npp)
            end do
            npp = npp - 1
            cycle
         end if
      else        
         part(1,ii) = rn
         part(2,ii) = th1
         part(3,ii) = acz
         part(4,ii) = (acx*v1+acy*v2)/rn
         part(5,ii) = r0*acy/rn
         ii = ii + 1
      end if
   end do

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine part3d_push
!
subroutine part3d_pmove(part,pp,ud,npp,dr,dz,sbufr,sbufl,rbufr,rbufl,ihole,pbuff,&
&xdim,npmax,nbmax,tag1,tag2,id,info)

   implicit none

   real, dimension(:,:), pointer, intent(inout) :: part
   real, dimension(:,:), intent(inout) :: sbufl,sbufr,rbufl,rbufr,pbuff
   integer(kind=LG), intent(inout) :: npp
   real, intent(in) :: dr, dz
   integer(kind=LG), intent(in) :: npmax, nbmax
   integer(kind=LG), dimension(:), intent(inout) :: ihole
   integer, intent(in) :: xdim
   class(parallel_pipe), pointer, intent(in) :: pp
   class(ufield), intent(in) :: ud
   integer, intent(in) :: tag1, tag2
   integer, intent(inout) :: id
   integer, dimension(:), intent(inout) :: info
! local data
   character(len=20), save :: sname = "part3d_pmove"
   real, dimension(4) :: edges
   integer, dimension(2) :: noff         
   integer, dimension(2) :: jsr, jsl, jss
   integer :: n1, n2, kstrt, nvpy, nvpz
   integer :: nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
   integer, parameter :: iy = 1, iz = 3
   integer :: ierr, ic, js, ks, mnblok, i, n, m, my, mz, moff, nvp, iter
   integer :: npr, nps, npt, kl, kr, j, j1, j2, nter, mter
   integer :: nbsize
   integer :: itermax
   integer, dimension(10) :: istatus
   integer, dimension(4) :: ibflg, iwork, msid
   integer, dimension(2) :: kb
   real(kind=DB), dimension(2) :: bflg, work
   real :: an, xt

   ierr = 0
   n1 = ud%get_nd(1); n2 = ud%get_nd(2)
   noff = ud%get_noff()
   if (noff(1) == 0) then
      edges(1) = noff(1)*dr
      edges(2) = edges(1) + (ud%get_ndp(1) + 0.5)*dr
   else
      edges(1) = (noff(1) + 0.5)*dr
      edges(2) = edges(1) + ud%get_ndp(1)*dr
   end if
   edges(3) = noff(2)*dz
   edges(4) = edges(3) + ud%get_ndp(2)*dz
   kstrt = pp%getkstrt()
   nvpy = pp%getlnvp()
   nvpz = pp%getnstage()
   mreal = pp%getmreal()
   mint = pp%getmint()
   lworld = pp%getlworld()
   lgrp = pp%getlgrp()

   ks = (kstrt - 1)/nvpy
   js = kstrt - nvpy*ks - 2
   ks = ks - 1
   nbsize = xdim*nbmax
   do j = 1, 9
      info(j) = 0
   end do
!buffer outgoing particles, first in y then in z direction
   ic = iz
   nvp = nvpy*nvpz
   an = float(n2)
   n = 2
   kl = kstrt - nvpy
   if (kl >= 1) then
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,tag1,lworld,msid(1),ierr)
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2) = nps/xdim
      jss(2) = npp + jsl(2)
      if (jss(2) <= npmax) then
         do j = 1, jsl(2)
            do i = 1, xdim
               part(i,j+npp) = rbufl(i,j)
            end do
         end do
         npp = jss(2)
      else
         write (erstr,*) 'particle overflow', jss(2)
         call write_dbg(cls_name, sname, cls_level, erstr)
         info(1) = jss(2)
         return
      end if
   end if
  
   iter = 2
   nter = 0
   mter = 0
   jsl(1) = 0
   jsr(1) = 0
   jss(2) = 0
   do j = 1, npp
      xt = part(ic,j)
! particles going backward, not going to happen
      if (xt < edges(2*n-1)) then
         jss(2) = 1
         write (erstr,*) 'Error: particles move to the previous stage'
         call write_dbg(cls_name, sname, cls_level, erstr)
         exit
! particles going forward
      else if (xt >= edges(2*n)) then
         if (jsr(1) < nbmax) then
            jsr(1) = jsr(1) + 1
            do i = 1, xdim
               pbuff(i,jsr(1)) = part(i,j)
            end do
            ihole(jsr(1)) = j
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
   if (nps > 0) then
      info(1) = nps
      write (erstr,*) 'particle buffer full'
      call write_err(cls_name//sname//erstr)
      return
   end if
! copy particle buffers
   iter = iter + 2
   mter = mter + 1
   kr = kstrt + nvpy
   if (kr <= nvp) then
! send particles
      call MPI_ISEND(pbuff,xdim*jsr(1),mreal,kr-1,tag2,lworld,id,ierr)
   end if
! fill the holes
   npt = npp
   do j = 1, jsr(1)
      j1 = ihole((jsr(1)-j+1))
      if (j1 < npt) then
         do i = 1, xdim
            part(i,j1) = part(i,npt)
         end do
         npt = npt - 1
      else
         npt = npt - 1
      end if
   end do
   npp = npt
   itermax = 20000
! buffer outgoing particles, first in y then in z direction
   n = 1
   ic = iy
   nvp = nvpy
   an = float(n1)
   iter = 2
   nter = 0
   do
   mter = 0
   kb(2) = 1 + ks
   kb(1) = 1 + js
   jsl(1) = 0
   jsr(1) = 0
   jss(2) = 0
   do j = 1, npp
      xt = part(ic,j)
! particles going down or backward
      if (xt < edges(2*n-1)) then
         if (jsl(1) < nbmax) then
            jsl(1) = jsl(1) + 1
            do i = 1, xdim
               sbufl(i,jsl(1)) = part(i,j)
            end do
            ihole(jsl(1)+jsr(1)) = j
         else
            jss(2) = 1
            exit
         end if
! particles going up or forward
      else if (xt >= edges(2*n)) then
         if (jsr(1) < nbmax) then
            jsr(1) = jsr(1) + 1
            do i = 1, xdim
               sbufr(i,jsr(1)) = part(i,j)
            end do
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
   kb(1) = 1 + js
   kb(2) = 1 + ks
! get particles from below and above or back and front
   kl = kb(n)
   kb(n) = kl + 1
   if (kb(n) >= nvp) kb(n) = kb(n) - nvp
   kr = kb(1) + nvpy*kb(2) + 1
   kb(n) = kl - 1
   if (kb(n) < 0) kb(n) = kb(n) + nvp
   kl = kb(1) + nvpy*kb(2) + 1
! post receive
   call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lworld,msid(1),ierr)
   call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lworld,msid(2),ierr)
! send particles
   call MPI_ISEND(sbufr,xdim*jsr(1),mreal,kr-1,iter-1,lworld,msid(3),ierr)
   call MPI_ISEND(sbufl,xdim*jsl(1),mreal,kl-1,iter,lworld,msid(4),ierr)
! wait for particles to arrive
   call MPI_WAIT(msid(1),istatus,ierr)
   call MPI_GET_COUNT(istatus,mreal,nps,ierr)
   jsl(2) = nps/xdim
   call MPI_WAIT(msid(2),istatus,ierr)
   call MPI_GET_COUNT(istatus,mreal,nps,ierr)
   jsr(2) = nps/xdim
! check if particles must be passed further
   nps = 0
! check if any particles coming from above or front belong here
   jsl(1) = 0
   jsr(1) = 0
   jss(2) = 0
   do j = 1, jsr(2)
      if (rbufr(ic,j) < edges(2*n-1)) jsl(1) = jsl(1) + 1
      if (rbufr(ic,j) >= edges(2*n)) jsr(1) = jsr(1) + 1
   end do
   if (jsr(1) /= 0) then
      write (erstr,*) 'Info:',jsr(1),' particles returning above'
      call write_dbg(cls_name, sname, cls_level, erstr)
   end if
! check if any particles coming from below or back belong here
   do j = 1, jsl(2)
      if (rbufl(ic,j) >= edges(2*n)) jsr(1) = jsr(1) + 1
      if (rbufl(ic,j) < edges(2*n-1)) jss(2) = jss(2) + 1
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
      kb(2) = 1 + ks
      kb(1) = 1 + js
! first check particles coming from above or front
      jsl(1) = 0
      jsr(1) = 0
      jss(2) = 0
      do j = 1, jsr(2)
         xt = rbufr(ic,j)
! particles going down or back
         if (xt < edges(2*n-1)) then
            jsl(1) = jsl(1) + 1
            do i = 1, xdim
               sbufl(i,jsl(1)) = rbufr(i,j)
            end do
! particles going up or front, should not happen
         else if (xt >= edges(2*n)) then
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
! next check particles coming from below or back
      jss(2) = 0
      do j = 1, jsl(2)
         xt = rbufl(ic,j)
! particles going up or front
         if (xt >= edges(2*n)) then
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
         elseif (xt < edges(2*n-1)) then
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
       info(1) = ierr
       return
    end if
! distribute incoming particles from buffers
! distribute particles coming from below or back into holes
      jss(2) = min0(jss(1),jsl(2))
      do j = 1, jss(2)
         do i = 1, xdim
            part(i,ihole(j)) = rbufl(i,j)
         end do
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
            do i = 1, xdim
               part(i,ihole(j+jsl(2))) = rbufr(i,j)
            end do
         else
! no more holes
! distribute remaining particles from below or back into bottom
            do i = 1, xdim
               part(i,j+npp) = rbufl(i,j+jss(1))
            end do
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
               do i = 1, xdim
                  part(i,ihole(j2)) = part(i,j1)
               end do
            end if
         else
! no more holes
! distribute remaining particles from above or front into bottom
            do i = 1, xdim
               part(i,j+npp) = rbufr(i,j+jss(1))
            end do
         end if
      end do
      if (jss(2) > 0) then
         npp = npp - jsr(2)
      else
         npp = npp + jsr(2)
      end if
      jss(1) = 0
! check if any particles have to be passed further
      info(5+n) = max0(info(5+n),mter)
      if (ibflg(2) <= 0) exit
      write (erstr,*) 'Info: particles being passed further = ', ibflg(2)
      call write_dbg(cls_name, sname, cls_level, erstr)
      if (ibflg(3) > 0) ibflg(3) = 1
      if (iter >= itermax) then
         ierr = -((iter-2)/2)
         write (erstr,*) 'Iteration overflow, iter = ', ierr
         call write_dbg(cls_name, sname, cls_level, erstr)
         info(1) = ierr
         if (nter > 0) then
            write (erstr,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
            call write_dbg(cls_name, sname, cls_level, erstr)
         end if
         call write_dbg(cls_name, sname, cls_level, 'ends')
         return
      end if
   end do
! check if buffer overflowed and more particles remain to be checked
      if (ibflg(3) <= 0) exit
      nter = nter + 1
      info(3+n) = nter
      write (erstr,*) "new loop, nter=", nter 
      call write_dbg(cls_name, sname, cls_level, erstr)
  end do
! information
   if (nter > 0) then
      write (erstr,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      call write_dbg(cls_name, sname, cls_level, erstr)
   end if
   call write_dbg(cls_name, sname, cls_level, 'ends')
   return
end subroutine part3d_pmove
!
end module part3d_lib