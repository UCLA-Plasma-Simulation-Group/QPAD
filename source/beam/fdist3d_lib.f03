module fdist3d_lib

use sysutil_module
use random
use param

implicit none

character(len=10), private, save :: cls_name = 'fdist3d_lib'
integer, private, save :: cls_level = 2

contains

!
subroutine beam_dist000(x,p,q,qm,edges,npp,dr,dz,nps,vtx,vty,vtz,vdx,&
&vdy,vdz,npx,npy,npz,rmax,zmin,zmax,npmax,sigx,&
&sigy,sigz,x0,y0,z0,cx,cy,lquiet,ierr)

   implicit none

   integer, intent(in) :: npx,npy,npz
   integer, intent(inout) :: ierr
   integer(kind=LG), intent(inout) :: npp
   integer(kind=LG), intent(in) :: npmax,nps
   real, intent(in) :: dr,dz,qm,sigx,sigy,sigz,x0,y0,z0
   real, intent(in) :: vtx,vty,vtz,vdx,vdy,vdz,rmax,zmin,zmax
   ! real, dimension(:,:), intent(inout) :: part
   real, dimension(:,:), intent(inout) :: x, p
   real, dimension(:), intent(inout) :: q
   real, dimension(4), intent(in) :: edges
   real, dimension(3), intent(in) :: cx,cy
   logical, intent(in) :: lquiet
! local data
   character(len=20), save :: sname = "beam_dist000"
   integer(kind=LG) :: i, np, npt
   real :: tempx,tempy,tempr,tempxx,tempyy,tempz,tvtx,tvty,tvtz
   real :: borderlz, borderz, sigmar0, trunc = 5.0

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ierr = 0

   npt = 1

   i = 1

   sigmar0 = sqrt(sigx**2+sigy**2)
   borderlz = max((z0-trunc*sigz),zmin)
   borderz = min((z0+trunc*sigz),zmax)
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
         if ((tempx**2+tempy**2) > trunc**2) then
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
            x(1,npt) = tempx
            x(2,npt) = tempy
            x(3,npt) = tempz
            p(1,npt) = tvtx
            p(2,npt) = tvty
            p(3,npt) = tvtz
            q(npt) = qm
            npt = npt + 1
         else
            ierr = ierr + 1
         end if
      end if
      i = i + 1
      if (lquiet) then
         if (npt < npmax) then
            tempx = 2.0 * x0 - tempx + 2.0*tempxx
            tempy = 2.0 * y0 - tempy + 2.0*tempyy
            tempr = sqrt(tempx**2+tempy**2)
            if ((tempr >= edges(1)) .and. (tempr < edges(2)) .and.&
            &(tempz >= edges(3)) .and. (tempz < edges(4))) then
               x(1,npt) = tempx
               x(2,npt) = tempy
               x(3,npt) = tempz
               p(1,npt) = -tvtx
               p(2,npt) = -tvty
               p(3,npt) =  tvtz
               q(npt) = qm
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
end subroutine beam_dist000

subroutine beam_dist001(x,p,q,qm,edges,npp,dr,dz,nps,vtx,vty,vtz,vdx,&
&vdy,vdz,npr,npth,npz,rmax,zmin,zmax,npmax,sigx,&
&sigy,sigz,x0,y0,z0,cx,cy,lquiet,ierr)

   implicit none

   integer, intent(in) :: npr,npth,npz
   integer, intent(inout) :: ierr
   integer(kind=LG), intent(inout) :: npp
   integer(kind=LG), intent(in) :: npmax,nps
   real, intent(in) :: dr,dz,qm,sigx,sigy,sigz,x0,y0,z0
   real, intent(in) :: vtx,vty,vtz,vdx,vdy,vdz,rmax,zmin,zmax
   real, dimension(:,:), intent(inout) :: x, p
   real, dimension(:), intent(inout) :: q
   real, dimension(4), intent(in) :: edges
   real, dimension(3), intent(in) :: cx,cy
   logical, intent(in) :: lquiet
! local data
   character(len=20), save :: sname = "beam_dist001"
   integer(kind=LG) :: i, np, npt
   real :: tempr,tempth,tempx,tempy,tempxx,tempyy,tempz,tvtx,tvty,tvtz
   real :: borderlz, borderz, qm_amp, sigr, trunc = 5.0

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ierr = 0

   npt = 1

   i = 1

   borderlz = max((z0-trunc*sigz),zmin)
   borderz  = min((z0+trunc*sigz),zmax)
   np = npr*npth*npz

   sigr = max( sigx, sigy )

   do
      if (i > np) exit
      do
         ! tempz = z0+sigz*ranorm()
         call random_number(tempz)
         tempz = z0 + 2.0 * trunc * sigz * (tempz-0.5)
         if (tempz < borderz .and. tempz > borderlz) then
            exit
         end if
      end do

      tempxx = -cx(1)*(tempz-z0)**2-cx(2)*(tempz-z0)-cx(3)
      tempyy = -cy(1)*(tempz-z0)**2-cy(2)*(tempz-z0)-cy(3)

      call random_number(tempth)
      tempth = tempth * 2.0 * pi
      call random_number(tempr)
      tempr  = trunc * sigr * tempr
      tempx  = tempr * cos(tempth)
      tempy  = tempr * sin(tempth)
      qm_amp = trunc**2 * sigr * sigz * tempr*&
         exp( -0.5 * ( (tempx/sigx)**2 + (tempy/sigy)**2 + ((tempz-z0)/sigz)**2 ) )

      ! add offset
      tempx = tempx + x0 + tempxx
      tempy = tempy + y0 + tempyy
      tempr  = sqrt( tempx**2 + tempy**2 )

      tvtx = vtx*ranorm() + vdx
      tvty = vty*ranorm() + vdy
      tvtz = vtz*ranorm() + vdz
      tvtz = sqrt(tvtz*tvtz-1-tvtx*tvtx-tvty*tvty)

      if ( (tempr >= edges(1)) .and. (tempr < edges(2)) .and.&
           (tempz >= edges(3)) .and. (tempz < edges(4)) ) then
         if (npt < npmax) then
            x(1,npt) = tempx
            x(2,npt) = tempy
            x(3,npt) = tempz
            p(1,npt) = tvtx
            p(2,npt) = tvty
            p(3,npt) = tvtz
            q(npt) = qm * qm_amp
            npt = npt + 1
         else
            ierr = ierr + 1
         end if
      end if
      i = i + 1
      if (lquiet) then
         if (npt < npmax) then
            tempx = 2.0 * x0 - tempx + 2.0*tempxx
            tempy = 2.0 * y0 - tempy + 2.0*tempyy
            tempr = sqrt(tempx**2+tempy**2)
            if ((tempr >= edges(1)) .and. (tempr < edges(2)) .and.&
                (tempz >= edges(3)) .and. (tempz < edges(4))) then
               x(1,npt) = tempx
               x(2,npt) = tempy
               x(3,npt) = tempz
               p(1,npt) = -tvtx
               p(2,npt) = -tvty
               p(3,npt) =  tvtz
               q(npt) = qm * qm_amp
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
end subroutine beam_dist001
!

subroutine beam_dist002(x,p,q,qm,edges,npp,dr,dz,nps,vtx,vty,vtz,vdx,&
&vdy,vdz,npx,npy,npz,rmax,zmin,zmax,npmax,sigx,&
&sigy,x0,y0,z0,cx,cy,zf,lquiet,ierr)

   implicit none

   integer, intent(in) :: npx,npy,npz
   integer, intent(inout) :: ierr
   integer(kind=LG), intent(inout) :: npp
   integer(kind=LG), intent(in) :: npmax,nps
   real, intent(in) :: dr,dz,qm,sigx,sigy,x0,y0,z0
   real, intent(in) :: vtx,vty,vtz,vdx,vdy,vdz,rmax,zmin,zmax
   real, dimension(:,:), intent(inout) :: x, p
   real, dimension(:), intent(inout) :: q
   real, dimension(4), intent(in) :: edges
   real, dimension(3), intent(in) :: cx,cy
   real, dimension(:), intent(inout) :: zf
   logical, intent(in) :: lquiet
! local data
   character(len=20), save :: sname = "beam_dist002"
   integer(kind=LG) :: i, np, npt
   real :: tempx,tempy,tempr,tempxx,tempyy,tempz,tvtx,tvty,tvtz
   real :: tag, sigmar0, trunc = 5.0
   integer :: nz, j

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ierr = 0

   npt = 1

   i = 1

   sigmar0 = sqrt(sigx**2+sigy**2)
   np = npx*npy*npz

   nz = size(zf)
   zf = abs(zf)
   do j = 2, nz
      zf(j) = zf(j-1) + zf(j)
   enddo
   zf = zf - zf(1)
   zf = zf / zf(nz)

   do while (i <= np)

      call random_number(tag)
      j = 1
      do while (zf(j) <= tag)
         j = j + 1
      enddo
      tempz = ( real(j-1) + (tag - zf(j-1)) / (zf(j) - zf(j-1)) ) * dz + zmin

      tempxx = -cx(1)*(tempz-z0)**2-cx(2)*(tempz-z0)-cx(3)
      tempyy = -cy(1)*(tempz-z0)**2-cy(2)*(tempz-z0)-cy(3)
      do
         tempx = ranorm()
         tempy = ranorm()
         if ((tempx**2+tempy**2) > trunc**2) then
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
            x(1,npt) = tempx
            x(2,npt) = tempy
            x(3,npt) = tempz
            p(1,npt) = tvtx
            p(2,npt) = tvty
            p(3,npt) = tvtz
            q(npt) = qm
            npt = npt + 1
         else
            ierr = ierr + 1
         end if
      end if
      i = i + 1
      if (lquiet) then
         if (npt < npmax) then
            tempx = 2.0 * x0 - tempx + 2.0*tempxx
            tempy = 2.0 * y0 - tempy + 2.0*tempyy
            tempr = sqrt(tempx**2+tempy**2)
            if ((tempr >= edges(1)) .and. (tempr < edges(2)) .and.&
            &(tempz >= edges(3)) .and. (tempz < edges(4))) then
               x(1,npt) = tempx
               x(2,npt) = tempy
               x(3,npt) = tempz
               p(1,npt) = -tvtx
               p(2,npt) = -tvty
               p(3,npt) =  tvtz
               q(npt) = qm
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
end subroutine beam_dist002

end module fdist3d_lib
