module part3d_comm

use sysutil
use param
use mpi
use parallel_pipe_class
use ufield_class

implicit none

character(len=10), private, save :: cls_name = 'part3d'
integer, private, save :: cls_level = 2
character(len=128), private :: erstr

contains

subroutine pmove_part3d(x,p,q,pp,ud,npp,dr,dz,sbufr,sbufl,rbufr,rbufl,ihole,pbuff,&
 xdim,npmax,nbmax,tag1,tag2,id,info)

   implicit none

   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:), pointer, intent(inout) :: q
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
   character(len=20), save :: sname = "pmove_part3d"
   real, dimension(4) :: edges
   integer, dimension(2) :: noff
   integer, dimension(2) :: jsr, jsl, jss
   integer :: n1, n2, kstrt, nvpy, nvpz
   integer :: lgrp, mreal, mint, lworld
   integer, parameter :: iy = 1, iz = 3
   integer :: ierr, ic, js, ks, i, n, nvp, iter
   integer :: nps, npt, kl, kr, j, j1, j2, nter, mter
   integer :: nbsize
   integer :: itermax
   integer, dimension(10) :: istatus
   integer, dimension(4) :: ibflg, iwork, msid
   integer, dimension(2) :: kb
   real :: an, xt

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'move 3D particles' )

   ierr = 0
   n1 = ud%get_nd(1); n2 = ud%get_nd(2)
   noff = ud%get_noff()
   if (noff(1) == 0) then
      edges(1) = 0.0
      edges(2) = edges(1) + ud%get_ndp(1) * dr
   else
      edges(1) = noff(1) * dr
      edges(2) = edges(1) + ud%get_ndp(1) * dr
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
            ! do i = 1, xdim
            !    part(i,j+npp) = rbufl(i,j)
            ! end do
            x(1,j+npp) = rbufl(1,j)
            x(2,j+npp) = rbufl(2,j)
            x(3,j+npp) = rbufl(3,j)
            p(1,j+npp) = rbufl(4,j)
            p(2,j+npp) = rbufl(5,j)
            p(3,j+npp) = rbufl(6,j)
            q(j+npp)   = rbufl(7,j)
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
      ! xt = part(ic,j)
      xt = x(3,j)
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
            ! do i = 1, xdim
            !    pbuff(i,jsr(1)) = part(i,j)
            ! end do
            pbuff(1,jsr(1)) = x(1,j)
            pbuff(2,jsr(1)) = x(2,j)
            pbuff(3,jsr(1)) = x(3,j)
            pbuff(4,jsr(1)) = p(1,j)
            pbuff(5,jsr(1)) = p(2,j)
            pbuff(6,jsr(1)) = p(3,j)
            pbuff(7,jsr(1)) = q(j)
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
         ! do i = 1, xdim
         !    part(i,j1) = part(i,npt)
         ! end do
         x(1,j1) = x(1,npt)
         x(2,j1) = x(2,npt)
         x(3,j1) = x(3,npt)
         p(1,j1) = p(1,npt)
         p(2,j1) = p(2,npt)
         p(3,j1) = p(3,npt)
         q(j1)   = q(npt)
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
      ! xt = part(ic,j)
      xt = sqrt( x(1,j)*x(1,j) + x(2,j)*x(2,j) )
! particles going down or backward
      if (xt < edges(2*n-1)) then
         if (jsl(1) < nbmax) then
            jsl(1) = jsl(1) + 1
            ! do i = 1, xdim
            !    sbufl(i,jsl(1)) = part(i,j)
            ! end do
            sbufl(1,jsl(1)) = x(1,j)
            sbufl(2,jsl(1)) = x(2,j)
            sbufl(3,jsl(1)) = x(3,j)
            sbufl(4,jsl(1)) = p(1,j)
            sbufl(5,jsl(1)) = p(2,j)
            sbufl(6,jsl(1)) = p(3,j)
            sbufl(7,jsl(1)) = q(j)
            ihole(jsl(1)+jsr(1)) = j
         else
            jss(2) = 1
            exit
         end if
! particles going up or forward
      else if (xt >= edges(2*n)) then
         if (jsr(1) < nbmax) then
            jsr(1) = jsr(1) + 1
            ! do i = 1, xdim
            !    sbufr(i,jsr(1)) = part(i,j)
            ! end do
            sbufr(1,jsr(1)) = x(1,j)
            sbufr(2,jsr(1)) = x(2,j)
            sbufr(3,jsr(1)) = x(3,j)
            sbufr(4,jsr(1)) = p(1,j)
            sbufr(5,jsr(1)) = p(2,j)
            sbufr(6,jsr(1)) = p(3,j)
            sbufr(7,jsr(1)) = q(j)
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
      ! if (rbufr(ic,j) < edges(2*n-1)) jsl(1) = jsl(1) + 1
      ! if (rbufr(ic,j) >= edges(2*n)) jsr(1) = jsr(1) + 1
      xt = sqrt( rbufr(1,j)*rbufr(1,j) + rbufr(2,j)* rbufr(2,j) )
      if (xt < edges(2*n-1)) jsl(1) = jsl(1) + 1
      if (xt >= edges(2*n)) jsr(1) = jsr(1) + 1
   end do
   if (jsr(1) /= 0) then
      write (erstr,*) 'Info:',jsr(1),' particles returning above'
      call write_dbg(cls_name, sname, cls_level, erstr)
   end if
! check if any particles coming from below or back belong here
   do j = 1, jsl(2)
      ! if (rbufl(ic,j) >= edges(2*n)) jsr(1) = jsr(1) + 1
      ! if (rbufl(ic,j) < edges(2*n-1)) jss(2) = jss(2) + 1
      xt = sqrt( rbufl(1,j)*rbufl(1,j) + rbufl(2,j)* rbufl(2,j) )
      if (xt >= edges(2*n)) jsr(1) = jsr(1) + 1
      if (xt < edges(2*n-1)) jss(2) = jss(2) + 1
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
         ! xt = rbufr(ic,j)
         xt = sqrt( rbufr(1,j)*rbufr(1,j) + rbufr(2,j)*rbufr(2,j) )
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
         ! xt = rbufl(ic,j)
         xt = sqrt( rbufl(1,j)*rbufl(1,j) + rbufl(2,j)*rbufl(2,j) )
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
         ! do i = 1, xdim
         !    part(i,ihole(j)) = rbufl(i,j)
         ! end do
         x(1,ihole(j)) = rbufl(1,j)
         x(2,ihole(j)) = rbufl(2,j)
         x(3,ihole(j)) = rbufl(3,j)
         p(1,ihole(j)) = rbufl(4,j)
         p(2,ihole(j)) = rbufl(5,j)
         p(3,ihole(j)) = rbufl(6,j)
         q(ihole(j))   = rbufl(7,j)
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
            x(1,ihole(j+jsl(2))) = rbufr(1,j)
            x(2,ihole(j+jsl(2))) = rbufr(2,j)
            x(3,ihole(j+jsl(2))) = rbufr(3,j)
            p(1,ihole(j+jsl(2))) = rbufr(4,j)
            p(2,ihole(j+jsl(2))) = rbufr(5,j)
            p(3,ihole(j+jsl(2))) = rbufr(6,j)
            q(ihole(j+jsl(2)))   = rbufr(7,j)
         else
! no more holes
! distribute remaining particles from below or back into bottom
            ! do i = 1, xdim
            !    part(i,j+npp) = rbufl(i,j+jss(1))
            ! end do
            x(1,j+npp) = rbufl(1,j+jss(1))
            x(2,j+npp) = rbufl(2,j+jss(1))
            x(3,j+npp) = rbufl(3,j+jss(1))
            p(1,j+npp) = rbufl(4,j+jss(1))
            p(2,j+npp) = rbufl(5,j+jss(1))
            p(3,j+npp) = rbufl(6,j+jss(1))
            q(j+npp)   = rbufl(7,j+jss(1))
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
               !    part(i,ihole(j2)) = part(i,j1)
               ! end do
               x(1,ihole(j2)) = x(1,j1)
               x(2,ihole(j2)) = x(2,j1)
               x(3,ihole(j2)) = x(3,j1)
               p(1,ihole(j2)) = p(1,j1)
               p(2,ihole(j2)) = p(2,j1)
               p(3,ihole(j2)) = p(3,j1)
               q(ihole(j2))   = q(j1)  
            end if
         else
! no more holes
! distribute remaining particles from above or front into bottom
            ! do i = 1, xdim
            !    part(i,j+npp) = rbufr(i,j+jss(1))
            ! end do
            x(1,j+npp) = rbufr(1,j+jss(1))
            x(2,j+npp) = rbufr(2,j+jss(1))
            x(3,j+npp) = rbufr(3,j+jss(1))
            p(1,j+npp) = rbufr(4,j+jss(1))
            p(2,j+npp) = rbufr(5,j+jss(1))
            p(3,j+npp) = rbufr(6,j+jss(1))
            q(j+npp)   = rbufr(7,j+jss(1))
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

   call stop_tprof( 'move 3D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')
   return
end subroutine pmove_part3d

end module part3d_comm