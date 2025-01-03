! fdist3d class for QPAD

module fdist3d_class

use parallel_pipe_class
use field_class
use ufield_class
use input_class
use part3d_lib
use param
use sysutil


implicit none

private

public :: fdist3d, fdist3d_wrap, fdist3d_000, fdist3d_001, fdist3d_002

type, abstract :: fdist3d

   private

   class(parallel_pipe), pointer, public :: p => null()
!
! ndprof = profile type
   integer :: npf, npmax
   logical :: evol = .true.
   real :: dx, dz
   real :: z0

   contains

   generic :: new => init_fdist3d
   generic :: del => end_fdist3d
   generic :: dist => dist3d
   generic :: dp => deposit_fdist3d
   procedure(ab_init_fdist3d), deferred, private :: init_fdist3d
   procedure, private :: end_fdist3d
   procedure(ab_dist3d), deferred, private :: dist3d
   procedure(ab_deposit_fdist3d), deferred, private :: deposit_fdist3d
   procedure :: getnpf, getnpmax, getevol, getdx, getdz, getz0

end type

abstract interface
!
subroutine ab_dist3d(this,part3d,npp,noff,ndp)
   import fdist3d
   import ufield
   import LG
   implicit none
   class(fdist3d), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: part3d
   integer(kind=LG), intent(inout) :: npp
   ! class(ufield), intent(in), pointer :: ud
   integer, intent(in), dimension(2) :: noff, ndp
end subroutine ab_dist3d
!
subroutine ab_init_fdist3d(this,input,i)
   import fdist3d
   import input_json
   implicit none
   class(fdist3d), intent(inout) :: this
   type(input_json), intent(inout), pointer :: input
   integer, intent(in) :: i
end subroutine ab_init_fdist3d
!
subroutine ab_deposit_fdist3d(this,q)
   import fdist3d
   import field
   implicit none
   class(fdist3d), intent(inout) :: this
   class(field), intent(inout) :: q
end subroutine ab_deposit_fdist3d
!
end interface

type fdist3d_wrap
   class(fdist3d), allocatable :: p
end type fdist3d_wrap
!
type, extends(fdist3d) :: fdist3d_000
! Tri Gaussian profile (the same particle charge)
   private

   integer :: npx, npy, npz
   real :: qm, sigx, sigy, sigz,rmax,zmin,zmax
   real :: bcx, bcy, bcz, sigvx, sigvy, sigvz
   real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
   logical :: quiet

   contains

   procedure, private :: init_fdist3d => init_fdist3d_000
   procedure, private :: deposit_fdist3d => deposit_fdist3d_000
   procedure, private :: dist3d => dist3d_000

end type fdist3d_000
!
type, extends(fdist3d) :: fdist3d_001
! Tri Gaussian profile (varying particle charge distributed in r and theta)
   private

   integer :: npr, npth, npz
   real :: qm, sigx, sigy, sigz, rmax, zmin, zmax
   real :: bcx, bcy, bcz, sigvx, sigvy, sigvz
   real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
   logical :: quiet

   contains

   procedure, private :: init_fdist3d => init_fdist3d_001
   procedure, private :: deposit_fdist3d => deposit_fdist3d_001
   procedure, private :: dist3d => dist3d_001

end type fdist3d_001

type, extends(fdist3d) :: fdist3d_002
! Bi Gaussian and piecewise in z (the same particle charge)
   private

   integer :: npx, npy, npz
   real :: qm, sigx, sigy, rmax, zmin, zmax
   real :: bcx, bcy, bcz, sigvx, sigvy, sigvz
   real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
   real, dimension(:), allocatable :: fz, z
   logical :: quiet

   contains

   procedure, private :: init_fdist3d => init_fdist3d_002
   procedure, private :: deposit_fdist3d => deposit_fdist3d_002
   procedure, private :: dist3d => dist3d_002

end type fdist3d_002

character(len=10), save :: cls_name = 'fdist3d'
integer, save :: cls_level = 2
character(len=128), save :: erstr

contains
!
function getnpf(this)

   implicit none

   class(fdist3d), intent(in) :: this
   integer :: getnpf

   getnpf = this%npf

end function getnpf
!
function getnpmax(this)

   implicit none

   class(fdist3d), intent(in) :: this
   integer(kind=LG) :: getnpmax

   getnpmax = this%npmax

end function getnpmax
!
function getevol(this)

   implicit none

   class(fdist3d), intent(in) :: this
   logical :: getevol

   getevol = this%evol

end function getevol
!
function getdx(this)

   implicit none

   class(fdist3d), intent(in) :: this
   real :: getdx

   getdx = this%dx

end function getdx
!
function getdz(this)

   implicit none

   class(fdist3d), intent(in) :: this
   real :: getdz

   getdz = this%dz

end function getdz
!
function getz0(this)

   implicit none

   class(fdist3d), intent(in) :: this
   real :: getz0

   getz0 = this%z0

end function getz0
!
subroutine end_fdist3d(this)

   implicit none

   class(fdist3d), intent(inout) :: this
   character(len=18), save :: sname = 'end_fdist3d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_fdist3d
!
subroutine init_fdist3d_000(this,input,i)

   implicit none

   class(fdist3d_000), intent(inout) :: this
   type(input_json), intent(inout), pointer :: input
   integer, intent(in) :: i
! local data
   integer :: npf,npx,npy,npz,npmax
   real :: qm,sigx,sigy,sigz,bcx,bcy,bcz,sigvx,sigvy,sigvz
   real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
   logical :: quiet, evol
   real :: min, max, cwp, n0
   real :: alx, alz, dr, dz
   integer :: nr, nz, num_modes
   character(len=20) :: sn,s1
   character(len=18), save :: sname = 'init_fdist3d_000:'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%p => input%pp

   write (sn,'(I3.3)') i
   s1 = 'beam('//trim(sn)//')'
   call input%get('simulation.n0',n0)
   call input%get('simulation.grid(1)',nr)
   call input%get('simulation.grid(2)',nz)
   call input%get('simulation.max_mode',num_modes)
   cwp=5.32150254*1e9/sqrt(n0)
   call input%get('simulation.box.r(1)',min)
   call input%get('simulation.box.r(2)',max)
   call input%get(trim(s1)//'.center(1)',bcx)
   alx = (max-min)
   dr=alx/real(nr)
   this%rmax = max
   call input%get(trim(s1)//'.center(2)',bcy)
   call input%get('simulation.box.z(1)',min)
   call input%get('simulation.box.z(2)',max)
   call input%get(trim(s1)//'.center(3)',bcz)
   bcz = bcz -min
   alz = (max-min)
   dz=alz/real(nz)
   this%z0 = min
   this%zmin = 0.0
   this%zmax = max-min
   call input%get(trim(s1)//'.profile',npf)
   call input%get(trim(s1)//'.np(1)',npx)
   call input%get(trim(s1)//'.np(2)',npy)
   call input%get(trim(s1)//'.np(3)',npz)
   call input%get(trim(s1)//'.q',qm)
   call input%get(trim(s1)//'.sigma(1)',sigx)
   call input%get(trim(s1)//'.sigma(2)',sigy)
   call input%get(trim(s1)//'.sigma(3)',sigz)
   call input%get(trim(s1)//'.sigma_v(1)',sigvx)
   call input%get(trim(s1)//'.sigma_v(2)',sigvy)
   call input%get(trim(s1)//'.sigma_v(3)',sigvz)
   call input%get(trim(s1)//'.centroid_x(1)',cx1)
   call input%get(trim(s1)//'.centroid_x(2)',cx2)
   call input%get(trim(s1)//'.centroid_x(3)',cx3)
   call input%get(trim(s1)//'.centroid_y(1)',cy1)
   call input%get(trim(s1)//'.centroid_y(2)',cy2)
   call input%get(trim(s1)//'.centroid_y(3)',cy3)
   call input%get(trim(s1)//'.quiet_start',quiet)
   call input%get(trim(s1)//'.gamma',gamma)
   call input%get(trim(s1)//'.peak_density',np)
   call input%get(trim(s1)//'.npmax',npmax)
   call input%get(trim(s1)//'.evolution',evol)
   this%npf = npf
   this%dx = dr
   this%dz = dz
   this%npx = npx
   this%npy = npy
   this%npz = npz
   this%npmax = npmax
   qm = qm/abs(qm)*np*(2*pi)**1.5*sigx*sigy*sigz
   qm = qm/dr/dz/dr
   qm = qm/npx
   qm = qm/npy
   qm = qm/npz
   this%qm = qm
   this%bcx = bcx
   this%bcy = bcy
   this%bcz = bcz
   this%sigx = sigx
   this%sigy = sigy
   this%sigz = sigz
   this%sigvx = sigvx
   this%sigvy = sigvy
   this%sigvz = sigvz
   this%cx1 = cx1
   this%cx2 = cx2
   this%cx3 = cx3
   this%cy1 = cy1
   this%cy2 = cy2
   this%cy3 = cy3
   this%gamma = gamma
   this%np = np
   this%quiet = quiet
   this%evol = evol

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist3d_000
!
subroutine dist3d_000(this,part3d,npp,noff,ndp)

   implicit none

   class(fdist3d_000), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: part3d
   integer(kind=LG), intent(inout) :: npp
   integer, intent(in), dimension(2) :: noff, ndp
! local data1
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
   real, dimension(:,:), pointer :: pt => null()
   integer :: npx, npy, npz
   real :: vtx, vty, vtz, vdx, vdy, vdz,dr,dz
   real :: sigx, sigy, sigz, x0, y0, z0, rmax, zmin, zmax
   real, dimension(3) :: cx, cy
   real, dimension(4) :: edges
   integer(kind=LG) :: nps, npmax
   logical :: lquiet = .false.
   integer :: idimp, ierr
   character(len=18), save :: sname = 'dist3d_000'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ierr = 0; nps = 1
   npx = this%npx; npy = this%npy; npz = this%npz
   pt => part3d
   vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
   vdx = 0.0; vdy = 0.0; vdz = this%gamma
   sigx = this%sigx; sigy = this%sigy; sigz = this%sigz
   x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
   cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
   lquiet = this%quiet
   idimp = size(part3d,1); npmax = size(part3d,2)
   dr = this%dx; dz= this%dz
   rmax = this%rmax
   zmax = this%zmax
   zmin = this%zmin
   if (noff(1) == 0) then
      edges(1) = 0.0
      edges(2) = edges(1) + ndp(1) * dr
   else
      edges(1) = noff(1) * dr
      edges(2) = edges(1) + ndp(1) * dr
   end if
   edges(3) = noff(2)*dz + zmin
   edges(4) = edges(3) + ndp(2)*dz

   call beam_dist000(pt,this%qm,edges,npp,this%dx,this%dz,nps,vtx,vty,vtz,vdx,vdy,&
   &vdz,npx,npy,npz,rmax,zmin,zmax,idimp,npmax,sigx,sigy,sigz,&
   &x0,y0,z0,cx,cy,lquiet,ierr)

   if (ierr /= 0) then
      write (erstr,*) 'beam_dist000 error'
      call write_err(cls_name//sname//erstr)
   endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine dist3d_000
!
subroutine deposit_fdist3d_000(this,q)

   implicit none

   class(fdist3d_000), intent(inout) :: this
   class(field), intent(inout) :: q
! local data
   class(ufield), dimension(:), pointer :: q_re, q_im
   real, dimension(:,:,:), pointer :: q0
   real :: r, z, dr, dz
   integer :: i, j, noff1, noff2, n1p, n2p
   real :: np, sigx, sigz, sigx2, sigz2
   real :: bcz

   q_im => null()
   q_re => q%get_rf_re()
   q0 => q_re(0)%get_f2()
   noff1 = q_re(0)%get_noff(1)
   noff2 = q_re(0)%get_noff(2)
   n1p = q_re(0)%get_ndp(1)
   n2p = q_re(0)%get_ndp(2)
   dr = this%dx
   dz = this%dz
   if (abs(this%qm) < 1.0e-6 ) then
      np = 0.0
   else
      np = this%np*this%qm/abs(this%qm)
   end if
   sigx = this%sigx
   sigz = this%sigz
   bcz = this%bcz
   sigx2 = 0.5/sigx**2
   sigz2 = 0.5/sigz**2

   do i = 1, n1p
      r = (i + noff1 - 1.0) * dr
      do j = 1, n2p+1
         z = (j + noff2) * dz - bcz
         q0(1,i,j) = np*exp(-r**2*sigx2-z**2*sigz2)
      end do
   end do

end subroutine deposit_fdist3d_000

subroutine init_fdist3d_001(this,input,i)

   implicit none

   class(fdist3d_001), intent(inout) :: this
   type(input_json), intent(inout), pointer :: input
   integer, intent(in) :: i
! local data
   integer :: npf,npr,npth,npz,npmax
   real :: qm,sigx,sigy,sigz,bcx,bcy,bcz,sigvx,sigvy,sigvz
   real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
   logical :: quiet, evol
   real :: min, max, cwp, n0
   real :: alx, alz, dr, dz
   integer :: nr, nz, num_modes
   character(len=20) :: sn,s1
   character(len=18), save :: sname = 'init_fdist3d_001:'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%p => input%pp

   write (sn,'(I3.3)') i
   s1 = 'beam('//trim(sn)//')'
   call input%get('simulation.n0',n0)
   call input%get('simulation.grid(1)',nr)
   call input%get('simulation.grid(2)',nz)
   call input%get('simulation.max_mode',num_modes)
   cwp=5.32150254*1e9/sqrt(n0)
   call input%get('simulation.box.r(1)',min)
   call input%get('simulation.box.r(2)',max)
   call input%get(trim(s1)//'.center(1)',bcx)
   alx = (max-min)
   dr=alx/real(nr)
   this%rmax = max
   call input%get(trim(s1)//'.center(2)',bcy)
   call input%get('simulation.box.z(1)',min)
   call input%get('simulation.box.z(2)',max)
   call input%get(trim(s1)//'.center(3)',bcz)
   bcz = bcz -min
   alz = (max-min)
   dz=alz/real(nz)
   this%z0 = min
   this%zmin = 0.0
   this%zmax = max-min
   call input%get(trim(s1)//'.profile',npf)
   call input%get(trim(s1)//'.np(1)',npr)
   call input%get(trim(s1)//'.np(2)',npth)
   call input%get(trim(s1)//'.np(3)',npz)
   call input%get(trim(s1)//'.q',qm)
   call input%get(trim(s1)//'.sigma(1)',sigx)
   call input%get(trim(s1)//'.sigma(2)',sigy)
   call input%get(trim(s1)//'.sigma(3)',sigz)
   call input%get(trim(s1)//'.sigma_v(1)',sigvx)
   call input%get(trim(s1)//'.sigma_v(2)',sigvy)
   call input%get(trim(s1)//'.sigma_v(3)',sigvz)
   call input%get(trim(s1)//'.centroid_x(1)',cx1)
   call input%get(trim(s1)//'.centroid_x(2)',cx2)
   call input%get(trim(s1)//'.centroid_x(3)',cx3)
   call input%get(trim(s1)//'.centroid_y(1)',cy1)
   call input%get(trim(s1)//'.centroid_y(2)',cy2)
   call input%get(trim(s1)//'.centroid_y(3)',cy3)
   call input%get(trim(s1)//'.quiet_start',quiet)
   call input%get(trim(s1)//'.gamma',gamma)
   call input%get(trim(s1)//'.peak_density',np)
   call input%get(trim(s1)//'.npmax',npmax)
   call input%get(trim(s1)//'.evolution',evol)
   this%npf   = npf
   this%dx    = dr
   this%dz    = dz
   this%npr   = npr
   this%npth  = npth
   this%npz   = npz
   this%npmax = npmax
   this%qm    = qm/abs(qm)*2.0*np / (dr*dr*dz*npth*npr*npz)
   this%bcx   = bcx
   this%bcy   = bcy
   this%bcz   = bcz
   this%sigx  = sigx
   this%sigy  = sigy
   this%sigz  = sigz
   this%sigvx = sigvx
   this%sigvy = sigvy
   this%sigvz = sigvz
   this%cx1   = cx1
   this%cx2   = cx2
   this%cx3   = cx3
   this%cy1   = cy1
   this%cy2   = cy2
   this%cy3   = cy3
   this%gamma = gamma
   this%np    = np
   this%quiet = quiet
   this%evol  = evol

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist3d_001
!
subroutine dist3d_001(this,part3d,npp,noff,ndp)

   implicit none

   class(fdist3d_001), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: part3d
   integer(kind=LG), intent(inout) :: npp
   integer, intent(in), dimension(2) :: noff, ndp
! local data1
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
   real, dimension(:,:), pointer :: pt => null()
   integer :: npr, npth, npz
   real :: vtx, vty, vtz, vdx, vdy, vdz,dr,dz
   real :: sigx, sigy, sigz, x0, y0, z0, rmax, zmin, zmax
   real, dimension(3) :: cx, cy
   real, dimension(4) :: edges
   integer(kind=LG) :: nps, npmax
   logical :: lquiet = .false.
   integer :: idimp, ierr
   character(len=18), save :: sname = 'dist3d_001'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ierr = 0; nps = 1
   npr = this%npr; npth = this%npth; npz = this%npz
   pt => part3d
   vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
   vdx = 0.0; vdy = 0.0; vdz = this%gamma
   sigx = this%sigx; sigy = this%sigy; sigz = this%sigz
   x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
   cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
   lquiet = this%quiet
   idimp = size(part3d,1); npmax = size(part3d,2)
   dr = this%dx; dz= this%dz
   rmax = this%rmax
   zmax = this%zmax
   zmin = this%zmin
   if (noff(1) == 0) then
      edges(1) = 0.0
      edges(2) = edges(1) + ndp(1) * dr
   else
      edges(1) = noff(1) * dr
      edges(2) = edges(1) + ndp(1)*dr
   end if
   edges(3) = noff(2)*dz+zmin
   edges(4) = edges(3) + ndp(2)*dz

   call beam_dist001(pt,this%qm,edges,npp,this%dx,this%dz,nps,vtx,vty,vtz,vdx,vdy,&
   &vdz,npr,npth,npz,rmax,zmin,zmax,idimp,npmax,sigx,sigy,sigz,&
   &x0,y0,z0,cx,cy,lquiet,ierr)

   if (ierr /= 0) then
      write (erstr,*) 'beam_dist001 error'
      call write_err(cls_name//sname//erstr)
   endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine dist3d_001
!
subroutine deposit_fdist3d_001(this,q)

   implicit none

   class(fdist3d_001), intent(inout) :: this
   class(field), intent(inout) :: q
! local data
   class(ufield), dimension(:), pointer :: q_re, q_im
   real, dimension(:,:,:), pointer :: q0
   real :: r, z, dr, dz
   integer :: i, j, noff1, noff2, n1p, n2p
   real :: np, sigx, sigz, sigx2, sigz2
   real :: bcz

   q_im => null()
   q_re => q%get_rf_re()
   q0 => q_re(0)%get_f2()
   noff1 = q_re(0)%get_noff(1)
   noff2 = q_re(0)%get_noff(2)
   n1p = q_re(0)%get_ndp(1)
   n2p = q_re(0)%get_ndp(2)
   dr = this%dx
   dz = this%dz
   if (abs(this%qm) < 1.0e-6 ) then
      np = 0.0
   else
      np = this%np*this%qm/abs(this%qm)
   end if
   sigx = this%sigx
   sigz = this%sigz
   bcz = this%bcz
   sigx2 = 0.5/sigx**2
   sigz2 = 0.5/sigz**2

   do i = 1, n1p
      r = (i + noff1 - 1.0) * dr
      do j = 1, n2p+1
         z = (j + noff2) * dz - bcz
         q0(1,i,j) = np*exp(-r**2*sigx2-z**2*sigz2)
      end do
   end do

end subroutine deposit_fdist3d_001
!
subroutine init_fdist3d_002(this,input,i)

   implicit none

   class(fdist3d_002), intent(inout) :: this
   type(input_json), intent(inout), pointer :: input
   integer, intent(in) :: i
! local data
   integer :: npf,npx,npy,npz,npmax
   real :: qm,sigx,sigy,bcx,bcy,bcz,sigvx,sigvy,sigvz
   real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
   logical :: quiet, evol
   real :: min, max, cwp, n0
   real :: alx, alz, dr, dz, sumz
   integer :: nr, nz, num_modes, ii
   character(len=20) :: sn,s1
   character(len=18), save :: sname = 'init_fdist3d_002:'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%p => input%pp

   write (sn,'(I3.3)') i
   s1 = 'beam('//trim(sn)//')'
   call input%get('simulation.n0',n0)
   call input%get('simulation.grid(1)',nr)
   call input%get('simulation.grid(2)',nz)
   call input%get('simulation.max_mode',num_modes)
   cwp=5.32150254*1e9/sqrt(n0)
   call input%get('simulation.box.r(1)',min)
   call input%get('simulation.box.r(2)',max)
   call input%get(trim(s1)//'.center(1)',bcx)
   alx = (max-min)
   dr=alx/real(nr)
   this%rmax = max
   call input%get(trim(s1)//'.center(2)',bcy)
   call input%get('simulation.box.z(1)',min)
   call input%get('simulation.box.z(2)',max)
   call input%get(trim(s1)//'.center(3)',bcz)
   bcz = bcz -min
   alz = (max-min)
   dz=alz/real(nz)
   this%z0 = min
   this%zmin = 0.0
   this%zmax = max-min
   call input%get(trim(s1)//'.profile',npf)
   call input%get(trim(s1)//'.np(1)',npx)
   call input%get(trim(s1)//'.np(2)',npy)
   call input%get(trim(s1)//'.np(3)',npz)
   call input%get(trim(s1)//'.q',qm)
   call input%get(trim(s1)//'.sigma(1)',sigx)
   call input%get(trim(s1)//'.sigma(2)',sigy)
   call input%get(trim(s1)//'.sigma_v(1)',sigvx)
   call input%get(trim(s1)//'.sigma_v(2)',sigvy)
   call input%get(trim(s1)//'.sigma_v(3)',sigvz)
   call input%get(trim(s1)//'.centroid_x(1)',cx1)
   call input%get(trim(s1)//'.centroid_x(2)',cx2)
   call input%get(trim(s1)//'.centroid_x(3)',cx3)
   call input%get(trim(s1)//'.centroid_y(1)',cy1)
   call input%get(trim(s1)//'.centroid_y(2)',cy2)
   call input%get(trim(s1)//'.centroid_y(3)',cy3)
   call input%get(trim(s1)//'.quiet_start',quiet)
   call input%get(trim(s1)//'.gamma',gamma)
   call input%get(trim(s1)//'.peak_density',np)
   call input%get(trim(s1)//'.npmax',npmax)
   call input%get(trim(s1)//'.piecewise_fz',this%fz)
   call input%get(trim(s1)//'.piecewise_z',this%z)
   call input%get(trim(s1)//'.evolution',evol)

   sumz = 0.0
   do ii = 2, size(this%z)
      if (this%z(ii) <= this%z(ii-1)) then
         write (erstr,*) 'Piecewise_z is not monotonically increasing'
         call write_err(cls_name//sname//erstr)
      endif
      sumz = sumz + (this%fz(ii) + this%fz(ii-1)) * (this%z(ii) - this%z(ii-1)) * 0.5
   end do
   this%z = this%z/dz
   this%npf = npf
   this%dx = dr
   this%dz = dz
   this%npx = npx
   this%npy = npy
   this%npz = npz
   this%npmax = npmax
   qm = qm/abs(qm)*np*sigx*sigy*sumz
   qm = qm/dr/dz/dr
   qm = qm/npx
   qm = qm/npy
   qm = qm/npz
   this%qm = qm
   this%bcx = bcx
   this%bcy = bcy
   this%bcz = bcz
   this%sigx = sigx
   this%sigy = sigy
   this%sigvx = sigvx
   this%sigvy = sigvy
   this%sigvz = sigvz
   this%cx1 = cx1
   this%cx2 = cx2
   this%cx3 = cx3
   this%cy1 = cy1
   this%cy2 = cy2
   this%cy3 = cy3
   this%gamma = gamma
   this%np = np
   this%quiet = quiet
   this%evol = evol

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist3d_002

subroutine dist3d_002(this,part3d,npp,noff,ndp)

   implicit none

   class(fdist3d_002), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: part3d
   integer(kind=LG), intent(inout) :: npp
   integer, intent(in), dimension(2) :: noff, ndp
! local data1
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
   real, dimension(:,:), pointer :: pt => null()
   integer :: npx, npy, npz, nzf, nz, i, j
   real :: vtx, vty, vtz, vdx, vdy, vdz,dr,dz
   real :: sigx, sigy, x0, y0, z0, rmax, zmin, zmax
   real, dimension(3) :: cx, cy
   real, dimension(4) :: edges
   integer(kind=LG) :: nps, npmax
   logical :: lquiet = .false.
   integer :: idimp, ierr
   real, dimension(:), allocatable :: zf
   character(len=18), save :: sname = 'dist3d_002'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ierr = 0; nps = 1
   npx = this%npx; npy = this%npy; npz = this%npz
   nzf = size(this%z)
   nz = (this%zmax - this%zmin) / this%dz
   if ((this%z(1)>=nz) .or. (this%z(nzf)<=0)) then
      npp = 0
      return
   endif
   pt => part3d
   vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
   vdx = 0.0; vdy = 0.0; vdz = this%gamma
   sigx = this%sigx; sigy = this%sigy
   x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
   cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
   lquiet = this%quiet
   idimp = size(part3d,1); npmax = size(part3d,2)
   dr = this%dx; dz= this%dz
   rmax = this%rmax
   zmax = this%zmax
   zmin = this%zmin
   if (noff(1) == 0) then
      edges(1) = 0.0
      edges(2) = edges(1) + ndp(1) * dr
   else
      edges(1) = noff(1) * dr
      edges(2) = edges(1) + ndp(1) * dr
   end if
   edges(3) = noff(2)*dz + zmin
   edges(4) = edges(3) + ndp(2)*dz

   allocate(zf(nz))
   zf = 0.0
   do i = 1, nz
      do j = 2, nzf
         if ((i >= this%z(j-1)) .and. (i < this%z(j))) then
            zf(i) = this%fz(j) + (this%fz(j-1)-this%fz(j)) /&
               (this%z(j-1)-this%z(j)) * (real(i)-this%z(j))
            exit
         endif
      enddo
   enddo

   call beam_dist002(pt,this%qm,edges,npp,this%dx,this%dz,nps,vtx,vty,vtz,vdx,vdy,&
   &vdz,npx,npy,npz,rmax,zmin,zmax,idimp,npmax,sigx,sigy,&
   &x0,y0,z0,cx,cy,zf,lquiet,ierr)

   if (ierr /= 0) then
      write (erstr,*) 'beam_dist002 error'
      call write_err(cls_name//sname//erstr)
   endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine dist3d_002

subroutine deposit_fdist3d_002(this,q)

   implicit none

   class(fdist3d_002), intent(inout) :: this
   class(field), intent(inout) :: q
! local data
   class(ufield), dimension(:), pointer :: q_re, q_im
   real, dimension(:,:,:), pointer :: q0
   real :: r, z, dr, dz, fz
   integer :: i, j, k, noff1, noff2, n1p, n2p, nzf
   real :: np, sigx, sigx2

   q_im => null()
   q_re => q%get_rf_re()
   q0 => q_re(0)%get_f2()
   noff1 = q_re(0)%get_noff(1)
   noff2 = q_re(0)%get_noff(2)
   n1p = q_re(0)%get_ndp(1)
   n2p = q_re(0)%get_ndp(2)
   nzf = size(this%z)
   dr = this%dx
   dz = this%dz
   if (abs(this%qm) < 1.0e-6 ) then
      np = 0.0
   else
      np = this%np*this%qm/abs(this%qm)
   end if
   sigx = this%sigx
   sigx2 = 0.5/sigx**2

   do i = 1, n1p
      r = (i + noff1 - 1.0) * dr
      do j = 1, n2p+1
         z = (j + noff2) + this%z0 / dz
         if (z < this%z(1) .or. z > this%z(nzf)) then
            fz = 0.0
         else
            do k = 2, nzf
               if ((z >= this%z(k-1)) .and. (z < this%z(k))) then
                  fz = this%fz(k) + (this%fz(k-1)-this%fz(k)) /&
                     (this%z(k-1)-this%z(k)) * (z-this%z(k))
                  exit
               endif
            enddo
         endif
         q0(1,i,j) = np*exp(-r**2*sigx2) * fz
      end do
   end do

end subroutine deposit_fdist3d_002

end module fdist3d_class