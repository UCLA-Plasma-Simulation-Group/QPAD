! fdist3d class for QPAD

! *TODO* fdist3d class should be higher than part3d class
module fdist3d_class

use parallel_pipe_class
use field_class
use ufield_class
use input_class
use fdist3d_lib
use param
use sysutil
use hdf5
use hdf5io_class
use iso_c_binding

implicit none

private

public :: fdist3d, fdist3d_wrap
public :: fdist3d_000, fdist3d_001, fdist3d_002, fdist3d_100

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
subroutine ab_dist3d(this,x,p,q,npp,noff,ndp,s)
   import fdist3d
   import ufield
   import LG
   implicit none
   class(fdist3d), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:,:), pointer, intent(inout), optional :: s
   real, dimension(:), pointer, intent(inout) :: q
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

type, extends(fdist3d) :: fdist3d_100
! External particle imported from OSIRIS
   private

   ! integer :: npx, npy, npz
   real :: rmax, zmin, zmax
   real :: os_n0, os_fraction
   real, dimension(3) :: os_ctr, bctr
   character(len=:), allocatable :: filename
   logical :: quiet

   contains

   procedure, private :: init_fdist3d => init_fdist3d_100
   procedure, private :: deposit_fdist3d => deposit_fdist3d_100
   procedure, private :: dist3d => dist3d_100

end type fdist3d_100

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
subroutine dist3d_000(this,x,p,q,npp,noff,ndp,s)

   implicit none

   class(fdist3d_000), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:,:), pointer, intent(inout), optional :: s
   real, dimension(:), pointer, intent(inout) :: q
   integer(kind=LG), intent(inout) :: npp
   integer, intent(in), dimension(2) :: noff, ndp
! local data1
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
   ! real, dimension(:,:), pointer :: pt => null()
   integer :: npx, npy, npz, i
   real :: vtx, vty, vtz, vdx, vdy, vdz,dr,dz
   real :: sigx, sigy, sigz, x0, y0, z0, rmax, zmin, zmax
   real, dimension(3) :: cx, cy
   real, dimension(4) :: edges
   integer(kind=LG) :: nps, npmax
   logical :: lquiet = .false.
   integer :: ierr
   character(len=18), save :: sname = 'dist3d_000'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ierr = 0; nps = 1
   npx = this%npx; npy = this%npy; npz = this%npz
   ! pt => part3d
   vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
   vdx = 0.0; vdy = 0.0; vdz = this%gamma
   sigx = this%sigx; sigy = this%sigy; sigz = this%sigz
   x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
   cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
   lquiet = this%quiet
   npmax = size(x,2)
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

   call beam_dist000(x,p,q,this%qm,edges,npp,this%dx,this%dz,nps,vtx,vty,vtz,vdx,vdy,&
   &vdz,npx,npy,npz,rmax,zmin,zmax,npmax,sigx,sigy,sigz,&
   &x0,y0,z0,cx,cy,lquiet,ierr)

   if (present(s)) then
      do i = 1, npp
         s(1,i) = 0.0
         s(2,i) = 0.0
         s(3,i) = 1.0
      enddo
   endif

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
subroutine dist3d_001(this,x,p,q,npp,noff,ndp,s)

   implicit none

   class(fdist3d_001), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:,:), pointer, intent(inout), optional :: s
   real, dimension(:), pointer, intent(inout) :: q
   integer(kind=LG), intent(inout) :: npp
   integer, intent(in), dimension(2) :: noff, ndp
! local data1
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
   ! real, dimension(:,:), pointer :: pt => null()
   integer :: npr, npth, npz, i
   real :: vtx, vty, vtz, vdx, vdy, vdz,dr,dz
   real :: sigx, sigy, sigz, x0, y0, z0, rmax, zmin, zmax
   real, dimension(3) :: cx, cy
   real, dimension(4) :: edges
   integer(kind=LG) :: nps, npmax
   logical :: lquiet = .false.
   integer :: ierr
   character(len=18), save :: sname = 'dist3d_001'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ierr = 0; nps = 1
   npr = this%npr; npth = this%npth; npz = this%npz
   ! pt => part3d
   vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
   vdx = 0.0; vdy = 0.0; vdz = this%gamma
   sigx = this%sigx; sigy = this%sigy; sigz = this%sigz
   x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
   cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
   lquiet = this%quiet
   npmax = size(x,2)
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

   call beam_dist001(x,p,q,this%qm,edges,npp,this%dx,this%dz,nps,vtx,vty,vtz,vdx,vdy,&
   &vdz,npr,npth,npz,rmax,zmin,zmax,npmax,sigx,sigy,sigz,&
   &x0,y0,z0,cx,cy,lquiet,ierr)

   if (present(s)) then
      do i = 1, npp
         s(1,i) = 0.0
         s(2,i) = 0.0
         s(3,i) = 1.0
      enddo
   endif

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

subroutine dist3d_002(this,x,p,q,npp,noff,ndp,s)

   implicit none

   class(fdist3d_002), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:,:), pointer, intent(inout), optional :: s
   real, dimension(:), pointer, intent(inout) :: q
   integer(kind=LG), intent(inout) :: npp
   integer, intent(in), dimension(2) :: noff, ndp
! local data1
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
   ! real, dimension(:,:), pointer :: pt => null()
   integer :: npx, npy, npz, nzf, nz, i, j
   real :: vtx, vty, vtz, vdx, vdy, vdz,dr,dz
   real :: sigx, sigy, x0, y0, z0, rmax, zmin, zmax
   real, dimension(3) :: cx, cy
   real, dimension(4) :: edges
   integer(kind=LG) :: nps, npmax
   logical :: lquiet = .false.
   integer :: ierr
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
   ! pt => part3d
   vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
   vdx = 0.0; vdy = 0.0; vdz = this%gamma
   sigx = this%sigx; sigy = this%sigy
   x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
   cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
   lquiet = this%quiet
   npmax = size(x,2)
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

   call beam_dist002(x,p,q,this%qm,edges,npp,this%dx,this%dz,nps,vtx,vty,vtz,vdx,vdy,&
   &vdz,npx,npy,npz,rmax,zmin,zmax,npmax,sigx,sigy,&
   &x0,y0,z0,cx,cy,zf,lquiet,ierr)

   if (present(s)) then
      do i = 1, npp
         s(1,i) = 0.0
         s(2,i) = 0.0
         s(3,i) = 1.0
      enddo
   endif

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

subroutine init_fdist3d_100(this,input,i)

   implicit none

   class(fdist3d_100), intent(inout) :: this
   type(input_json), intent(inout), pointer :: input
   integer, intent(in) :: i

   ! local data
   real :: min, max
   integer :: nr, nz
   character(len=20) :: s1
   character(len=32), save :: sname = 'init_fdist3d_100:'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%p => input%pp

   s1 = 'beam('//num2str(i)//')'

   ! regular params
   call input%get(trim(s1)//'.profile', this%npf)
   call input%get(trim(s1)//'.npmax', this%npmax)
   call input%get(trim(s1)//'.evolution', this%evol)
   ! call input%get(trim(s1)//'.quiet_start',quiet)
   ! call input%get('simulation.n0', n0)
   call input%get('simulation.grid(1)', nr)
   call input%get('simulation.box.r(1)', min)
   call input%get('simulation.box.r(2)', max)
   this%dx = (max - min) / real(nr)
   this%rmax = max
   call input%get('simulation.grid(2)', nz)
   call input%get('simulation.box.z(1)', min)
   call input%get('simulation.box.z(2)', max)
   this%dz = (max - min) / real(nz)
   this%z0 = min
   this%zmin = 0.0
   this%zmax = max - min

   ! params for external particle import
   call input%get(trim(s1)//'.os_center(1)', this%os_ctr(1))
   call input%get(trim(s1)//'.os_center(2)', this%os_ctr(2))
   call input%get(trim(s1)//'.os_center(3)', this%os_ctr(3))
   call input%get(trim(s1)//'.center(1)', this%bctr(1))
   call input%get(trim(s1)//'.center(2)', this%bctr(2))
   call input%get(trim(s1)//'.center(3)', this%bctr(3))
   ! may not need dx in OSIRIS
   ! call input%get(trim(s1)//'.os_dx(1)', this%os_dx(1))
   ! call input%get(trim(s1)//'.os_dx(2)', this%os_dx(2))
   ! call input%get(trim(s1)//'.os_dx(3)', this%os_dx(3))
   call input%get(trim(s1)//'.os_n0', this%os_n0)
   call input%get(trim(s1)//'.os_fraction', this%os_fraction)
   call input%get(trim(s1)//'.filename', this%filename)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist3d_100
!
subroutine dist3d_100(this,x,p,q,npp,noff,ndp,s)

   use iso_c_binding

   implicit none

   class(fdist3d_100), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:,:), pointer, intent(inout), optional :: s
   real, dimension(:), pointer, intent(inout) :: q
   integer(kind=LG), intent(inout) :: npp
   integer, intent(in), dimension(2) :: noff, ndp

   ! local data
   real :: rbuf, xconv, qconv
   real, dimension(2) :: redge, zedge
   integer :: i, pp, ptrcur, ierr
   character(len=18), save :: sname = 'dist3d_100'

   integer(hsize_t), dimension(2) :: dims, maxdims
   integer(hsize_t), dimension(2) :: chunk_size, offset
   integer(hid_t) :: file_id, grp_id, treal, dtype_id
   integer(hid_t), dimension(p_x_dim) :: xdset_id
   integer(hid_t), dimension(p_p_dim) :: pdset_id
   integer(hid_t), dimension(p_s_dim) :: sdset_id
   integer(hid_t) :: qdset_id, fspace_id, mspace_id
   real, dimension(:,:), pointer :: xbuf, pbuf, sbuf
   real, dimension(:), pointer :: qbuf
   logical :: has_spin = .false.
   integer, parameter :: real_kind_8 = kind(1.0d0)

   call write_dbg(cls_name, sname, cls_level, 'starts')

   if ( present(s) ) has_spin = .true.

   ! conversion factor from OSIRIS to QPAD
   xconv = 1.0 / sqrt( this%os_n0 )
   qconv = this%os_n0 / ( pi * this%os_fraction )

   allocate( xbuf(p_cache_size, p_x_dim) )
   allocate( pbuf(p_cache_size, p_p_dim) )
   allocate( qbuf(p_cache_size) )
   if (has_spin) allocate( sbuf(p_cache_size, p_s_dim) )

   treal = detect_precision()

   if (noff(1) == 0) then
      redge(1) = 0.0
   else
      redge(1) = noff(1) * this%dx
   endif
   redge(2) = redge(1) + ndp(1) * this%dx
   zedge(1) = noff(2) * this%dz + this%zmin
   zedge(2) = zedge(1) + ndp(2) * this%dz

   call h5open_f( ierr )
   call h5fopen_f( this%filename, H5F_ACC_RDONLY_F, file_id, ierr )
   call h5gopen_f( file_id, '/', grp_id, ierr )

   ! get metadata of dataset to be read
   call h5dopen_f( grp_id, 'x1', xdset_id(1), ierr )
   call h5dget_space_f( xdset_id(1), fspace_id, ierr )
   call h5dget_type_f( xdset_id(1), dtype_id, ierr )
   call h5sget_simple_extent_dims_f( fspace_id, dims, maxdims, ierr )
   call h5sclose_f( fspace_id, ierr )
   call h5dclose_f( xdset_id(1), ierr )

   ! create memory dataspace
   call h5screate_simple_f( 2, (/int(p_cache_size, hsize_t), 1_hsize_t/), &
      mspace_id, ierr )

   ! open all the datasets to be read
   call h5dopen_f( grp_id, 'x2', xdset_id(1), ierr )
   call h5dopen_f( grp_id, 'x3', xdset_id(2), ierr )
   call h5dopen_f( grp_id, 'x1', xdset_id(3), ierr )
   call h5dopen_f( grp_id, 'p2', pdset_id(1), ierr )
   call h5dopen_f( grp_id, 'p3', pdset_id(2), ierr )
   call h5dopen_f( grp_id, 'p1', pdset_id(3), ierr )
   call h5dopen_f( grp_id, 'q', qdset_id, ierr )
   if ( has_spin ) then
      call h5dopen_f( grp_id, 's2', sdset_id(1), ierr )
      call h5dopen_f( grp_id, 's3', sdset_id(2), ierr )
      call h5dopen_f( grp_id, 's1', sdset_id(3), ierr )
   endif

   pp = 0
   offset = 0
   chunk_size(1) = 0
   chunk_size(2) = 1
   do ptrcur = 1, dims(1), p_cache_size

      ! check if last copy of table and set np
      if( ptrcur + p_cache_size > dims(1) ) then
        chunk_size(1) = dims(1) - ptrcur + 1
        call h5sset_extent_simple_f( mspace_id, 2, (/int(chunk_size(1), hsize_t), 1_hsize_t/), &
         (/int(chunk_size(1), hsize_t), 1_hsize_t/), ierr )
      else
        chunk_size(1) = p_cache_size
      endif

      if ( pp + chunk_size(1) > this%npmax ) then
         call write_err('Insufficient memory allocated. "npmax" is too small.')
      endif

      ! read chunk from datasets x1, x2, x3
      do i = 1, p_x_dim
         call h5dget_space_f( xdset_id(i), fspace_id, ierr )
         call h5sselect_hyperslab_f( fspace_id, H5S_SELECT_SET_F, offset, &
            chunk_size, ierr) 
         call h5dread_f( xdset_id(i), h5kind_to_type(real_kind_8, H5_REAL_KIND), &
            xbuf(1,i), chunk_size, ierr, mspace_id, fspace_id )
         call h5sclose_f( fspace_id, ierr )
      enddo

      ! read chunk from datasets p1, p2, p3
      do i = 1, p_p_dim
         call h5dget_space_f( pdset_id(i), fspace_id, ierr )
         call h5sselect_hyperslab_f( fspace_id, H5S_SELECT_SET_F, offset, &
            chunk_size, ierr) 
         call h5dread_f( pdset_id(i), h5kind_to_type(real_kind_8, H5_REAL_KIND), &
            pbuf(1,i), chunk_size, ierr, mspace_id, fspace_id )
         call h5sclose_f( fspace_id, ierr )
      enddo

      ! read chunk from dataset q
      call h5dget_space_f( qdset_id, fspace_id, ierr )
      call h5sselect_hyperslab_f( fspace_id, H5S_SELECT_SET_F, offset, &
         chunk_size, ierr) 
      call h5dread_f( qdset_id, h5kind_to_type(real_kind_8, H5_REAL_KIND), &
         qbuf, chunk_size, ierr, mspace_id, fspace_id )
      call h5sclose_f( fspace_id, ierr )

      ! read chunk from datasets s1, s2, s3
      if ( has_spin ) then
      do i = 1, p_s_dim
         call h5dget_space_f( sdset_id(i), fspace_id, ierr )
         call h5sselect_hyperslab_f( fspace_id, H5S_SELECT_SET_F, offset, &
            chunk_size, ierr) 
         call h5dread_f( sdset_id(i), h5kind_to_type(real_kind_8, H5_REAL_KIND), &
            sbuf(1,i), chunk_size, ierr, mspace_id, fspace_id )
         call h5sclose_f( fspace_id, ierr )
      enddo
      endif

      ! convert and store particle data
      do i = 1, chunk_size(1)

         ! centralize particle position
         xbuf(i,1) = ( xbuf(i,1) - this%os_ctr(1) ) * xconv
         xbuf(i,2) = ( xbuf(i,2) - this%os_ctr(2) ) * xconv
         xbuf(i,3) = ( xbuf(i,3) - this%os_ctr(3) ) * xconv

         ! shift in QPAD space
         xbuf(i,1) = this%bctr(1) + xbuf(i,1)
         xbuf(i,2) = this%bctr(2) + xbuf(i,2)
         xbuf(i,3) = this%bctr(3) - xbuf(i,3)

         rbuf = sqrt( xbuf(i,1)**2 + xbuf(i,2)**2 )
         if ( rbuf >= redge(1) .and. rbuf < redge(2) .and. &
              xbuf(i,3) >= zedge(1) .and. xbuf(i,3) < zedge(2) ) then
            pp = pp + 1
            x(1,pp) = xbuf(i,1)
            x(2,pp) = xbuf(i,2)
            x(3,pp) = xbuf(i,3)
            p(1,pp) = pbuf(i,1)
            p(2,pp) = pbuf(i,2)
            p(3,pp) = pbuf(i,3)
            q(pp) = qbuf(i) * qconv
            if (has_spin) then
               s(1,pp) = sbuf(i,1)
               s(2,pp) = sbuf(i,2)
               s(3,pp) = sbuf(i,3)
            endif
         endif

      enddo

      offset(1) = offset(1) + chunk_size(1)

   enddo

   npp = pp

   call h5sclose_f( mspace_id, ierr )
   call h5dclose_f( xdset_id(1), ierr )
   call h5dclose_f( xdset_id(2), ierr )
   call h5dclose_f( xdset_id(3), ierr )
   call h5dclose_f( pdset_id(1), ierr )
   call h5dclose_f( pdset_id(2), ierr )
   call h5dclose_f( pdset_id(3), ierr )
   call h5dclose_f( qdset_id, ierr )
   if (has_spin) then
      call h5dclose_f( sdset_id(1), ierr )
      call h5dclose_f( sdset_id(2), ierr )
      call h5dclose_f( sdset_id(3), ierr )
   endif

   call h5gclose_f( grp_id, ierr )
   call h5fclose_f( file_id, ierr )
   call h5close_f( ierr )

   deallocate( xbuf, pbuf, qbuf )
   if (has_spin) deallocate( sbuf )

   if (ierr /= 0) then
      call write_err(cls_name//sname//'beam_dist100 error')
   endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine dist3d_100

!
subroutine deposit_fdist3d_100(this,q)

   implicit none

   class(fdist3d_100), intent(inout) :: this
   class(field), intent(inout) :: q

   call write_err( 'Free-stream deposition is not supported for external particles.' )   

end subroutine deposit_fdist3d_100

end module fdist3d_class