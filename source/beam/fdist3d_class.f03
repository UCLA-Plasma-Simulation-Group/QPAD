! fdist3d class for QPAD

module fdist3d_class

use parallel_pipe_class
use ufield_class
use input_class
use part3d_lib
use param
use system

   
implicit none

private

public :: fdist3d, fdist3d_000

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
   procedure(ab_init_fdist3d), deferred, private :: init_fdist3d
   procedure, private :: end_fdist3d
   procedure(ab_dist3d), deferred, private :: dist3d
   procedure :: getnpf, getnpmax, getevol, getdx, getdz, getz0
            
end type

abstract interface
!
subroutine ab_dist3d(this,part3d,npp,ud)
   import fdist3d
   import ufield
   implicit none
   class(fdist3d), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: part3d
   integer, intent(inout) :: npp
   class(ufield), intent(in), pointer :: ud
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
end interface
!
type, extends(fdist3d) :: fdist3d_000
! Tri Gaussian profile (the same particle charge)
   private

   integer :: npx, npy, npz
   real :: qm, sigx, sigy, sigz
   real :: bcx, bcy, bcz, sigvx, sigvy, sigvz
   real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
   logical :: quiet

   contains

   procedure, private :: init_fdist3d => init_fdist3d_000
   procedure, private :: dist3d => dist3d_000

end type fdist3d_000
!
character(len=10), save :: class = 'fdist3d'
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

end function getdex
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

end function getdex
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
   call input%get(trim(s1)//'.center(2)',bcy)
   call input%get('simulation.box.z(1)',min)
   call input%get('simulation.box.z(2)',max)
   call input%get(trim(s1)//'.center(3)',bcz)
   bcz = bcz -min
   alz = (max-min) 
   dz=alz/real(nz)
   this%z0 = min
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
   qm = qm*np*(2*3.1415926535897932)**1.5*sigx*sigy*sigz
   qm = qm*nr
   qm = qm*nz
   qm = qm/(npx*alx*alz) 
   qm = qm/npy
   qm = qm/npz
   this%qm = qm
   this%bcx = bcx/dr
   this%bcy = bcy/dr
   this%bcz = bcz/dz
   this%sigx = sigx/dr
   this%sigy = sigy/dr
   this%sigz = sigz/dz
   this%sigvx = sigvx
   this%sigvy = sigvy
   this%sigvz = sigvz
   this%cx1 = cx1*dz*dz/dr
   this%cx2 = cx2*dz/dr
   this%cx3 = cx3/dr
   this%cy1 = cy1*dz*dz/dr
   this%cy2 = cy2*dz/dr
   this%cy3 = cy3/dr
   this%gamma = gamma
   this%np = np
   this%quiet = quiet
   this%evol = evol
   this%dex = dr

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist3d_000
!
subroutine dist3d_000(this,part3d,npp,ud)

   implicit none
   
   class(fdist3d_000), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: part3d
   integer, intent(inout) :: npp
   class(ufield), intent(in), pointer :: ud
! local data1
! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
   real, dimension(:,:), pointer :: pt => null()
   integer :: npx, npy, npz, n1, n2, ipbc
   real :: vtx, vty, vtz, vdx, vdy, vdz
   real :: sigx, sigy, sigz, x0, y0, z0
   real, dimension(3) :: cx, cy
   real, dimension(4) :: edges
   integer, dimension(2) :: noff
   integer :: nps=1
   logical :: lquiet = .false.
   integer :: idimp, npmax, ierr = 0
   character(len=18), save :: sname = 'dist3d_000'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   
   npx = this%npx; npy = this%npy; npz = this%npz
   n1 = ud%get_nd(1); n2 = ud%get_nd(2)
   pt => part3d
   vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
   vdx = 0.0; vdy = 0.0; vdz = this%gamma
   sigx = this%sigx; sigy = this%sigy; sigz = this%sigz
   x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
   cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
   lquiet = this%quiet
   idimp = size(part3d,1); npmax = size(part3d,2)
   noff = ud%get_noff()
   if (noff(1) == 0) then
      edges(1) = noff(1)
      edges(2) = edges(1) + ud%get_ndp(1)
      edges(2) = edges(2) - 0.5
   else
      edges(1) = noff(1)
      edges(1) = edges(1) - 0.5
      edges(2) = edges(1) + ud%get_ndp(1)
   end if
   edges(3) = noff(2)
   edges(4) = edges(3) + ud%get_ndp(2)         
   
   call beam_dist000(pt,this%qm,edges,npp,nps,vtx,vty,vtz,vdx,vdy,&
   &vdz,npx,npy,npz,n1,n2,idimp,npmax,sigx,sigy,sigz,&
   &x0,y0,z0,cx,cy,lquiet,ierr)

   if (ierr /= 0) then
      write (erstr,*) 'beam_dist000 error'
      call write_err(cls_name//sname//erstr)
   endif
   
   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine dist3d_000
!
end module fdist3d_class