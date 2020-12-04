module fdist2d_class

use parallel_module
use options_class
use ufield_class
use input_class
use param
use sysutil_module


implicit none

private

public :: fdist2d, fdist2d_wrap, fdist2d_000, fdist2d_012

type, abstract :: fdist2d

   ! private

   ! ndprof = profile type
   integer :: npf, npmax
   real :: dex

   integer :: noff, nr, nrp

   integer :: ppc1, ppc2, nmode, num_theta
   real :: qm, den
   character(len=:), allocatable :: long_prof
   real, dimension(:), allocatable :: s, fs

   contains
   generic :: new => init_fdist2d
   generic :: del => end_fdist2d
   generic :: dist => dist2d
   procedure(ab_init_fdist2d), deferred, private :: init_fdist2d
   procedure, private :: end_fdist2d
   procedure(ab_dist2d), deferred, private :: dist2d
   procedure :: getnpf, getnpmax, getdex, get_density

end type fdist2d

abstract interface

subroutine ab_dist2d( this, x, p, gamma, q, psi, npp, s )
   import fdist2d
   import ufield
   import LG
   implicit none
   class(fdist2d), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:), pointer, intent(inout) :: gamma, q, psi
   integer(kind=LG), intent(inout) :: npp
   ! class(ufield), intent(in), pointer :: ud
   real, intent(in) :: s
end subroutine ab_dist2d

subroutine ab_init_fdist2d(this,input,opts,i,sect)
   import fdist2d
   import input_json
   import options
   implicit none
   class(fdist2d), intent(inout) :: this
   type(input_json), intent(inout) :: input
   type(options), intent(in) :: opts
   integer, intent(in) :: i
   character(len=*), intent(in) :: sect
end subroutine ab_init_fdist2d

end interface

type fdist2d_wrap
   class(fdist2d), allocatable :: p
end type fdist2d_wrap

type, extends(fdist2d) :: fdist2d_000
! Transeversely uniform profile with uniform or piecewise longitudinal profile
   private
! xppc, yppc = particle per cell in x and y directions
   ! integer :: ppc1, ppc2, nmode
   ! real :: qm, den
   ! character(len=:), allocatable :: long_prof
   ! real, dimension(:), allocatable :: s, fs

   contains
   procedure, private :: init_fdist2d => init_fdist2d_000
   procedure, private :: dist2d => dist2d_000

end type fdist2d_000

type, extends(fdist2d) :: fdist2d_012
! hollow channel with f(r) profile
   private
! xppc, yppc = particle per cell in x and y directions
   ! integer :: ppc1, ppc2, nmode
   ! real :: qm, den
   ! real :: cx, cy
   real, dimension(:), allocatable :: r, fr
   ! character(len=:), allocatable :: long_prof
   ! real, dimension(:), allocatable :: s, fs

   contains
   procedure, private :: init_fdist2d => init_fdist2d_012
   procedure, private :: dist2d => dist2d_012

end type fdist2d_012

character(len=20), parameter :: cls_name = "fdist2d"
integer, parameter :: cls_level = 2
character(len=128) :: erstr

contains

function getnpf(this)

   implicit none

   class(fdist2d), intent(in) :: this
   integer :: getnpf

   getnpf = this%npf

end function getnpf

function getnpmax(this)

   implicit none

   class(fdist2d), intent(in) :: this
   integer :: getnpmax

   getnpmax = this%npmax

end function getnpmax

function getdex(this)

   implicit none

   class(fdist2d), intent(in) :: this
   real :: getdex

   getdex = this%dex

end function getdex

function get_density( this, s )

   implicit none

   class(fdist2d), intent(in) :: this
   real, intent(in) :: s

   real :: get_density
   integer :: i, prof_len

   get_density = 1.0
   if ( trim(this%long_prof) == 'piecewise' ) then

      prof_len = size( this%fs )
      if ( s <= this%s(1) .or. s > this%s(prof_len) ) then
         call write_err( 'The longitudinal position is out of the bound!' )
      endif

      do i = 2, prof_len

         if ( this%s(i) < this%s(i-1) ) then
            call write_err( 's is not monotonically increasing!' )
         endif

         if ( s <= this%s(i) ) then
            get_density = this%fs(i-1) + ( this%fs(i) - this%fs(i-1) ) &
               / ( this%s(i) - this%s(i-1) ) * ( s - this%s(i-1) )
            exit
         endif

      enddo
   endif

end function get_density

subroutine end_fdist2d(this)

   implicit none

   class(fdist2d), intent(inout) :: this

   character(len=18), save :: sname = 'end_fdist2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_fdist2d

subroutine init_fdist2d_000(this,input,opts,i,sect)

   implicit none

   class(fdist2d_000), intent(inout) :: this
   type(input_json), intent(inout) :: input
   type(options), intent(in) :: opts
   integer, intent(in) :: i
   character(len=*), intent(in) :: sect

   ! local data
   integer :: npf,ppc1,ppc2,n1,nmode,num_theta, xtra
   integer(kind=LG) :: npmax
   real :: qm, den, lr, ur
   character(len=20) :: sn,s1
   character(len=18), save :: sname = 'init_fdist2d_000'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%noff = opts%get_noff(1)
   this%nr   = opts%get_nd(1)
   this%nrp  = opts%get_ndp(1)

   ! write (sn,'(I3.3)') i
   ! s1 = 'species('//trim(sn)//')'
   s1 = trim(sect) // '(' // num2str(i) // ')'
   call input%get('simulation.grid(1)',n1)
   call input%get('simulation.box.r(1)',lr)
   call input%get('simulation.box.r(2)',ur)
   call input%get('simulation.max_mode',nmode)
   call input%get(trim(s1)//'.profile',npf)
   call input%get(trim(s1)//'.ppc(1)',ppc1)
   call input%get(trim(s1)//'.ppc(2)',ppc2)
   call input%get(trim(s1)//'.num_theta', num_theta)
   call input%get(trim(s1)//'.q',qm)
   call input%get(trim(s1)//'.density',den)
   call input%get(trim(s1)//'.longitudinal_profile',this%long_prof)
   if (trim(this%long_prof) == 'piecewise') then
      call input%get(trim(s1)//'.piecewise_density',this%fs)
      call input%get(trim(s1)//'.piecewise_s',this%s)
   end if
   this%dex = (ur - lr)/real(n1)
   this%npf = npf
   this%nmode = nmode
   this%ppc1 = ppc1
   this%ppc2 = ppc2
   this%num_theta = num_theta
   this%qm = qm
   this%den = den
   xtra = 10
   this%npmax = this%nrp * ppc1 * ppc2 * num_theta * xtra
   call write_dbg(cls_name, sname, cls_level, 'ends')
end subroutine init_fdist2d_000

subroutine dist2d_000( this, x, p, gamma, q, psi, npp, s )
   implicit none
   class(fdist2d_000), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:), pointer, intent(inout) :: gamma, q, psi
   integer(kind=LG), intent(inout) :: npp
   real, intent(in) :: s
   ! local data
   character(len=18), save :: sname = 'dist2d_000'
   integer(kind=LG) :: nps
   integer :: nr, nrp, n_theta, ppc1, ppc2, i1, i2, noff, i, j
   real :: qm, den_temp
   integer :: prof_l
   real :: r1, theta, dr, dtheta

   call write_dbg(cls_name, sname, cls_level, 'starts')

   nr   = this%nr
   nrp  = this%nrp
   noff = this%noff
   n_theta = this%num_theta
   
   ppc1 = this%ppc1; ppc2 = this%ppc2
   dr = this%dex
   dtheta = 2.0 * pi / n_theta
   den_temp = 1.0
   if (noff+nrp == nr .and. nrp > 2) nrp = nrp - 2
   if (trim(this%long_prof) == 'piecewise') then
      prof_l = size(this%fs)
      if (s<this%s(1) .or. s>this%s(prof_l)) then
         write (erstr,*) 'The s is out of the bound!'
         call write_err(trim(erstr))
         return
      end if
      do i = 2, prof_l
         if (this%s(i) < this%s(i-1)) then
            write (erstr,*) 's is not monotonically increasing!'
            call write_err(trim(erstr))
            return
         end if
         if (s<=this%s(i)) then
            den_temp = this%fs(i-1) + (this%fs(i)-this%fs(i-1))/&
            &(this%s(i)-this%s(i-1))*(s-this%s(i-1))
            exit
         end if
      end do
   end if
   qm = den_temp*this%den*this%qm/abs(this%qm)/real(ppc1)/real(ppc2)/real(n_theta)
   
   ! initialize the particle positions
   nps = 0
   do j = 1, n_theta
      do i = 1, nrp
         do i1 = 1, ppc1
            r1 = (i1 - 0.5)/ppc1 + i - 1 + noff
            do i2 = 1, ppc2
               nps = nps + 1
               theta = ( (i2 - 0.5) / ppc2 + j - 1.5 ) * dtheta
               x(1,nps) = r1 * dr * cos(theta)
               x(2,nps) = r1 * dr * sin(theta)
               p(1,nps) = 0.0
               p(2,nps) = 0.0
               p(3,nps) = 0.0
               gamma(nps) = 1.0
               psi(nps) = 1.0
               q(nps) = qm*r1
            enddo
         enddo
      enddo
   enddo

   npp = nps

   call write_dbg(cls_name, sname, cls_level, 'ends')
end subroutine dist2d_000
!
subroutine init_fdist2d_012(this,input,opts,i,sect)

   implicit none

   class(fdist2d_012), intent(inout) :: this
   type(input_json), intent(inout) :: input
   type(options), intent(in) :: opts
   integer, intent(in) :: i
   character(len=*), intent(in) :: sect
! local data
   integer :: npf,ppc1,ppc2,n1,nmode, num_theta, xtra
   integer(kind=LG) :: npmax
   real :: qm, den, lr, ur
   character(len=20) :: sn,s1
   character(len=18), save :: sname = 'init_fdist2d_012:'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%noff = opts%get_noff(1)
   this%nr   = opts%get_nd(1)
   this%nrp  = opts%get_ndp(1)

   ! write (sn,'(I3.3)') i
   ! s1 = 'species('//trim(sn)//')'
   s1 = trim(sect) // '(' // num2str(i) // ')'
   call input%get('simulation.grid(1)',n1)
   call input%get('simulation.box.r(1)',lr)
   call input%get('simulation.box.r(2)',ur)
   call input%get('simulation.max_mode',nmode)
   call input%get(trim(s1)//'.profile',npf)
   call input%get(trim(s1)//'.ppc(1)',ppc1)
   call input%get(trim(s1)//'.ppc(2)',ppc2)
   call input%get(trim(s1)//'.q',qm)
   call input%get(trim(s1)//'.density',den)
   call input%get(trim(s1)//'.longitudinal_profile',this%long_prof)
   if (trim(this%long_prof) == 'piecewise') then
      call input%get(trim(s1)//'.piecewise_density',this%fs)
      call input%get(trim(s1)//'.piecewise_s',this%s)
   end if
   call input%get(trim(s1)//'.piecewise_radial_density',this%fr)
   call input%get(trim(s1)//'.piecewise_r',this%r)

   this%dex = (ur - lr)/real(n1)
   this%npf = npf
   this%nmode = nmode
   this%ppc1 = ppc1
   this%ppc2 = ppc2
   this%num_theta = num_theta
   this%qm = qm
   this%den = den
   xtra = 10
   this%npmax = this%nrp * ppc1 * ppc2 * num_theta * xtra
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist2d_012
!
subroutine dist2d_012(this,x,p,gamma,q,psi,npp,s)
   implicit none
   class(fdist2d_012), intent(inout) :: this
   real, dimension(:,:), pointer, intent(inout) :: x, p
   real, dimension(:), pointer, intent(inout) :: gamma, q, psi
   integer(kind=LG), intent(inout) :: npp
   real, intent(in) :: s
! local data
   character(len=18), save :: sname = 'dist2d_012:'
   integer(kind=LG) :: nps
   integer :: nr, nrp, ppc1, ppc2, i1, i2, noff, i, j, n_theta
   real :: qm, den_temp
   integer :: prof_l, ii
   real :: r1, theta, dtheta, dr, rr

   call write_dbg(cls_name, sname, cls_level, 'starts')

   nr   = this%nr
   nrp  = this%nrp
   noff = this%noff
   n_theta = this%num_theta

   ppc1 = this%ppc1
   ppc2 = this%ppc2
   dr = this%dex
   dtheta = 2.0 * pi / n_theta
   den_temp = 1.0

   if (trim(this%long_prof) == 'piecewise') then
      prof_l = size(this%fs)
      if (prof_l /= size(this%s)) then
         write (erstr,*) 'The piecewise_density and s array have different sizes!'
         call write_err(trim(erstr))
         return
      end if
      if (s<this%s(1) .or. s>this%s(prof_l)) then
         write (erstr,*) 'The s is out of the bound!'
         call write_err(trim(erstr))
         return
      end if
      do i = 2, prof_l
         if (this%s(i) < this%s(i-1)) then
            write (erstr,*) 's is not monotonically increasing!'
            call write_err(trim(erstr))
            return
         end if
         if (s<=this%s(i)) then
            den_temp = this%fs(i-1) + (this%fs(i)-this%fs(i-1))/&
            &(this%s(i)-this%s(i-1))*(s-this%s(i-1))
            exit
         end if
      end do
   end if
   qm = den_temp*this%den*this%qm/abs(this%qm)/real(ppc1)/real(ppc2)/real(n_theta)
   prof_l = size(this%fr)
   if (prof_l /= size(this%r)) then
      write (erstr,*) 'The piecewise_radial_density and r array have different sizes!'
      call write_err(trim(erstr))
      return
   end if

   ! initialize the particle positions
   nps = 0
   do j = 1, n_theta
      do i = 1, nrp
         do i1 = 1, ppc1
            rr = (i1 - 0.5)/ppc1 + i - 1 + noff
            r1 = rr*dr
            if (r1<this%r(1) .or. r1>this%r(prof_l)) then
               cycle
            end if
            do ii = 2, prof_l
               if (this%r(ii) <= this%r(ii-1)) then
                  write (erstr,*) 'r is not monotonically increasing!'
                  call write_err(trim(erstr))
                  return
               end if
               if (r1<=this%r(ii)) then
                  den_temp = this%fr(ii-1) + (this%fr(ii)-this%fr(ii-1))/&
                  &(this%r(ii)-this%r(ii-1))*(r1-this%r(ii-1))
                  exit
               end if
            end do
            do i2 = 1, ppc2
               nps = nps + 1
               theta = ( (i2 - 0.5) / ppc2 + j - 1.5 ) * dtheta
               x(1,nps) = r1 * cos(theta)
               x(2,nps) = r1 * sin(theta)
               p(1,nps) = 0.0
               p(2,nps) = 0.0
               p(3,nps) = 0.0
               gamma(nps) = 1.0
               psi(nps) = 1.0
               q(nps) = qm * den_temp * rr
            enddo
         enddo
      enddo
   enddo

   npp = nps

   call write_dbg(cls_name, sname, cls_level, 'ends')
end subroutine dist2d_012
!
end module fdist2d_class