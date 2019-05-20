! part3d class for QPAD

module part3d_class

use param
use sys
use parallel_pipe_class
use field_class
use ufield_class
use fdist3d_class
use hdf5io_class
use part3d_lib
use mpi
         
implicit none

private

public :: part3d

type part3d

   private

!
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! xdim = dimension of the particle coordinates
! nbmax = size of buffer for passing particles between processors
! npp = number of particles in current partition
! npmax = maximum number of particles in each partition
! part(:,:) = initial particle coordinates
!         
   class(parallel_pipe), pointer, public :: pp => null()
   real :: qbm, dt, dx, dz
   real :: z0
   integer(kind=LG) :: npmax, nbmax, npp
   integer :: xdim
   real, dimension(:,:), pointer :: part => null(), pbuff => null()
   
   contains
   
   generic :: new => init_part3d
   generic :: del => end_part3d
   generic :: push => partpush
   generic :: pmv => pmove
   generic :: qdp => qdeposit  
   generic :: wr => writehdf5_part3d
   generic :: wrst => writerst_part3d
   generic :: rrst => readrst_part3d
   procedure, private :: init_part3d
   procedure, private :: end_part3d
   procedure, private :: partpush
   procedure, private :: pmove         
   procedure, private :: qdeposit, writehdf5_part3d
   procedure, private :: writerst_part3d, readrst_part3d
   procedure :: getnpp
            
end type 

save      

character(len=20), parameter :: cls_name = "part3d"
integer, parameter :: cls_level = 2
character(len=128) :: erstr

! buffer data for particle managers
real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
integer(kind=LG), dimension(:), allocatable :: ihole

contains
!
subroutine init_part3d(this,pp,pf,fd,qbm,dt,xdim)

   implicit none
   
   class(part3d), intent(inout) :: this
   class(parallel_pipe), intent(in), pointer :: pp
   class(fdist3d), intent(inout) :: pf
   class(field), pointer, intent(in) :: fd
   real, intent(in) :: qbm, dt
   integer, intent(in) :: xdim
! local data
   character(len=18), save :: sname = 'init_part3d'
   integer :: noff, nxyp, nx, prof, npmax, nbmax
   class(ufield), pointer :: ud
            
   call write_dbg(cls_name, sname, cls_level, 'starts')
   this%pp => pp
   this%qbm = qbm
   this%dt = dt
   this%xdim = xdim
   npmax = pf%getnpmax()
   this%dx = pf%getdx()
   this%dz = pf%getdz()
   this%npmax = npmax
   nbmax = int(0.01*this%npmax)
   this%nbmax = nbmax
   this%npp = 0
   prof = pf%getnpf()
   allocate(this%part(xdim,npmax),this%pbuff(xdim,nbmax))
   ud => fd%get_rf_re(0)
   call pf%dist(this%part,this%npp,ud)
   this%z0 = pf%getz0()
   if (.not. allocated(sbufl)) then
      allocate(sbufl(xdim,nbmax),sbufr(xdim,nbmax))
      allocate(rbufl(xdim,nbmax),rbufr(xdim,nbmax))
      allocate(ihole(nbmax*2))
   end if   
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_part3d
!
subroutine end_part3d(this)
    
   implicit none
   
   class(part3d), intent(inout) :: this
   character(len=18), save :: sname = 'end_part3d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   deallocate(this%part,this%pbuff)
   call write_dbg(cls_name, sname, cls_level, 'ends')
   return
   
end subroutine end_part3d
!      
subroutine qdeposit(this,q)

   implicit none
   
   class(part3d), intent(inout) :: this
   class(field), intent(in) :: q
! local data
   character(len=18), save :: sname = 'qdeposit'
   class(ufield), dimension(:), pointer :: q_re => null(), q_im => null()

   call write_dbg(cls_name, sname, cls_level, 'starts')
   
   q_re => q%get_rf_re()
   q_im => q%get_rf_im()
   
   call part3d_qdeposit(this%part,this%npp,this%dx,this%dz,q_re,q_im,q%get_num_modes())         
   
   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine qdeposit
!
subroutine partpush(this,ef,bf)
      
   implicit none
   
   class(part3d), intent(inout) :: this
   class(field), intent(in) :: ef, bf
! local data
   character(len=18), save :: sname = 'partpush'
   class(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
   class(ufield), dimension(:), pointer :: bf_re => null(), bf_im => null()
   
   call write_dbg(cls_name, sname, cls_level, 'starts')
   
   ef_re => ef%get_rf_re()
   ef_im => ef%get_rf_im()
   bf_re => bf%get_rf_re()
   bf_im => bf%get_rf_im()
   
   call part3d_push(this%part,this%npp,this%dx,this%dz,this%xdim,this%dt,this%qbm,&
   &ef_re,ef_im,bf_re,bf_im,ef%get_num_modes())

   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine partpush
!
subroutine pmove(this,fd,rtag,stag,sid)

   implicit none
   
   class(part3d), intent(inout) :: this
   class(field), intent(in) :: fd
   integer, intent(in) :: rtag, stag
   integer, intent(inout) :: sid
! local data
   character(len=18), save :: sname = 'pmove:'
   class(ufield), pointer :: ud
   integer, dimension(9) :: info


   call write_dbg(cls_name, sname, cls_level, 'starts')

   ud => fd%get_rf_re(0)

   call part3d_pmove(this%part,this%pp,ud,this%npp,this%dx,this%dz,sbufr,sbufl,&
   &rbufr,rbufl,ihole,this%pbuff,this%xdim,this%npmax,this%nbmax,rtag,stag,sid,info)

   if (info(1) /= 0) then
      write (erstr,*) 'part3d_pmove error'
      call write_err(cls_name//sname//erstr)
   endif
   
   if (this%pp%getstageid() == this%pp%getnstage() - 1) sid = MPI_REQUEST_NULL
   
   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine pmove
!      
function getnpp(this)

   implicit none

   class(part3d), intent(in) :: this
   integer(kind=LG) :: getnpp
   
   getnpp = this%npp

end function getnpp
!
subroutine writehdf5_part3d(this,file,dspl,rtag,stag,id)

   implicit none
   
   class(part3d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
   integer, intent(in) :: dspl, rtag, stag
   integer, intent(inout) :: id
! local data
   character(len=18), save :: sname = 'writehdf5_part3d'
   integer :: ierr = 0

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call pwpart_pipe(this%pp,file,this%part,this%npp,dspl,&
   &this%z0,rtag,stag,id,ierr)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_part3d
!            
subroutine writerst_part3d(this,file)

   implicit none
   
   class(part3d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'writerst_part3d'
   integer :: ierr = 0

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call wpart(this%pp,file,this%part,this%npp,1,ierr)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writerst_part3d
!            
subroutine readrst_part3d(this,file)

   implicit none
   
   class(part3d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'readrst_part3d'
   integer :: ierr = 0

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call rpart(this%pp,file,this%part,this%npp,ierr)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine readrst_part3d
!            
end module part3d_class