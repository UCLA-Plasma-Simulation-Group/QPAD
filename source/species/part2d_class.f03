! part2d class for QPAD

module part2d_class

use param
use sys
use parallel_pipe_class
use field_class
use ufield_class
use fdist2d_class
use hdf5io_class
use part2d_lib
         
implicit none

private

public :: part2d

type part2d
   
   private
!
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! xdim = dimension of the particle coordinates
! nbmax = size of buffer for passing particles between processors
! npp = number of particles in current partition
! npmax = maximum number of particles in each partition
! part(:,:) = initial particle coordinates
!         
   real :: qbm, dt, dex
   integer :: xdim
   integer(kind=LG) :: npmax, nbmax, npp = 0
   real, dimension(:,:), pointer :: part => null()
   class(parallel_pipe), pointer :: pp => null()
   
   contains
   
   generic :: new => init_part2d
   generic :: renew => renew_part2d
   generic :: del => end_part2d
   generic :: qdp => qdeposit
   generic :: amjdp => amjdeposit
   generic :: push => partpush
   generic :: pmv => pmove
   generic :: extpsi => extractpsi
   generic :: psend => pipesend_part2d
   generic :: precv => piperecv_part2d
   generic :: wr => writehdf5_part2d
   procedure, private :: init_part2d, renew_part2d
   procedure, private :: end_part2d
   procedure, private :: qdeposit
   procedure, private :: amjdeposit
   procedure, private :: partpush
   procedure, private :: pmove
   procedure, private :: extractpsi
   procedure, private :: pipesend_part2d
   procedure, private :: piperecv_part2d, writehdf5_part2d
                     
end type 

save

character(len=20), parameter :: cls_name = "part2d"
integer, parameter :: cls_level = 2
character(len=128) :: erstr
real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
integer(kind=LG), dimension(:), allocatable :: ihole

contains
!
subroutine init_part2d(this,pp,pf,fd,qbm,dt,xdim,s)

   implicit none
   
   class(part2d), intent(inout) :: this
   class(parallel_pipe), intent(in), pointer :: pp
   class(fdist2d), intent(inout) :: pf
   class(field), intent(in) :: fd
   real, intent(in) :: qbm, dt, s
   integer, intent(in) :: xdim

! local data
   character(len=18), save :: sname = 'init_part2d'
   integer :: xtras, noff, nxyp, nx, npmax, nbmax
   class(ufield), pointer :: ud

   call write_dbg(cls_name, sname, cls_level, 'starts')
   this%qbm = qbm
   this%dt = dt
   this%xdim = xdim
   npmax = pf%getnpmax()
   this%dex = pf%getdex()
   this%npmax = npmax
   nbmax = max(int(0.01*npmax),100)
   this%nbmax = nbmax
   ud => fd%get_rf_re(0)
   this%pp => pp

   allocate(this%part(xdim,npmax))
   call pf%dist(this%part,this%npp,ud,s)
   if (.not. allocated(sbufl)) then
      allocate(sbufl(xdim,nbmax),sbufr(xdim,nbmax))
      allocate(rbufl(xdim,nbmax),rbufr(xdim,nbmax))
      allocate(ihole(nbmax*2))
   end if

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_part2d
!
subroutine end_part2d(this)
    
   implicit none
   
   class(part2d), intent(inout) :: this
   character(len=18), save :: sname = 'end_part2d'
   
   call write_dbg(cls_name, sname, cls_level, 'starts')
   deallocate(this%part)
   call write_dbg(cls_name, sname, cls_level, 'ends')
      
end subroutine end_part2d
!
subroutine renew_part2d(this,pf,fd,s)

   implicit none
   
   class(part2d), intent(inout) :: this
   class(fdist2d), intent(inout) :: pf
   class(field), intent(in) :: fd
   real, intent(in) :: s
! local data
   character(len=18), save :: sname = 'renew_part2d'
   integer :: noff, prof
   class(ufield), pointer :: ud
         
   call write_dbg(cls_name, sname, cls_level, 'starts')
   
   ud => fd%get_rf_re(0)
   call pf%dist(this%part,this%npp,ud,s)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine renew_part2d
!      
subroutine qdeposit(this,q)
! deposit the charge density (rho - Jz)     
      
   implicit none
   
   class(part2d), intent(in) :: this
   class(field), intent(in) :: q
! local data
   character(len=18), save :: sname = 'qdeposit'
   class(ufield), dimension(:), pointer :: q_re => null(), q_im => null()
            
   call write_dbg(cls_name, sname, cls_level, 'starts')
   
   q_re => q%get_rf_re()
   q_im => q%get_rf_im()
   
   call part2d_qdeposit(this%part,this%npp,this%dex,q_re,q_im,q%get_num_modes())         
   
   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine qdeposit
!      
subroutine amjdeposit(this,ef,bf,cu,amu,dcu)
! deposit the current, acceleration and momentum flux      
      
   implicit none
   
   class(part2d), intent(inout) :: this
   class(field), intent(in) :: cu, amu, dcu
   class(field), intent(in) :: ef, bf
! local data
   character(len=18), save :: sname = 'amjdeposit'
   class(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
   class(ufield), dimension(:), pointer :: bf_re => null(), bf_im => null()
   class(ufield), dimension(:), pointer :: cu_re => null(), cu_im => null()
   class(ufield), dimension(:), pointer :: dcu_re => null(), dcu_im => null()
   class(ufield), dimension(:), pointer :: amu_re => null(), amu_im => null()

   call write_dbg(cls_name, sname, cls_level, 'starts')
   
   ef_re => ef%get_rf_re()
   ef_im => ef%get_rf_im()
   bf_re => bf%get_rf_re()
   bf_im => bf%get_rf_im()
   cu_re => cu%get_rf_re()
   cu_im => cu%get_rf_im()
   dcu_re => dcu%get_rf_re()
   dcu_im => dcu%get_rf_im()
   amu_re => amu%get_rf_re()
   amu_im => amu%get_rf_im()
   
   call part2d_amjdeposit(this%part,this%npp,this%dex,this%dt,this%qbm,&
   &ef_re,ef_im,bf_re,bf_im,cu_re,cu_im,dcu_re,dcu_im,amu_re,amu_im,&
   &ef%get_num_modes())

   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine amjdeposit
!      
subroutine partpush(this,ef,bf)
      
   implicit none
   
   class(part2d), intent(inout) :: this
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
   
   call part2d_push(this%part,this%npp,this%dex,this%xdim,this%dt,this%qbm,&
   &ef_re,ef_im,bf_re,bf_im,ef%get_num_modes())

   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine partpush
!      
subroutine pmove(this,fd)
      
   implicit none
         
   class(part2d), intent(inout) :: this
   class(field), intent(in) :: fd
! local data   
   character(len=18), save :: sname = 'pmove'
   class(ufield), pointer :: ud

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ud => fd%get_rf_re(0)
   call part2d_pmove(this%part,this%pp,this%npp,this%dex,this%xdim,this%npmax,&
   &this%nbmax,ud,sbufl,sbufr,rbufl,rbufr,ihole)
   
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine pmove
!      
subroutine extractpsi(this,psi)
      
   implicit none
   
   class(part2d), intent(inout) :: this
   class(field), intent(in) :: psi
! local data
   character(len=18), save :: sname = 'extractpsi'
   class(ufield), dimension(:), pointer :: psi_re => null(), psi_im => null()
   
   call write_dbg(cls_name, sname, cls_level, 'starts')

   psi_re => psi%get_rf_re()
   psi_im => psi%get_rf_im()
   call part2d_extractpsi(this%part,this%npp,this%dex,this%qbm,psi_re,psi_im,psi%get_num_modes())

   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine extractpsi
!
subroutine pipesend_part2d(this)
      
   implicit none
         
   class(part2d), intent(inout) :: this

   return
!          integer, intent(in) :: tag
!          integer, intent(inout) :: id
! ! local data
!          character(len=18), save :: sname = 'pipesend_part2d:'
!          integer :: des, ierr
         
         
!          call this%err%werrfl2(class//sname//' started')
         
!          des = this%p%getidproc()+this%p%getlnvp()
         
!          if (des >= this%p%getnvp()) then
!             id = MPI_REQUEST_NULL         
!             call this%err%werrfl2(class//sname//' ended')
!             return
!          endif
         
!          call this%pcb()
                  
!          call MPI_ISEND(this%part,this%npp*this%xdim,this%p%getmreal(),&
!          &des,tag,this%p%getlworld(),id,ierr)

! ! check for errors
!          if (ierr /= 0) then
!             write (erstr,*) 'MPI_ISEND failed'
!             call this%err%equit(class//sname//erstr); return
!          endif

!          call this%err%werrfl2(class//sname//' ended')
         
end subroutine pipesend_part2d
!      
subroutine piperecv_part2d(this)
      
   implicit none
         
   class(part2d), intent(inout) :: this

   return
!          class(ufield2d), pointer, intent(in) :: fd
!          integer, intent(in) :: tag
! ! local data
!          character(len=18), save :: sname = 'piperecv_part2d:'
!          integer, dimension(10) :: istat
!          integer :: nps, id, des, ierr
         
         
!          call this%err%werrfl2(class//sname//' started')

!          des = this%p%getidproc()-this%p%getlnvp()
         
!          if (des < 0) then
!             call this%err%werrfl2(class//sname//' ended')
!             return
!          endif

!          call MPI_IRECV(this%part,this%npmax*this%xdim,this%p%getmreal(),&
!          &des,tag,this%p%getlworld(),id,ierr)

!          call MPI_WAIT(id,istat,ierr)
         
!          call MPI_GET_COUNT(istat,this%p%getmreal(),nps,ierr)

!          this%npp = nps/this%xdim
         
!          call this%pcp(fd)
         
! ! check for errors
!          if (ierr /= 0) then
!             write (erstr,*) 'MPI failed'
!             call this%err%equit(class//sname//erstr); return
!          endif
         
!          call this%err%werrfl2(class//sname//' ended')
         
end subroutine piperecv_part2d
!      
subroutine writehdf5_part2d(this,file)

   implicit none
         
   class(part2d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'writehdf5_part2d'
   integer :: ierr

   call write_dbg(cls_name, sname, cls_level, 'starts')

   call pwpart(this%pp,file,this%part,this%npp,1,ierr)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_part2d
!      
end module part2d_class