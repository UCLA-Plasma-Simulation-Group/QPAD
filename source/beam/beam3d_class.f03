! beam3d class for QPAD

module beam3d_class

use parallel_pipe_class
use param
use system
use fdist3d_class
use grid_class
use field_src_class
use field_class
use field_e_class
use field_b_class
use part3d_class
use hdf5io_class
use mpi

implicit none

private

public :: beam3d

type beam3d

   private

   class(parallel_pipe), pointer, public :: pp => null()
   class(part3d), pointer :: pd => null()
   class(field_rho), allocatable :: q
   class(fdist3d), pointer :: pf => null()
   logical :: evol
   contains
   
   generic :: new => init_beam3d
   generic :: del => end_beam3d
   generic :: push => push_beam3d
   generic :: qdp => qdeposit_beam3d
   generic :: wr => writehdf5_beam3d
   generic :: wrq => writeq_beam3d
   generic :: wrst => writerst_beam3d       
   generic :: rrst => readrst_beam3d       
   procedure, private :: init_beam3d
   procedure, private :: end_beam3d
   procedure, private :: push_beam3d
   procedure, private :: qdeposit_beam3d, writehdf5_beam3d
   procedure, private :: writerst_beam3d, readrst_beam3d
   procedure, private :: writeq_beam3d
end type

save      

character(len=10) :: cls_name = 'beam3d'
integer, parameter :: cls_level = 2
character(len=128) :: erstr

contains
!
subroutine init_beam3d(this,pp,fd,gd,part_shape,pf,qbm,dt,xdim,smooth_type,smooth_order)

   implicit none
   
   class(beam3d), intent(inout) :: this
   class(grid), intent(in), pointer :: gd
   class(field), intent(in), pointer :: fd
   class(parallel_pipe), intent(in), pointer :: pp
   class(fdist3d), intent(inout), target :: pf
   real, intent(in) :: qbm, dt
   integer, intent(in) :: xdim, part_shape
   integer, intent(in), optional :: smooth_type, smooth_order
! local data
   character(len=18), save :: sname = 'init_beam3d'
   integer :: id, num_modes, ierr
   integer, dimension(10) :: istat
   real :: dr, dxi
            
   call write_dbg(cls_name, sname, cls_level, 'starts')
   this%pp => pp
   this%pf => pf
   dr = fd%get_dr()
   dxi = fd%get_dxi()
   num_modes = fd%get_num_modes()

   allocate(this%pd,this%q)
   this%evol = pf%getevol()
   if ( present(smooth_type) .and. present(smooth_order) ) then
      call this%q%new(pp,gd,dr,dxi,num_modes,part_shape,smooth_type,smooth_order)
   else
      call this%q%new(pp,gd,dr,dxi,num_modes,part_shape)
   endif
   call this%pd%new(pp,pf,fd,qbm,dt,xdim)
   call this%pd%pmv(this%q,1,1,id)
   call MPI_WAIT(id,istat,ierr)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_beam3d
!
subroutine end_beam3d(this)
    
   implicit none
   
   class(beam3d), intent(inout) :: this
   character(len=18), save :: sname = 'end_beam3d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%pd%del()
   call this%q%del()
   call write_dbg(cls_name, sname, cls_level, 'ends')
            
end subroutine end_beam3d
!      
subroutine qdeposit_beam3d(this,q)
! deposit the charge density      

   implicit none
   
   class(beam3d), intent(inout) :: this
   class(field_rho), intent(inout) :: q
! local data
   character(len=18), save :: sname = 'qdeposit_beam3d'
            
   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'deposit 3D particles' )

   if (.not. this%evol) then
      call this%q%as(0.0)
      call this%pf%dp(this%q)
      call this%q%copy_gc_f2()
      call add_f2( this%q, q )
   else
      call this%q%as(0.0)
      call this%pd%qdp(this%q)
      call this%q%acopy_gc_f2()
      call this%q%copy_gc_f2()
      call add_f2( this%q, q )
   end if
   
   call stop_tprof( 'deposit 3D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine qdeposit_beam3d
!
! subroutine qdeposit_beam3d(this,id1,id2,id3,tag1,tag2)
! ! deposit the charge density      

!    implicit none
   
!    class(beam3d), intent(inout) :: this
!    integer, intent(inout) :: id1, id2, id3, tag1, tag2
! ! local data
!    character(len=18), save :: sname = 'qdeposit_beam3d'
!    integer, dimension(10) :: istat
!    integer :: ierr
            
!    call this%err%werrfl2(class//sname//' started')
!    call this%q%as(0.0)
!    call MPI_WAIT(id1,istat,ierr)
!    call MPI_WAIT(id3,istat,ierr)
!    call this%pd%qdp(this%q%getrs())
!    call this%q%ag(tag1,tag1,id1)
!    call this%q%pcg(tag2,tag2,id2,id3)    
           
!    call this%err%werrfl2(class//sname//' ended')
   
! end subroutine qdeposit_beam3d
!     
subroutine push_beam3d(this,ef,bf,rtag,stag,sid)

   implicit none
   
   class(beam3d), intent(inout) :: this
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf
   integer, intent(in) :: rtag, stag
   integer, intent(inout) :: sid         
! local data
   character(len=18), save :: sname = 'partpush'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'push 3D particles' )

   if (.not. this%evol) then
      call write_dbg(cls_name, sname, cls_level, 'ends')
      call stop_tprof( 'push 3D particles' )
      return
   end if
   call this%pd%push(ef,bf)
   call this%pd%pmv(ef,rtag,stag,sid)

   call stop_tprof( 'push 3D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')
   
end subroutine push_beam3d
!
subroutine writehdf5_beam3d(this,file,dspl,rtag,stag,id)

   implicit none
   
   class(beam3d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
   integer, intent(in) :: dspl, rtag, stag
   integer, intent(inout) :: id
! local data
   character(len=18), save :: sname = 'writehdf5_beam3d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%pd%wr(file,dspl,rtag,stag,id)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_beam3d
!            
subroutine writerst_beam3d(this,file)
   implicit none
   
   class(beam3d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'writerst_beam3d'
   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%pd%wrst(file)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writerst_beam3d
!       
subroutine writeq_beam3d(this,file,rtag,stag,id)
   implicit none
   
   class(beam3d), intent(inout) :: this
   class(hdf5file), intent(in), dimension(:) :: file
   integer, intent(in) :: rtag, stag
   integer, intent(inout) :: id
! local data
   character(len=18), save :: sname = 'writeq_beam3d:'
   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%q%write_hdf5(file,1,rtag,stag,id)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writeq_beam3d
!
subroutine readrst_beam3d(this,file)
   implicit none
   
   class(beam3d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'readrst_beam3d'
   
   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%pd%rrst(file)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine readrst_beam3d
!            
end module beam3d_class