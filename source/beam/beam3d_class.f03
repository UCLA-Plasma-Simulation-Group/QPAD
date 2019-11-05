! beam3d class for QPAD

module beam3d_class

use parallel_pipe_class
use param
use sys
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

   generic :: new  => init_beam3d
   generic :: del  => end_beam3d
   generic :: push => push_beam3d
   generic :: qdp  => qdeposit_beam3d, qdeposit_beam3d_finish
   generic :: wr   => writehdf5_beam3d
   generic :: wrq  => writeq_beam3d
   generic :: wrst => writerst_beam3d
   generic :: rrst => readrst_beam3d
   procedure, private :: init_beam3d
   procedure, private :: end_beam3d
   procedure, private :: push_beam3d
   procedure, private :: qdeposit_beam3d, qdeposit_beam3d_finish, writehdf5_beam3d
   procedure, private :: writerst_beam3d, readrst_beam3d
   procedure, private :: writeq_beam3d
end type

save

character(len=10) :: cls_name = 'beam3d'
integer, parameter :: cls_level = 2
character(len=128) :: erstr

contains
!
subroutine init_beam3d(this,pp,gd,max_mode,part_shape,pf,qbm,dt,xdim,smooth_type,smooth_order)

   implicit none

   class(beam3d), intent(inout) :: this
   class(grid), intent(in), pointer :: gd
   class(parallel_pipe), intent(in), pointer :: pp
   class(fdist3d), intent(inout), target :: pf
   real, intent(in) :: qbm, dt
   integer, intent(in) :: xdim, part_shape, max_mode
   integer, intent(in), optional :: smooth_type, smooth_order
! local data
   character(len=32), save :: sname = 'init_beam3d'
   integer :: id, num_modes, ierr
   integer, dimension(10) :: istat

   call write_dbg(cls_name, sname, cls_level, 'starts')
   this%pp => pp
   this%pf => pf

   allocate(this%pd,this%q)
   this%evol = pf%getevol()
   if ( present(smooth_type) .and. present(smooth_order) ) then
      call this%q%new(pp,gd,max_mode,part_shape,smooth_type,smooth_order)
   else
      call this%q%new(pp,gd,max_mode,part_shape)
   endif
   call this%pd%new(pp,gd,pf,qbm,dt,xdim)
   call this%pd%pmv(this%q,1,1,id)
   call MPI_WAIT(id,istat,ierr)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_beam3d
!
subroutine end_beam3d(this)

   implicit none

   class(beam3d), intent(inout) :: this
   character(len=32), save :: sname = 'end_beam3d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%pd%del()
   call this%q%del()
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_beam3d
!
! subroutine qdeposit_beam3d(this,q)
! ! deposit the charge density

!    implicit none

!    class(beam3d), intent(inout) :: this
!    class(field_rho), intent(inout) :: q
! ! local data
!    character(len=32), save :: sname = 'qdeposit_beam3d'

!    call write_dbg(cls_name, sname, cls_level, 'starts')

!    if (.not. this%evol) then
!       call this%q%as(0.0)
!       call this%pf%dp(this%q)
!       call this%q%copy_gc_f2()
!       call add_f2( this%q, q )
!    else
!       call this%q%as(0.0)
!       call this%pd%qdp(this%q)
!       call this%q%acopy_gc_f2()
!       call this%q%copy_gc_f2()
!       call add_f2( this%q, q )
!    end if

!    call write_dbg(cls_name, sname, cls_level, 'ends')

! end subroutine qdeposit_beam3d

subroutine qdeposit_beam3d(this,q,tag,sid)
! deposit the charge density

   implicit none

   class(beam3d), intent(inout) :: this
   class(field_rho), intent(inout) :: q
   integer, intent(in) :: tag
   integer, intent(inout) :: sid
! local data
   character(len=32), save :: sname = 'qdeposit_beam3d'
   integer, dimension(10) :: istat
   integer :: ierr

   call write_dbg(cls_name, sname, cls_level, 'starts')

   call this%q%as(0.0)
   call MPI_WAIT(sid, istat, ierr)

   if (.not. this%evol) then
      call this%pf%dp(this%q)
      call this%q%copy_gc_f2()
   else
      call this%pd%qdp(this%q)
      call this%q%acopy_gc_f2()
      call this%q%copy_gc_f2()
   endif

   call add_f2( this%q, q )
   call this%q%pipe_gc_send(tag, sid)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine qdeposit_beam3d

subroutine qdeposit_beam3d_finish(this,tag)
! finish deposit the charge density

   implicit none

   class(beam3d), intent(inout) :: this
   integer, intent(in) :: tag
! local data
   character(len=32), save :: sname = 'qdeposit_beam3d_finish'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%q%pipe_gc_recv(tag)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine qdeposit_beam3d_finish

subroutine push_beam3d(this,ef,bf,rtag,stag,sid)

   implicit none

   class(beam3d), intent(inout) :: this
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf
   integer, intent(in) :: rtag, stag
   integer, intent(inout) :: sid
! local data
   character(len=32), save :: sname = 'partpush'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   if (.not. this%evol) then
      call write_dbg(cls_name, sname, cls_level, 'ends')
      return
   end if
   call this%pd%push(ef,bf)
   call this%pd%pmv(ef,rtag,stag,sid)

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
   character(len=32), save :: sname = 'writehdf5_beam3d'

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
   character(len=32), save :: sname = 'writerst_beam3d'
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
   character(len=32), save :: sname = 'writeq_beam3d:'
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
   character(len=32), save :: sname = 'readrst_beam3d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%pd%rrst(file)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine readrst_beam3d
!
end module beam3d_class