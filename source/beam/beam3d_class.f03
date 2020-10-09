module beam3d_class

use parallel_module
use param
use sysutil_module
use fdist3d_class
use options_class
use field_src_class
use field_class
use field_e_class
use field_b_class
use part3d_class
use hdf5io_class
use mpi
use part3d_comm

implicit none

private

public :: beam3d

type beam3d

   private

   class(part3d), pointer :: part => null()
   class(field_rho), allocatable :: q
   class(fdist3d), pointer :: pf => null()
   logical :: evol
   integer :: push_type
   contains

   procedure :: new  => init_beam3d
   procedure :: del  => end_beam3d
   procedure :: push => push_beam3d
   procedure :: wr   => writehdf5_beam3d
   procedure :: wrq  => writeq_beam3d
   procedure :: wrst => writerst_beam3d
   procedure :: rrst => readrst_beam3d

   generic :: qdp  => qdeposit_beam3d, qdeposit_beam3d_finish
   procedure, private :: qdeposit_beam3d, qdeposit_beam3d_finish
end type

save

character(len=10) :: cls_name = 'beam3d'
integer, parameter :: cls_level = 2

contains
!
subroutine init_beam3d( this, opts, max_mode, part_shape, pf, qbm, dt, &
   push_type, smooth_type, smooth_order, has_spin, amm )

   implicit none

   class(beam3d), intent(inout) :: this
   type(options), intent(in) :: opts
   class(fdist3d), intent(inout), target :: pf
   real, intent(in) :: qbm, dt
   integer, intent(in) :: push_type, part_shape, max_mode
   integer, intent(in) :: smooth_type, smooth_order
   logical, intent(in) :: has_spin
   real, intent(in) :: amm
! local data
   character(len=32), save :: sname = 'init_beam3d'
   integer :: id, ierr
   integer, dimension(10) :: istat

   call write_dbg(cls_name, sname, cls_level, 'starts')
   this%pf => pf
   this%push_type = push_type

   allocate(this%part,this%q)
   this%evol = pf%getevol()
   call this%q%new(opts,max_mode,part_shape,smooth_type,smooth_order)
   call this%part%new(opts,pf,qbm,dt,has_spin,amm)
   ! call this%part%pmv(this%q,1,1,id)
   call move_part3d_comm( this%part, 1, 1, id )
   call MPI_WAIT(id,istat,ierr)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_beam3d
!
subroutine end_beam3d(this)

   implicit none

   class(beam3d), intent(inout) :: this
   character(len=32), save :: sname = 'end_beam3d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%part%del()
   call this%q%del()
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_beam3d
!
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
      call this%part%qdeposit(this%q)
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

   select case ( this%push_type )
   case ( p_push_reduced )
      call this%part%push_reduced( ef, bf )
   case ( p_push_boris )
      call this%part%push_boris( ef, bf )
   end select

   call this%part%update_bound()
   ! call this%part%pmv(ef,rtag,stag,sid)
   call move_part3d_comm( this%part, rtag, stag, sid )

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
   call this%part%wr(file,dspl,rtag,stag,id)
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
   call this%part%wrst(file)
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
   call this%part%rrst(file)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine readrst_beam3d
!
end module beam3d_class