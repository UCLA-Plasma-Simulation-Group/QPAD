! species2d class for QPAD

module species2d_class

use parallel_pipe_class
use param
use sysutil
use grid_class
use fdist2d_class
use field_psi_class
use field_e_class
use field_b_class
use field_src_class
use field_class
use part2d_class
use hdf5io_class

implicit none

private

public :: species2d

type species2d

   private

   class(part2d), pointer :: pd => null()
   class(field_rho), allocatable :: q, qn
   class(field_jay), allocatable :: cu, amu
   class(field_djdxi), allocatable :: dcu
   class(fdist2d), pointer :: pf => null()
   class(parallel_pipe), pointer :: pp => null()

   contains

   generic :: new => init_species2d
   generic :: renew => renew_species2d
   generic :: del => end_species2d
   generic :: qdp => qdp_species2d
   generic :: amjdp => amjdp_species2d
   generic :: push => push_species2d
   generic :: extpsi => extpsi_species2d
   generic :: psend => psend_species2d
   generic :: precv => precv_species2d
   generic :: wr => writehdf5_species2d
   generic :: wrq => writeq_species2d
   generic :: cbq => cbq_species2d
   procedure, private :: init_species2d, renew_species2d
   procedure, private :: end_species2d
   procedure, private :: qdp_species2d
   procedure, private :: amjdp_species2d
   procedure, private :: push_species2d
   procedure, private :: extpsi_species2d
   procedure, private :: psend_species2d
   procedure, private :: precv_species2d, writehdf5_species2d
   procedure, private :: cbq_species2d, writeq_species2d

end type

save

character(len=10) :: cls_name = 'species2d'
integer, parameter :: cls_level = 2

contains
!
subroutine init_species2d(this,pp,gd,pf,part_shape,&
&num_modes,qbm,xdim,s,smooth_type,smooth_order)

   implicit none

   class(species2d), intent(inout) :: this
   class(parallel_pipe), intent(in), pointer :: pp
   class(grid), intent(in), pointer :: gd
   class(fdist2d), intent(inout), target :: pf
   real, intent(in) :: qbm, s
   integer, intent(in) :: xdim, part_shape, num_modes
   integer, intent(in), optional :: smooth_type, smooth_order
! local data
   real :: dt
   character(len=18), save :: sname = 'init_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%pf => pf
   this%pp => pp
   dt = gd%get_dxi()

   allocate(this%pd,this%q,this%qn,this%cu,this%amu,this%dcu)

   if ( present(smooth_type) .and. present(smooth_order) ) then
      call this%q%new(pp,gd,num_modes,part_shape,smooth_type,smooth_order)
      call this%qn%new(pp,gd,num_modes,part_shape,smooth_type,smooth_order)
      call this%cu%new(pp,gd,num_modes,part_shape,smooth_type,smooth_order)
      call this%dcu%new(pp,gd,num_modes,part_shape,smooth_type,smooth_order)
      call this%amu%new(pp,gd,num_modes,part_shape,smooth_type,smooth_order)
   else
      call this%q%new(pp,gd,num_modes,part_shape)
      call this%qn%new(pp,gd,num_modes,part_shape)
      call this%cu%new(pp,gd,num_modes,part_shape)
      call this%dcu%new(pp,gd,num_modes,part_shape)
      call this%amu%new(pp,gd,num_modes,part_shape)
   endif
   call this%pd%new(pp,pf,this%q,qbm,dt,xdim,s)

   this%qn = 0.0
   this%cu = 0.0
   call this%pd%qdp(this%qn)
   call this%qn%acopy_gc_f1( dir=p_mpi_forward )
   call this%qn%copy_gc_f1()
   this%q = this%qn
   if (pp%getstageid() == 0) then
      ! call this%q%smooth(this%q)
      call this%q%copy_slice(1,p_copy_1to2)
   end if
   ! this%qn = this%qn*(-1.0)
   call dot_f1( -1.0, this%qn )
   call write_dbg(cls_name, sname, cls_level, 'ends')
end subroutine init_species2d
!
subroutine end_species2d(this)

   implicit none

   class(species2d), intent(inout) :: this
   character(len=18), save :: sname = 'end_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%pd%del()
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_species2d
!
subroutine renew_species2d(this,s)

   implicit none

   class(species2d), intent(inout) :: this
   real, intent(in) :: s
! local data
   character(len=18), save :: sname = 'renew_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   call this%pd%renew(this%pf,this%qn,s)
   this%qn = 0.0
   call this%pd%qdp(this%qn)
   call this%qn%acopy_gc_f1( dir=p_mpi_forward )
   call this%qn%copy_gc_f1()
   this%q = this%qn
   if (this%pp%getstageid() == 0) then
      ! call this%q%smooth(this%q)
      call this%q%copy_slice(1,p_copy_1to2)
   end if
   ! this%qn = this%qn * (-1.0)
   call dot_f1( -1.0, this%qn )

   call write_dbg(cls_name, sname, cls_level, 'ends')
end subroutine renew_species2d
!
subroutine qdp_species2d(this,q)
! deposit the charge density

   implicit none

   class(species2d), intent(inout) :: this
   class(field_rho), intent(inout) :: q
! local data
   character(len=18), save :: sname = 'qdp_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%q = 0.0
   call this%pd%qdp(this%q)
   call this%q%acopy_gc_f1( dir=p_mpi_forward )
   call this%q%smooth_f1()
   call this%q%copy_gc_f1()
   ! q = this%q + q + this%qn
   call add_f1( this%q, q )
   call add_f1( this%qn, q )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine qdp_species2d
!
subroutine amjdp_species2d(this,ef,bf,cu,amu,dcu)
! deposit the current, acceleration and momentum flux

   implicit none

   class(species2d), intent(inout) :: this
   class(field_jay), intent(inout) :: cu, amu
   class(field_djdxi), intent(inout) :: dcu
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf
! local data
   character(len=18), save :: sname = 'amjdp_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%cu = 0.0
   this%dcu = 0.0
   this%amu = 0.0
   call this%pd%amjdp(ef,bf,this%cu,this%amu,this%dcu)
   call this%cu%acopy_gc_f1( dir=p_mpi_forward )
   call this%dcu%acopy_gc_f1( dir=p_mpi_forward )
   call this%amu%acopy_gc_f1( dir=p_mpi_forward )
   call this%cu%smooth_f1()
   call this%dcu%smooth_f1()
   call this%amu%smooth_f1()
   call this%cu%copy_gc_f1()
   call this%dcu%copy_gc_f1()
   call this%amu%copy_gc_f1()

   call add_f1( this%cu, cu )
   call add_f1( this%dcu, dcu )
   call add_f1( this%amu, amu )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine amjdp_species2d
!
subroutine push_species2d(this,ef,bf)

   implicit none

   class(species2d), intent(inout) :: this
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf
! local data
   character(len=18), save :: sname = 'push_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   call this%pd%push(ef,bf)
   call this%pd%pmv(this%q)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_species2d
!
subroutine extpsi_species2d(this,psi)

   implicit none

   class(species2d), intent(inout) :: this
   class(field_psi), intent(in) :: psi
! local data
   character(len=18), save :: sname = 'extpsi_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'extract psi' )

   call this%pd%extpsi(psi)

   call stop_tprof( 'extract psi' )
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine extpsi_species2d
!
subroutine psend_species2d(this, tag, id)

   implicit none

   class(species2d), intent(inout) :: this
   integer, intent(in) :: tag
   integer, intent(inout) :: id
! local data
   character(len=18), save :: sname = 'pipesend_part2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%pd%psend(tag,id)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine psend_species2d
!
subroutine precv_species2d(this, tag)

   implicit none

   class(species2d), intent(inout) :: this
   integer, intent(in) :: tag
! local data
   character(len=18), save :: sname = 'precv_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%pd%precv(tag)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine precv_species2d
!
subroutine writehdf5_species2d(this,file)

   implicit none

   class(species2d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'writehdf5_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   call this%pd%wr(file)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_species2d
!
subroutine writeq_species2d(this,files,rtag,stag,id)

   implicit none

   class(species2d), intent(inout) :: this
   class(hdf5file), dimension(:), intent(in) :: files
   integer, intent(in) :: rtag, stag
   integer, intent(inout) :: id
! local data
   character(len=18), save :: sname = 'writeq_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   call this%q%write_hdf5(files,1,rtag,stag,id)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writeq_species2d
!
subroutine cbq_species2d(this,pos)

   implicit none

   class(species2d), intent(inout) :: this
   integer, intent(in) :: pos
! local data
   character(len=18), save :: sname = 'cpq_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   ! call add_f1(this%q,this%cu,this%q,(/1/),(/3/),(/1/))
   call add_f1( this%cu, this%q, (/3/), (/1/) )
   call this%q%smooth_f1()
   call this%q%copy_slice(pos,p_copy_1to2)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine cbq_species2d
!
end module species2d_class