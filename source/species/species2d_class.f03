module species2d_class

use parallel_module
use param
use sysutil_module
use options_class
use fdist2d_class
use field_psi_class
use field_b_class
use field_e_class
use field_src_class
use field_class
use part2d_class
use part2d_comm
use hdf5io_class

implicit none

private

public :: species2d

type species2d

   private

   class(part2d), pointer :: part => null()
   class(field_rho), allocatable :: q, qn
   class(field_jay), allocatable :: cu, amu
   class(field_djdxi), allocatable :: dcu
   class(fdist2d), pointer :: pf => null()
   integer :: push_type

   contains

   procedure :: alloc => alloc_species2d
   procedure :: new   => init_species2d
   procedure :: renew => renew_species2d
   procedure :: del   => end_species2d
   procedure :: qdp   => qdp_species2d
   procedure :: amjdp => amjdp_species2d
   procedure :: edp   => edp_species2d
   procedure :: push  => push_species2d
   procedure :: epush => push_species2d_explicit
   procedure :: psend => psend_species2d
   procedure :: precv => precv_species2d
   procedure :: wr    => writehdf5_species2d
   procedure :: wrq   => writeq_species2d
   procedure :: cbq   => cbq_species2d
   procedure :: sort  => sort_species2d

end type

save

character(len=10) :: cls_name = 'species2d'
integer, parameter :: cls_level = 2

contains

subroutine alloc_species2d( this )

   implicit none

   class(species2d), intent(inout) :: this

   if ( .not. associated( this%part ) ) then
      allocate( part2d :: this%part )
   endif

end subroutine alloc_species2d

subroutine init_species2d( this, opts, pf, part_shape, max_mode, qbm, s, &
   push_type, smooth_type, smooth_order )

   implicit none

   class(species2d), intent(inout) :: this
   type(options), intent(in) :: opts
   class(fdist2d), intent(inout), target :: pf
   real, intent(in) :: qbm, s
   integer, intent(in) :: part_shape, max_mode, push_type
   integer, intent(in), optional :: smooth_type, smooth_order
   ! local data
   real :: dt
   logical :: ntd
   character(len=18), save :: sname = 'init_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%pf        => pf
   dt             = opts%get_dxi()
   this%push_type = push_type
   ntd            = pf%neutralized

   allocate( this%q, this%cu, this%amu, this%dcu )

   if ( present(smooth_type) .and. present(smooth_order) ) then
      call this%q%new(opts,max_mode,part_shape,smooth_type,smooth_order)
      call this%cu%new(opts,max_mode,part_shape,smooth_type,smooth_order)
      call this%dcu%new(opts,max_mode,part_shape,smooth_type,smooth_order)
      call this%amu%new(opts,max_mode,part_shape,smooth_type,smooth_order)
   else
      call this%q%new(opts,max_mode,part_shape)      
      call this%cu%new(opts,max_mode,part_shape)
      call this%dcu%new(opts,max_mode,part_shape)
      call this%amu%new(opts,max_mode,part_shape)
   endif
   call this%part%new(opts,pf,qbm,dt,s)

   this%q  = 0.0
   this%cu = 0.0
!    this%amu = 0.0
!    this%ve = 0.0
!    this%gam = 0.0

   ! deposit charge
   call this%part%qdeposit( this%q )
   call this%q%acopy_gc_f1( dir=p_mpi_forward )
   call this%q%copy_gc_f1()
   if (id_stage() == 0) call this%q%copy_slice( 1, p_copy_1to2 )

   ! initialize the background
   if ( ntd ) then

      allocate( this%qn )
      if ( present(smooth_type) .and. present(smooth_order) ) then
         call this%qn%new(opts,max_mode,part_shape,smooth_type,smooth_order)
      else
         call this%qn%new(opts,max_mode,part_shape)
      endif

      this%qn = this%q
      call dot_f1( -1.0, this%qn )

   endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_species2d

subroutine end_species2d(this)

   implicit none

   class(species2d), intent(inout) :: this
   character(len=18), save :: sname = 'end_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%part%del()
   call this%q%del() 
   call this%cu%del() 
   call this%dcu%del()
   call this%amu%del() 
   if ( allocated( this%qn ) ) call this%qn%del()
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_species2d

subroutine renew_species2d(this,s)

   implicit none

   class(species2d), intent(inout) :: this
   real, intent(in) :: s
   ! local data
   logical :: ntd
   character(len=18), save :: sname = 'renew_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ntd = this%pf%neutralized

   call this%part%renew(this%pf,s)

   ! renew the charge density
   this%q = 0.0
   call this%part%qdeposit( this%q )
   call this%q%acopy_gc_f1( dir=p_mpi_forward )
   call this%q%copy_gc_f1()
   if (id_stage() == 0) call this%q%copy_slice( 1, p_copy_1to2 )

   ! renew the background
   if ( ntd ) then
      this%qn = this%q
      call dot_f1( -1.0, this%qn )
   endif

   call write_dbg(cls_name, sname, cls_level, 'ends')
end subroutine renew_species2d

subroutine qdp_species2d(this,q, qn)
! deposit the charge density

   implicit none

   class(species2d), intent(inout) :: this
   class(field_rho), intent(inout) :: q
   class(field_rho), intent(inout) :: qn
   ! local data
   character(len=18), save :: sname = 'qdp_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%q = 0.0
!    write(2,*)  this%q%getresum(), "q 1"
   call this%part%qdeposit(this%q)
!    write(2,*)  this%q%getresum(), "q 2"
   call this%q%acopy_gc_f1( dir=p_mpi_forward )
   call this%q%smooth_f1()
   call this%q%copy_gc_f1()
   call add_f1( this%q, q )
!    write(2,*)  this%q%getresum(), "q"
!    write(2,*)  this%qn%getresum(), "qn"
   if ( this%pf%neutralized ) call add_f1( this%qn, q )
   if ( this%pf%neutralized ) call add_f1( this%qn, qn )
!    write(2,*)  q%getresum(), "q2"

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine qdp_species2d

subroutine amjdp_species2d( this, ef, bf, cu, amu, dcu )
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
   select case ( this%push_type )
      case ( p_push2_robust )
         call this%part%amjdeposit_robust( ef, bf, this%cu, this%amu, this%dcu )
      case ( p_push2_clamp )
         call this%part%amjdeposit_clamp( ef, bf, this%cu, this%amu, this%dcu )
      case ( p_push2_robust_subcyc )
         call this%part%amjdeposit_robust_subcyc( ef, bf, this%cu, this%amu, this%dcu )
   end select

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

subroutine edp_species2d( this, ef, bf, b_beam, cu, amu, dcu )
! deposit the current, acceleration and momentum flux

   implicit none

   class(species2d), intent(inout) :: this
   class(field_jay), intent(inout) :: cu, amu
   class(field_djdxi), intent(inout) :: dcu
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf
   class(field_b), intent(in) :: b_beam
   ! local data
   character(len=18), save :: sname = 'amjdp_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%cu = 0.0
   this%dcu = 0.0
   this%amu = 0.0
!    write(2,*) cu%getresum(),"cu species2d"
   call this%part%edeposit( ef, bf, b_beam, this%cu, this%amu, this%dcu )
!    write(2,*) this%cu%getresum(),"this%cu species2d"
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
!    write(2,*) cu%getresum(),"cu species2d end"
!    write(2,*) this%cu%getresum(),"this%cu species2d end"

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine edp_species2d

subroutine push_species2d(this,ef,bf)

   implicit none

   class(species2d), intent(inout) :: this
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf
   ! local data
   character(len=18), save :: sname = 'push_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   select case ( this%push_type )
      case ( p_push2_robust )
         call this%part%push_robust( ef, bf )
      case ( p_push2_clamp )
         call this%part%push_clamp( ef, bf )
      case ( p_push2_robust_subcyc )
         call this%part%push_robust_subcyc( ef, bf )
   end select

   call this%part%update_bound()
   call move_part2d_comm( this%part )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_species2d

subroutine push_species2d_explicit(this,ef,bf)

   implicit none

   class(species2d), intent(inout) :: this
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf 
   ! local data
   character(len=18), save :: sname = 'push_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   select case ( this%push_type )
      case ( p_push2_robust )
         call this%part%push_robust( ef, bf )
      case ( p_push2_clamp )
         call this%part%push_clamp( ef, bf )
      case ( p_push2_robust_subcyc )
         call this%part%push_robust_subcyc( ef, bf )
      case( p_push2_explicit )
         call this%part%expush( ef, bf)
   end select

   call this%part%update_bound()
   call move_part2d_comm( this%part )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_species2d_explicit

subroutine psend_species2d(this, tag, id)

   implicit none

   class(species2d), intent(inout) :: this
   integer, intent(in) :: tag
   integer, intent(inout) :: id
   ! local data
   character(len=18), save :: sname = 'pipesend_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%part%pipesend(tag,id)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine psend_species2d

subroutine precv_species2d(this, tag)

   implicit none

   class(species2d), intent(inout) :: this
   integer, intent(in) :: tag
   ! local data
   character(len=18), save :: sname = 'precv_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%part%piperecv(tag)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine precv_species2d

subroutine writehdf5_species2d(this,file)

   implicit none

   class(species2d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
   ! local data
   character(len=18), save :: sname = 'writehdf5_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   call this%part%wr(file)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_species2d

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

subroutine cbq_species2d(this,pos)

   implicit none

   class(species2d), intent(inout) :: this
   integer, intent(in) :: pos
   ! local data
   character(len=18), save :: sname = 'cpq_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call add_f1( this%cu, this%q, (/3/), (/1/) )
   call this%q%copy_slice(pos,p_copy_1to2)
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine cbq_species2d

subroutine sort_species2d( this, step2d )

   implicit none

   class(species2d), intent(inout) :: this
   integer, intent(in) :: step2d
   ! local data
   integer :: sort_freq, nrp, noff
   character(len=18), save :: sname = 'sort_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   sort_freq = this%pf%sort_freq
   nrp       = this%q%rf_re(0)%get_ndp(1)
   noff      = this%q%rf_re(0)%get_noff(1)

   if ( sort_freq == 0 ) return
   if ( mod(step2d, sort_freq) == 0 ) then
      call this%part%sort( nrp, noff )
   endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine sort_species2d

end module species2d_class