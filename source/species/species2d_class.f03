module species2d_class

use parallel_module
use param
use sysutil_module
use options_class
use fdist2d_class
use field_psi_class
use field_e_class
use field_b_class
use field_laser_class
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

   class(part2d), public, pointer :: part => null()
   class(field_rho), allocatable :: q, qn
   class(field_jay), allocatable :: cu, amu
   class(field_djdxi), allocatable :: dcu
   class(fdist2d), pointer :: pf => null()
   integer :: part_shape
   integer :: push_type

   contains

   procedure :: alloc => alloc_species2d
   procedure :: new   => init_species2d
   procedure :: renew => renew_species2d
   procedure :: del   => end_species2d
   procedure :: qdp   => qdp_species2d
   procedure :: deposit_chi => deposit_chi_species2d
   procedure :: amjdp => amjdp_species2d
   procedure :: push_u => push_u_species2d
   procedure :: push_x => push_x_species2d
   procedure :: interp_psi => interp_psi_species2d
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
   push_type, smooth_order )

   implicit none

   class(species2d), intent(inout) :: this
   type(options), intent(in) :: opts
   class(fdist2d), intent(inout), target :: pf
   real, intent(in) :: qbm, s
   integer, intent(in) :: part_shape, max_mode, push_type
   integer, intent(in), optional :: smooth_order
   ! local data
   real :: dt
   logical :: ntd
   integer :: smth_ord
   character(len=18), save :: sname = 'init_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%pf        => pf
   dt             = opts%get_dxi()
   this%push_type = push_type
   this%part_shape = part_shape
   ntd            = pf%neutralized

   allocate( this%q, this%cu, this%amu, this%dcu )

   smth_ord = 0
   if (present(smooth_order)) smth_ord = smooth_order

   call this%q%new(opts,max_mode,this%part_shape,smth_ord)
   call this%cu%new(opts,max_mode,this%part_shape,smth_ord)
   call this%dcu%new(opts,max_mode,this%part_shape,smth_ord)
   call this%amu%new(opts,max_mode,this%part_shape,smth_ord)
   call this%part%new(opts,pf,qbm,dt,s)

   this%q  = 0.0
   this%cu = 0.0

   ! deposit charge
   call this%part%qdeposit( this%q, this%part_shape )
   call this%q%acopy_gc_f1( dir=p_mpi_forward )
   call this%q%copy_gc_f1()
   if (id_stage() == 0) call this%q%copy_slice( 1, p_copy_1to2 )

   ! initialize the background
   if ( ntd ) then

      allocate( this%qn )
      smth_ord = 0
      if (present(smooth_order)) smth_ord = smooth_order
      call this%qn%new(opts,max_mode,part_shape,smth_ord)

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
   call this%part%qdeposit( this%q, this%part_shape )
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

subroutine qdp_species2d(this,q)
! deposit the charge density

  implicit none

  class(species2d), intent(inout) :: this
  class(field_rho), intent(inout) :: q
  ! local data
  character(len=18), save :: sname = 'qdp_species2d'

  call write_dbg(cls_name, sname, cls_level, 'starts')

  this%q = 0.0
  call this%part%qdeposit(this%q, this%part_shape)
  call this%q%acopy_gc_f1( dir=p_mpi_forward )
  call this%q%copy_gc_f1()
  call this%q%smooth()
  call add_f1( this%q, q )
  if ( this%pf%neutralized ) call add_f1( this%qn, q )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine qdp_species2d

subroutine deposit_chi_species2d( this, chi )
! deposit the plasma susceptibility

  implicit none

  class(species2d), intent(inout) :: this
  class(field), intent(inout) :: chi
  ! local data
  character(len=32), save :: sname = 'deposit_chi_species2d'

  call write_dbg(cls_name, sname, cls_level, 'starts')

  call this%part%deposit_chi( chi, this%part_shape )
  ! call this%q%acopy_gc_f1( dir=p_mpi_forward )
  ! call this%q%smooth_f1()
  ! call this%q%copy_gc_f1()
  ! call add_f1( this%q, q )
  ! if ( this%pf%neutralized ) call add_f1( this%qn, q )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine deposit_chi_species2d

subroutine amjdp_species2d( this, ef, bf, af, cu, amu, dcu, dt )
! deposit the current, acceleration and momentum flux

   implicit none

   class(species2d), intent(inout) :: this
   class(field_jay), intent(inout) :: cu, amu
   class(field_djdxi), intent(inout) :: dcu
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf
   class(field_laser), intent(in) :: af
   real, intent(in) :: dt
   ! local data
   character(len=18), save :: sname = 'amjdp_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%cu = 0.0
   this%dcu = 0.0
   this%amu = 0.0
   select case ( this%push_type )
      case (p_push2_std)
         call this%part%amjdeposit_std( ef, bf, this%cu, this%amu, this%dcu, dt, this%part_shape )
      case ( p_push2_robust )
         call this%part%amjdeposit_robust( ef, bf, this%cu, this%amu, this%dcu, dt, this%part_shape )
      case ( p_push2_std_pgc )
         call this%part%amjdeposit_std_pgc( ef, bf, af, this%cu, this%amu, this%dcu, dt, this%part_shape )
      case ( p_push2_robust_pgc )
         call this%part%amjdeposit_robust_pgc( ef, bf, af, this%cu, this%amu, this%dcu, dt, this%part_shape )
   end select

   call this%cu%acopy_gc_f1( dir=p_mpi_forward )
   call this%dcu%acopy_gc_f1( dir=p_mpi_forward )
   call this%amu%acopy_gc_f1( dir=p_mpi_forward )
   call this%cu%copy_gc_f1()
   call this%dcu%copy_gc_f1()
   call this%amu%copy_gc_f1()
   call this%cu%smooth()
   call this%dcu%smooth()
   call this%amu%smooth()

   call add_f1( this%cu, cu )
   call add_f1( this%dcu, dcu )
   call add_f1( this%amu, amu )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine amjdp_species2d

subroutine push_u_species2d(this, ef, bf, af, dt)

   implicit none

   class(species2d), intent(inout) :: this
   class(field_e), intent(in) :: ef
   class(field_b), intent(in) :: bf
   class(field_laser), intent(in) :: af
   real, intent(in) :: dt
   ! local data
   character(len=18), save :: sname = 'push_u_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   select case (this%push_type)
      case (p_push2_std)
         call this%part%push_u_std(ef, bf, dt, this%part_shape)
      case (p_push2_robust)
         call this%part%push_u_robust(ef, bf, dt, this%part_shape)
      case (p_push2_std_pgc)
         call this%part%push_u_std_pgc(ef, bf, af, dt, this%part_shape)
      case (p_push2_robust_pgc)
         call this%part%push_u_robust_pgc(ef, bf, af, dt, this%part_shape)
   end select

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_u_species2d

subroutine push_x_species2d(this, dt)

   implicit none

   class(species2d), intent(inout) :: this
   real, intent(in) :: dt
   ! local data
   character(len=18), save :: sname = 'push_x_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call this%part%push_x(dt)
   call this%part%update_bound()
   call move_part2d_comm( this%part )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_x_species2d

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
   character(len=18), save :: sname = 'cbq_species2d'

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

subroutine interp_psi_species2d(this, psi)

   implicit none

   class(species2d), intent(inout) :: this
   type(field_psi), intent(in) :: psi
   ! local data
   character(len=18), save :: sname = 'interp_psi_species2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')
   if (this%push_type == p_push2_std .or. this%push_type == p_push2_std_pgc) then
      call this%part%interp_psi(psi%get_rf_re(), psi%get_rf_im(), psi%get_max_mode(), this%part_shape)
   endif
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine interp_psi_species2d

end module species2d_class