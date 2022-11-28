module beam3d_class

use parallel_module
use param
use input_class
use sysutil_module
use options_class
use field_src_class
use field_class
use field_em_class
use part3d_class
use part3d_comm
use hdf5io_class
use mpi

use fdist3d_class
use fdist3d_std_class
use fdist3d_rnd_class
use fdist3d_file_class

implicit none

private

public :: beam3d

type beam3d

  ! particle dataset
  class(part3d), pointer :: part => null()

  ! charge density field of this beam
  class(field_rho), allocatable :: q

  ! beam profile
  integer :: pf_type
  class(fdist3d), allocatable :: pf

  ! if update particle phase space each time step
  logical :: evol

  ! type of the pusher
  integer :: push_type

  contains

  procedure :: alloc => alloc_beam3d
  procedure :: new   => init_beam3d
  procedure :: del   => end_beam3d
  procedure :: push  => push_beam3d
  procedure :: qdp   => qdeposit_beam3d
  procedure :: wr    => writehdf5_beam3d
  procedure :: wrq   => writeq_beam3d
  procedure :: wrst  => writerst_beam3d
  procedure :: rrst  => readrst_beam3d

end type

save

character(len=10) :: cls_name = 'beam3d'
integer, parameter :: cls_level = 2

contains

subroutine alloc_beam3d( this, input, opts, beam_id )

  implicit none

  class(beam3d), intent(inout) :: this
  type(input_json), intent(inout) :: input
  type(options), intent(in) :: opts
  integer, intent(in) :: beam_id
  ! local data
  character(len=:), allocatable :: read_str
  character(len=20) :: sect_name
  character(len=32), save :: sname = 'alloc_beam3d'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%part ) ) then
    allocate( part3d :: this%part )
  endif

  ! read beam profile type
  sect_name = 'beam(' // num2str(beam_id) // ')'
  this%pf_type = p_pf3d_std
  if ( input%found( trim(sect_name) // '.profile_type' ) ) then
    call input%get( trim(sect_name) // '.profile_type', read_str )
    select case ( trim(read_str) )
      case ( 'standard' )
        this%pf_type = p_pf3d_std
      case ( 'random' )
        this%pf_type = p_pf3d_rnd
      case ( 'file' )
        this%pf_type = p_pf3d_file
      case default
        call write_err( 'Invalid beam profile type!' )
    end select
  endif

  ! initialize beam profile
  select case ( this%pf_type )
    case ( p_pf3d_std )
      allocate( fdist3d_std :: this%pf )
    case ( p_pf3d_rnd )
      allocate( fdist3d_rnd :: this%pf )
    case ( p_pf3d_file )
      allocate( fdist3d_file :: this%pf )
  end select

  ! initialize beam profile. Note the initialization must be called in the allocation
  ! subroutine before initializing the particle manager.
  call this%pf%new( input, opts, beam_id )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_beam3d

subroutine init_beam3d( this, input, opts, max_mode, part_shape, dt, &
  smooth_type, smooth_order, beam_id )

  implicit none

  class(beam3d), intent(inout) :: this
  type(input_json), intent(inout) :: input
  type(options), intent(in) :: opts
  real, intent(in) :: dt
  integer, intent(in) :: part_shape, max_mode, beam_id
  integer, intent(in) :: smooth_type, smooth_order
  ! local data
  integer :: ierr, id
  character(len=:), allocatable :: read_str
  character(len=20) :: sect_name
  integer, dimension(MPI_STATUS_SIZE) :: istat
  character(len=32), save :: sname = 'init_beam3d'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  sect_name = 'beam(' // num2str(beam_id) // ')'

  ! read particle pusher type
  this%push_type = p_push3_reduced
  if ( input%found( trim(sect_name) // '.push_type' ) ) then
    call input%get( trim(sect_name) // '.push_type', read_str )
    select case ( trim(read_str) )
    case ( 'reduced' )
      this%push_type = p_push3_reduced
    case ( 'boris' )
      this%push_type = p_push3_boris
    case default
      call write_err( 'Invalid pusher type! Only "reduced" and "boris" are supported currently.' )
    end select
  endif

  this%evol = this%pf%evol

  ! initialize charge field
  allocate( this%q )
  call this%q%new( opts, max_mode, part_shape, smooth_type, smooth_order )

  ! initialize particles
  if ( this%pf%has_spin ) then
    call this%part%new( opts, this%pf%npmax, dt, this%pf%z0, this%pf%qbm, this%pf%amm )
  else
    call this%part%new( opts, this%pf%npmax, dt, this%pf%z0, this%pf%qbm )
  endif

  ! inject particles
  call this%pf%inject( this%part )
  call move_part3d_comm( this%part, 1, id )
  call mpi_wait( id, istat, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_beam3d

subroutine end_beam3d(this)

  implicit none

  class(beam3d), intent(inout) :: this
  character(len=32), save :: sname = 'end_beam3d'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call this%part%del()
  call this%q%del()
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_beam3d

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
  if ( id_stage() > 0 ) call this%q%pipe_recv( tag, 'forward', 'inner', 'add' )
  call MPI_WAIT(sid, istat, ierr)
  call this%part%qdeposit(this%q)
  call this%q%acopy_gc_f2()
  call this%q%copy_gc_f2()
  call this%q%pipe_send( tag, sid, 'forward', 'guard' )

  call add_f2( this%q, q )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine qdeposit_beam3d

subroutine push_beam3d( this, ef, bf, tag, sid )

  implicit none

  class(beam3d), intent(inout) :: this
  class(field_e), intent(in) :: ef
  class(field_b), intent(in) :: bf
  integer, intent(in) :: tag
  integer, intent(inout) :: sid
  ! local data
  character(len=32), save :: sname = 'push_beam3d'

  call write_dbg(cls_name, sname, cls_level, 'starts')

  if (.not. this%evol) then
    call write_dbg(cls_name, sname, cls_level, 'ends')
    return
  endif

  select case ( this%push_type )
  case ( p_push3_reduced )
    call this%part%push_reduced( ef, bf )
  case ( p_push3_boris )
    call this%part%push_boris( ef, bf )
  end select

  call this%part%update_bound()
  call move_part3d_comm( this%part, tag, sid )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_beam3d

subroutine writehdf5_beam3d( this, file, dspl, rtag, stag, id )

  implicit none

  class(beam3d), intent(inout) :: this
  class(hdf5file), intent(in) :: file
  integer, intent(in) :: dspl, rtag, stag
  integer, intent(inout) :: id
  ! local data
  character(len=32), save :: sname = 'writehdf5_beam3d'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call this%part%wr( file, dspl, rtag, stag, id )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_beam3d

subroutine writerst_beam3d( this, file )
  implicit none

  class(beam3d), intent(inout) :: this
  class(hdf5file), intent(in) :: file
  ! local data
  character(len=32), save :: sname = 'writerst_beam3d'
  call write_dbg(cls_name, sname, cls_level, 'starts')
  call this%part%wrst(file)
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writerst_beam3d

subroutine writeq_beam3d( this, file, rtag, stag, id )
  implicit none

  class(beam3d), intent(inout) :: this
  class(hdf5file), intent(in), dimension(:) :: file
  integer, intent(in) :: rtag, stag
  integer, intent(inout) :: id
  ! local data
  character(len=32), save :: sname = 'writeq_beam3d'
  call write_dbg(cls_name, sname, cls_level, 'starts')
  call this%q%write_hdf5(file, 1, rtag, stag, id)
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writeq_beam3d

subroutine readrst_beam3d( this, file )
  implicit none

  class(beam3d), intent(inout) :: this
  class(hdf5file), intent(in) :: file
  ! local data
  character(len=32), save :: sname = 'readrst_beam3d'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call this%part%rrst(file)
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine readrst_beam3d

end module beam3d_class