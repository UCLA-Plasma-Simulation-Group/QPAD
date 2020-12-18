module diagnostics_popas_class

use parallel_module
use sysutil_module
use param
use input_class
use options_class
use diagnostics_class
use sim_beams_class
use sim_beams_popas_class
use sim_plasma_class
use sim_fields_class
use field_class
use beam3d_class
use beam3d_popas_class
use species2d_class
use mpi

implicit none

private

public :: sim_diag_popas

integer, parameter :: p_tdiag_popas_emit = 91, p_tdiag_popas_esprd = 92

type, extends( sim_diag ) :: sim_diag_popas

  ! private

  ! class( diag_node ), pointer :: head => null()
  ! class( diag_node ), pointer :: diag => null()
  ! integer :: num_diag = 0
  ! integer :: ndump_gcd = 0 ! greatest common divisor of all the values of ndump
  ! logical :: has_vpotz = .false.
  ! logical :: has_vpott = .false.

  contains

  ! procedure :: alloc => alloc_sim_diag
  procedure :: new => init_sim_diag_popas
  ! procedure :: del => end_sim_diag_popas
  procedure :: run => run_sim_diag_popas
  ! generic :: add_diag => add_diag_cym, add_diag_raw, add_diag_rst

  procedure, private :: add_diag_popas
  ! procedure, private :: add_diag_cym, add_diag_raw, add_diag_rst
  ! procedure, private :: to_next, to_head, to_tail, is_tail
  ! procedure, private :: set_ndump_gcd
  ! procedure, private :: init_diag_beams
  ! procedure, private :: init_diag_species
  ! procedure, private :: init_diag_fields
  procedure, private :: init_diag_popas

end type sim_diag_popas

contains

subroutine init_diag_popas( this, input, beams )

  implicit none

  class( sim_diag_popas ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  class( sim_beams_popas ), intent(in), target :: beams

  ! local data
  integer :: nbeams, ndump
  integer :: i, j, k, m, n
  character(len=:), allocatable :: ss

  call input%get( 'simulation.nbeams', nbeams )

  ! add beam diagnostics
  do i = 1, nbeams
    call input%info( 'beam('//num2str(i)//').diag', n_children=m )
    do j = 1, m
      call input%get( 'beam('//num2str(i)//').diag'//'('//num2str(j)//').ndump', ndump )
      if ( ndump > 0 ) then
        call input%info( 'beam('//num2str(i)//').diag'//'('//num2str(j)//').name', n_children=n )
        do k = 1, n
          if ( allocated(ss) ) deallocate(ss)
          call input%get( 'beam('//num2str(i)//').diag'//'('//num2str(j)//').name'//'('//num2str(k)//')', ss )
          select case ( trim(ss) )
          case ( 'popas_emittance' )
            call this%add_diag_popas( &
              obj        = beams%beam(i), &
              dump_freq  = ndump, &
              type_label = 'popas_emittance', &
              pathname   = './Beam'//num2str(i)//'/POPAS/', &
              dataname   = 'emittance', &
              fid        = 100+i )
          case ( 'popas_ene_spread' )
            call this%add_diag_popas( &
              obj        = beams%beam(i), &
              dump_freq  = ndump, &
              type_label = 'popas_ene_spread', &
              pathname   = './Beam'//num2str(i)//'/POPAS/', &
              dataname   = 'ene_spread', &
              fid        = 200+i )
          end select
        enddo ! end of k
      endif
    enddo ! end of j
  enddo ! end of i

end subroutine init_diag_popas

subroutine init_sim_diag_popas( this, input, opts, fields, beams, plasma )

  implicit none

  class( sim_diag_popas ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  class( sim_fields ), intent(inout) :: fields
  class( sim_beams ), intent(in) :: beams
  class( sim_plasma ), intent(in) :: plasma

  ! local data
  integer :: ierr
  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag_popas'
  character(len=32), save :: sname = 'init_sim_diag_popas'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ! call initialization procedure of the super-class
  call this%sim_diag%new( input, opts, fields, beams, plasma )

  ! initialize POPAS diagnostics
  select type ( obj => beams )
  type is ( sim_beams_popas )
    call this%init_diag_popas( input, obj )
  end select

  call this%set_ndump_gcd()
  call mpi_barrier( comm_world(), ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_diag_popas

subroutine run_sim_diag_popas( this, tstep, dt )

  implicit none

  class( sim_diag_popas ), intent(inout) :: this
  integer, intent(in) :: tstep
  real, intent(in) :: dt

  ! local data
  integer :: stag, rtag, ierr
  integer, dimension(MPI_STATUS_SIZE) :: istat

  integer, save :: cls_level = 3
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'run_sim_diag'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ! run super-class procedure
  call this%sim_diag%run( tstep, dt )

  ! POPAS diagnostics is here!
  call this%to_head()
  if ( .not. associated(this%diag) ) return

  do
    if ( mod( tstep, this%diag%dump_freq ) == 0 ) then

      call this%diag%set_sim_time( tstep, tstep*dt )

      select type ( obj => this%diag%obj )
      type is ( beam3d_popas )

        select case ( trim(this%diag%type_label) )

        case ( 'popas_emittance' )

          rtag = ntag(); stag = rtag
          call mpi_wait( this%diag%id, istat, ierr )
          call obj%get_emittance( this%diag%fid, tstep, rtag, stag, this%diag%id )

        case ( 'popas_ene_spread' )

          rtag = ntag(); stag = rtag
          call mpi_wait( this%diag%id, istat, ierr )
          call obj%get_ene_spread( this%diag%fid, tstep, rtag, stag, this%diag%id )

        end select

      end select

    endif

    if ( .not. this%is_tail() ) then
      call this%to_next()
    else
      exit
    endif

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine run_sim_diag_popas

! add POPAS data for particles
subroutine add_diag_popas( this, obj, dump_freq, type_label, pathname, dataname, &
  fid )

  implicit none

  class( sim_diag_popas ), intent(inout) :: this
  class(*), intent(in) :: obj
  integer, intent(in) :: dump_freq, fid
  character(len=*), intent(in) :: pathname, dataname, type_label

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag_popas'
  character(len=32), save :: sname = 'add_diag_popas'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%head ) ) then
    allocate( this%head )
    call this%head%new( obj, dump_freq, num_files=1, id=MPI_REQUEST_NULL )
    this%diag => this%head
  else
    call this%to_tail()
    allocate( this%diag%next )
    this%diag => this%diag%next
    call this%diag%new( obj, dump_freq, num_files=1, id=MPI_REQUEST_NULL )
  endif

  this%diag%type_label = trim(type_label)
  this%diag%fid = fid

  ! the first processor of the last stage create output files
  if ( id_stage() == num_stages() - 1 .and. id_proc_loc() == 0 ) then
    call system( 'mkdir -p '//trim(pathname) )
    open( unit=fid, file=trim(pathname)//trim(dataname)//'.txt', &
      form='formatted', status='replace' )
    select case ( trim(type_label) )
    case ( 'popas_emittance' )
      write (fid,'(A8,8A22)') 'Step #', '<x^2>', '<px^2>', '<x*px>', 'emit_x', &
        '<y^2>', '<py^2>', '<y*py>', 'emit_y'
    case ( 'popas_ene_spread' )
      write (fid, '(A8,A22)') 'Step #', 'Energy_spread'
    end select
    flush( fid )
  endif

  this%num_diag = this%num_diag + 1

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_diag_popas

end module diagnostics_popas_class