module sysutil_module

! use, intrinsic :: iso_fortran_env
use mpi

implicit none

private

public :: init_stdout, write_stdout
public :: init_errors, end_errors, set_monitor, write_err, write_wrn, write_dbg
public :: num2str
public :: init_tprof, write_tprof, start_tprof, stop_tprof

integer, save :: class_monitor = 0
integer, save :: fid_err
integer, dimension(4), save :: itime
double precision, save :: dtime
logical, save :: is_master

! variables for timing
integer, parameter :: p_max_event = 128
double precision, dimension(p_max_event), save :: t_event
character(len=64), dimension(p_max_event), save :: name_event
integer, save :: num_event
logical, save :: if_timing = .false.

interface init_stdout
  module procedure init_stdout
end interface

interface write_stdout
  module procedure write_stdout
end interface

interface init_errors
  module procedure init_errors
end interface

interface end_errors
  module procedure end_errors
end interface

interface write_err
  module procedure write_err
end interface

interface write_wrn
  module procedure write_wrn
end interface

interface write_dbg
  module procedure write_dbg
end interface

interface num2str
  module procedure num2str_int
  module procedure num2str_real
end interface

interface init_tprof
  module procedure init_tprof
end interface

interface start_tprof
  module procedure start_tprof
end interface start_tprof

interface stop_tprof
  module procedure stop_tprof
end interface stop_tprof

interface write_tprof
  module procedure write_tprof
end interface write_tprof

contains

subroutine init_stdout( idproc )

  implicit none

  integer, intent(in) :: idproc

  if ( idproc == 0 ) then
    is_master = .true.
  else
    is_master = .false.
  endif

end subroutine init_stdout

subroutine write_stdout( msg )

  implicit none

  character(len=*), intent(in) :: msg

  if ( is_master ) then
    write( *, * ) trim(adjustl(msg))
  endif

end subroutine write_stdout

subroutine init_errors( eunit, idproc, monitor )

  implicit none

  integer, intent(in) :: eunit, idproc, monitor

  character(len=32) :: filename
  integer :: ierr

  fid_err = eunit
  class_monitor = monitor

  if ( idproc == 0 ) then
    call system( 'mkdir ./ELOG' )
  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  write( filename, '(A,I0.6,A)') './ELOG/elog-', idproc, '.log'
  open( unit=fid_err, file=trim(filename), form='formatted', status='replace' )

  call dtimer( dtime, itime, -1 )

end subroutine init_errors

subroutine end_errors()

  implicit none

  close( unit=fid_err )

end subroutine

subroutine write_err( estr )

  implicit none

  character(len=*), intent(in) :: estr
  
  integer :: ierr
  logical :: flag

  call dtimer( dtime, itime, 1 )
  write( fid_err, '(A, F12.3, A12, A)' ) 't = ', dtime, ', [ERROR] ', trim(adjustl(estr))
  
  call mpi_initialized( flag, ierr )
  if ( flag ) then
    call mpi_finalize( ierr )
  else
    stop
  endif

end subroutine write_err

subroutine set_monitor( monitor )

  implicit none

  integer, intent(in) :: monitor

  class_monitor = monitor

end subroutine set_monitor

subroutine write_wrn( wstr )

  implicit none

  character(len=*), intent(in) :: wstr

  call dtimer( dtime, itime, 1 )
  write( fid_err, '(A, F12.3, A12, A)' ) 't = ', dtime, ', [WARNING] ', trim(adjustl(wstr))

end subroutine write_wrn

subroutine write_dbg( clsname, sname, level, msg )

  implicit none

  character(len=*), intent(in) :: clsname, sname
  integer, intent(in) :: level
  character(len=*), intent(in), optional :: msg

  character(len=32) , save :: prefix
  character(len=128), save :: str
  integer :: i

  call dtimer( dtime, itime, 1 )

  if ( level <= class_monitor ) then

    if ( level >= 1 ) then
      i = 1
      prefix = ''
      do while ( i < (level-1)*3 )
        prefix( i:i+2 ) = '|  '
        i = i + 3
      enddo
      prefix(i:i+2) = '|--'
    else
      prefix = ''
    endif

    str = ''
    write( str, '(A, F12.3, A12, A, A, A, A)' ) &
      't = ', dtime, ', [DEBUG] ', trim(prefix), trim(clsname), ' -> ', trim(sname)
    if ( present(msg) ) then
      write( fid_err, * ) trim(str) // ': ' // trim(msg)
    else
      write( fid_err, * ) trim(str)
    endif
    flush( fid_err )

  endif

end subroutine write_dbg

function num2str_int( number, width ) result( str )

  implicit none

  integer, intent(in) :: number
  integer, intent(in), optional :: width
  character(len=:), allocatable :: str

  integer :: clen
  character(len=32) :: tmp_str, format_str

  if ( present(width) ) then
    write( format_str, * ) width
    write( format_str, * ) '(I'//trim(adjustl(format_str))//'.'&
    //trim(adjustl(format_str))//')'
    write( tmp_str, format_str ) number
  else
    write( tmp_str, * ) number
  endif
  clen = len( trim(adjustl(tmp_str)) )
  allocate( character(len=clen) :: str )
  str = trim(adjustl(tmp_str))

end function num2str_int

function num2str_real( number, prec ) result( str )

  implicit none

  real, intent(in) :: number
  integer, intent(in), optional :: prec
  character(len=:), allocatable :: str

  integer :: clen, prec_ = 4
  character(len=32) :: tmp_str, format_str

  if ( present(prec) ) prec_ = prec

  write( format_str, * ) prec_
  write( format_str, * ) '(ES20.'//trim(adjustl(format_str))//')'

  write( tmp_str, format_str ) number
  clen = len_trim( adjustl(tmp_str) )

  allocate( character(len=clen) :: str )

  str = trim(adjustl(tmp_str))

end function num2str_real

subroutine init_tprof( enable_timing )

  implicit none

  logical, intent(in) :: enable_timing

  t_event = 0.0
  num_event = 0

  call add_tprof( 'solve psi' )
  call add_tprof( 'solve bz' )
  call add_tprof( 'solve beam bt' )
  call add_tprof( 'solve plasma bt' )
  call add_tprof( 'solve ez' )
  call add_tprof( 'solve beam et' )
  call add_tprof( 'solve plasma et' )
  call add_tprof( 'solve plasma A' )
  call add_tprof( 'smooth' )
  call add_tprof( 'copy guard cells' )
  call add_tprof( 'copy & add guard cells' )
  call add_tprof( 'copy slices' )
  call add_tprof( 'pipeline' )
  call add_tprof( 'deposit 2D particles' )
  call add_tprof( 'push 2D particles' )
  call add_tprof( 'move 2D particles' )
  call add_tprof( 'deposit 3D particles' )
  call add_tprof( 'push 3D particles' )
  call add_tprof( 'move 3D particles' )
  call add_tprof( 'extract psi' )
  call add_tprof( 'arithmetics' )
  call add_tprof( 'write hdf5' )
  call add_tprof( 'total simulation time' ) ! this must be put in the last

  if_timing = enable_timing

end subroutine init_tprof

subroutine add_tprof( event )

  implicit none

  character(len=*), intent(in) :: event

  num_event = num_event + 1
  name_event(num_event) = trim(event)

end subroutine add_tprof

subroutine start_tprof( event )

  implicit none

  character(len=*), intent(in) :: event

  integer :: i

  if ( (num_event == 0) .or. (if_timing .eqv. .false.) ) return

  do i = 1, num_event
    if ( trim( name_event(i) ) == trim(event) ) then
      call dtimer( dtime, itime, 1 )
      t_event(i) = t_event(i) - dtime
      return
    endif
  enddo

end subroutine start_tprof

subroutine stop_tprof( event )

  implicit none

  character(len=*), intent(in) :: event

  integer :: i

  if ( (num_event == 0) .or. (if_timing .eqv. .false.) ) return

  do i = 1, num_event
    if ( trim( name_event(i) ) == trim(event) ) then
      call dtimer( dtime, itime, 1 )
      t_event(i) = t_event(i) + dtime
      return
    endif
  enddo

end subroutine stop_tprof

subroutine write_tprof()

  implicit none

  integer :: idproc, nproc, ierr, fid = 10, i
  character(len=32) :: filename
  double precision, dimension(:,:), pointer :: buf => null()
  double precision :: avg, max, min, std, lb, prct, whole

  if ( (num_event == 0) .or. (.not. if_timing) ) return

  call MPI_COMM_RANK( MPI_COMM_WORLD, idproc, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )
  if ( idproc == 0 ) then
    call system( 'mkdir ./TIMING' )
  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  write( filename, '(A,I0.6,A)') './TIMING/tprof-', idproc, '.out'
  open( unit=fid, file=trim(filename), form='formatted', status='replace' )
  write ( fid, * ) repeat( '=', 50 )
  write ( fid, '(A30, A20)' ) "EVENT", "ELAPSE TIME (s)"
  write ( fid, * ) repeat( '-', 50 )

  do i = 1, num_event
    write( fid, '(A30, E20.4)' ) trim(adjustl(name_event(i))), t_event(i)
  enddo

  write ( fid, * ) repeat( '=', 50 )
  close( unit=fid )

  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  fid = 11
  if ( idproc == 0 ) then
    open( unit=fid, file='./TIMING/stats-tprof.out', form='formatted', status='replace' )
    write ( fid, * ) repeat( '=', 120 )
    write ( fid, '(A25, A15, A10, 4A15)' ) &
      "EVENT", "AVG TIME (s)", "(%)", "MAX TIME (s)", "MIN TIME (s)", "STD ERR (s)", "BALANCE (%)"
    write ( fid, * ) repeat( '-', 120 )
  endif

  allocate( buf( nproc, num_event ) )

  do i = 1, num_event
    call MPI_GATHER( t_event(i), 1, MPI_DOUBLE_PRECISION, &
      buf(1,i), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  enddo

  if ( idproc == 0 ) then
    whole = sum( buf(:,num_event) )
    do i = 1, num_event
      avg  = sum(buf(:,i)) / nproc
      prct = avg * nproc / whole * 100.0
      max  = maxval(buf(:,i))
      min  = minval(buf(:,i))
      std  = sqrt( sum( (buf(:,i) - avg)**2 ) / nproc )
      lb   = std / avg * 100.0
      write( fid, '(A25, E15.4, F10.1, 3E15.4, F15.1)' ) &
        trim(adjustl(name_event(i))), avg, prct, max, min, std, lb
    enddo
    write ( fid, * ) repeat( '=', 120 )
    close( unit=fid )
  endif

end subroutine write_tprof

end module sysutil_module
