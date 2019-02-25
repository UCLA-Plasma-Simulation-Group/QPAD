module system

! use, intrinsic :: iso_fortran_env
use mpi

implicit none

private

public :: init_stdout, write_stdout
public :: init_errors, end_errors, write_err, write_wrn, write_dbg
public :: num2str

integer, save :: class_monitor = 0
integer, save :: fid_err, fid_stdout
integer, dimension(4), save :: itime
double precision, save :: dtime
logical, save :: is_master

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

  call dtimer( dtime, itime, 1 )
  write( fid_err, '(A, F12.3, A12, A)' ) 't = ', dtime, ', [ERROR] ', trim(adjustl(estr))
  stop

end subroutine write_err

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

end module system