module system

implicit none

private

public :: init_errors, end_errors, write_err, write_wrn, write_dbg

integer, save :: class_monitor = 0
integer, save :: fid
integer, dimension(4), save :: itime
double precision, save :: dtime

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

contains

subroutine init_errors( eunit, idproc, monitor )

  implicit none

  integer, intent(in) :: eunit, idproc, monitor

  character(len=32) :: filename

  fid = eunit
  class_monitor = monitor

  call system( 'mkdir ./ELOG' )

  write( filename, '(A,I0.6,A)') './ELOG/elog-', idproc, '.log'
  open( unit=fid, file=trim(filename), form='formatted', status='replace' )

  call dtimer( dtime, itime, -1 )

end subroutine init_errors

subroutine end_errors()

  implicit none

  close( unit=fid )

end subroutine

subroutine write_err( estr )

  implicit none

  character(len=*), intent(in) :: estr

  call dtimer( dtime, itime, 1 )
  write( fid, '(A, F12.3, A12, A)' ) 't = ', dtime, ', [ERROR] ', trim(adjustl(estr))
  stop

end subroutine write_err

subroutine write_wrn( wstr )

  implicit none

  character(len=*), intent(in) :: wstr

  call dtimer( dtime, itime, 1 )
  write( fid, '(A, F12.3, A12, A)' ) 't = ', dtime, ', [WARNING] ', trim(adjustl(wstr))

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
      write( fid, * ) trim(str) // ': ' // trim(msg)
    else
      write( fid, * ) trim(str)
    endif
    flush( fid )

  endif

end subroutine write_dbg

end module system