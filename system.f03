module system

implicit none

private

integer, save :: class_moniter = 0
integer, save :: fid
integer, dimension(4), save :: itime
double precision, save :: dtime

interface init_errors
  module procedure init_errors
end interface

interface end_errors
end interface

interface ERROR
  module procedure write_err
end interface

interface WARNING
  module procedure write_wrn
end interface WARNING

interface DEBUG
  module procedure write_dbg
end interface

public :: init_errors, end_errors, ERROR, WARNING, DEBUG

contains

subroutine init_errors( eunit, moniter )

  implicit none

  ! class( parallel ), intent(in), pointer :: prl
  integer, intent(in) :: eunit, moniter

  fid = eunit
  class_moniter = moniter

  call system( 'mkdir ./ELOG' )

  open( unit=fid, file='./ELOG/elog-000', form='formatted', status='replace' )

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

  character(len=32) :: prefix = ''
  character(len=128) :: str

  call dtimer( dtime, itime, 1 )

  if ( level <= class_moniter ) then

    if ( level == 0 ) then
      prefix = trim(prefix)
    elseif ( level == 1 ) then
      prefix = '|--'
      prefix = trim(prefix) // ' '
    elseif ( level > 1 ) then
      prefix( (level-1)*3+1:level ) = '|--'
      prefix = trim(prefix) // ' '
    endif

    write( str, '(A, F12.3, A12, A, A, A, A)' ) &
      't = ', dtime, ', [DEBUG] ', prefix, trim(clsname), ' -> ', trim(sname)
    if ( present(msg) ) then
      write( fid, * ) trim(str) // ': ' // trim(msg)
    else
      write( fid, * ) trim(str)
    endif
    flush( fid )

  endif

end subroutine write_dbg

end module system