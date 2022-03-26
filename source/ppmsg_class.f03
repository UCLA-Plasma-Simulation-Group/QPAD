module ppmsg_class

use parallel_module, only: ntag
use sysutil_module
use mpi

implicit none

type ppmsg

  integer :: tag = 0
  integer :: id = MPI_REQUEST_NULL

  contains

  procedure :: set_idle
  procedure :: get_tag
  procedure :: wait_task

end type ppmsg

contains

subroutine set_idle( this )
  implicit none
  class( ppmsg ), intent(inout) :: this
  this%id = MPI_REQUEST_NULL
end subroutine set_idle

subroutine get_tag( this )
  implicit none
  class( ppmsg ), intent(inout) :: this
  this%tag = ntag()
end subroutine get_tag

subroutine wait_task( this )
  implicit none
  class( ppmsg ), intent(inout) :: this
  integer :: ierr
  integer, dimension(MPI_STATUS_SIZE) :: stat
  call mpi_wait( this%id, stat, ierr )
  if ( ierr /= 0 ) then
    call write_err( 'MPI error!' )
  endif
end subroutine wait_task

end module ppmsg_class