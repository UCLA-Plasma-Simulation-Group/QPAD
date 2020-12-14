module part3d_popas_class

use param
use sysutil_module
use parallel_module
use part3d_class
use mpi

implicit none

private

public :: part3d_popas

type, extends( part3d ) :: part3d_popas

  contains

  ! new procedure
  procedure :: get_emittance => get_emittance_part3d_popas
  procedure :: get_ene_spread => get_ene_spread_part3d_popas

end type

save

character(len=32), parameter :: cls_name = "part3d_popas"
integer, parameter :: cls_level = 2

contains

subroutine get_emittance_part3d_popas( this, fid, tstep, recv_tag, send_tag, id )

  implicit none

  class( part3d_popas ), intent(inout) :: this
  integer, intent(in) :: fid, tstep, recv_tag, send_tag
  integer, intent(inout) :: id

  integer :: nprocs_loc, idproc_orig, idproc_dest, ierr, i
  real, dimension(2) :: avg_x2, avg_p2, avg_xp, emit
  real, dimension(11) :: msg, rslt
  real, dimension(2) :: x, p, x2, p2, xp
  real :: w, w_tot

  nprocs_loc  = num_procs_loc()
  idproc_orig = id_proc() - nprocs_loc
  idproc_dest = id_proc() + nprocs_loc

  ! receive message from the previous stage
  if ( idproc_orig >= 0 ) then

    call mpi_recv( msg, 11, p_dtype_real, idproc_orig, recv_tag, comm_world(), &
      MPI_STATUS_IGNORE, ierr )
    x = msg(1:2); x2 = msg(3:4)
    p = msg(5:6); p2 = msg(7:8)
    xp = msg(9:10); w_tot = msg(11)

  else

    x = 0.0; x2 = 0.0
    p = 0.0; p2 = 0.0
    xp = 0.0; w_tot = 0.0

  endif

  ! calculate the charge-weighted first and second moments, and the total weight
  ! of the local stage
  do i = 1, this%npp

    w = abs( this%q(i) )
    x = x + this%x(1:2,i) * w
    p = p + this%p(1:2,i) * w
    x2(1) = x2(1) + this%x(1,i) * this%x(1,i) * w
    x2(2) = x2(2) + this%x(2,i) * this%x(2,i) * w
    p2(1) = p2(1) + this%p(1,i) * this%p(1,i) * w
    p2(2) = p2(2) + this%p(2,i) * this%p(2,i) * w
    xp(1) = xp(1) + this%x(1,i) * this%p(1,i) * w
    xp(2) = xp(2) + this%x(2,i) * this%p(2,i) * w
    w_tot = w_tot + w

  enddo

  ! pack message for sending
  msg = (/x(1), x(2), x2(1), x2(2), p(1), p(2), p2(1), p2(2), xp(1), xp(2), w_tot/)

  ! send message to next stage
  ! the last stage calculates the results
  if ( idproc_dest < num_procs() ) then
    call mpi_isend( msg, 11, p_dtype_real, idproc_dest, send_tag, comm_world(), &
      id, ierr )
  else

    call mpi_reduce( msg, rslt, 11, p_dtype_real, MPI_SUM, id_proc_loc(), comm_loc(), ierr )

    ! the first processor in the last stage output data
    if ( id_proc_loc() == 0 ) then

      rslt(1:10) = rslt(1:10) / rslt(11)
      avg_x2(1) = rslt(3) - rslt(1) * rslt(1)
      avg_x2(2) = rslt(4) - rslt(2) * rslt(2)
      avg_p2(1) = rslt(7) - rslt(5) * rslt(5)
      avg_p2(2) = rslt(8) - rslt(6) * rslt(6)
      avg_xp(1) = rslt(9) - rslt(1) * rslt(5)
      avg_xp(2) = rslt(10)- rslt(2) * rslt(6)
      emit(1) = sqrt( avg_x2(1) * avg_p2(1) - avg_xp(1) * avg_xp(1) )
      emit(2) = sqrt( avg_x2(2) * avg_p2(2) - avg_xp(2) * avg_xp(2) )

      write( fid, '(I8.8,8D22.14)' ) tstep, avg_x2(1), avg_p2(1), &
        avg_xp(1), emit(1), avg_x2(2), avg_p2(2), avg_xp(2), emit(2)
      flush( fid )
    endif

  endif

end subroutine get_emittance_part3d_popas

subroutine get_ene_spread_part3d_popas( this, fid, tstep, recv_tag, send_tag, id )

  implicit none

  class( part3d_popas ), intent(inout) :: this
  integer, intent(in) :: tstep, recv_tag, send_tag, fid
  integer, intent(inout) :: id

  integer :: nprocs_loc, idproc_orig, idproc_dest, ierr, i
  real, dimension(3) :: msg, rslt
  real :: w, w_tot, gam, gam2, gam_tmp, ene_spread

  nprocs_loc  = num_procs_loc()
  idproc_orig = id_proc() - nprocs_loc
  idproc_dest = id_proc() + nprocs_loc

  ! receive message from the previous stage
  if ( idproc_orig >= 0 ) then

    call mpi_recv( msg, 3, p_dtype_real, idproc_orig, recv_tag, comm_world(), &
      MPI_STATUS_IGNORE, ierr )
    gam   = msg(1)
    gam2  = msg(2)
    w_tot = msg(3)

  else
    gam   = 0.0
    gam2  = 0.0
    w_tot = 0.0
  endif

  ! calculate the charge-weighted first and second moments, and the total weight
  ! of the local stage
  do i = 1, this%npp

    w = abs( this%q(i) )
    gam_tmp = sqrt( 1.0 + this%p(1,i)*this%p(1,i) + this%p(2,i)*this%p(2,i) &
      + this%p(3,i)*this%p(3,i) )
    gam     = gam + gam_tmp * w
    gam2    = gam2 + gam_tmp * gam_tmp * w
    w_tot   = w_tot + w

  enddo

  ! send message to next stage
  ! the last stage calculates the results
  msg = (/gam, gam2, w_tot/)

  if ( idproc_dest < num_procs() ) then
    call mpi_isend( msg, 3, p_dtype_real, idproc_dest, recv_tag, comm_world(), &
      id, ierr )
  else

    call mpi_reduce( msg, rslt, 3, p_dtype_real, MPI_SUM, id_proc_loc(), comm_loc(), ierr )

    ! the first processor in the last stage output data
    if ( id_proc_loc() == 0 ) then

      rslt(1:2) = rslt(1:2) / rslt(3)
      ene_spread = sqrt( rslt(2) - rslt(1) * rslt(1) ) / rslt(1)
    
      write (fid,'(I8.8,D22.14)') tstep, ene_spread
      flush(fid)

    endif

  endif

end subroutine get_ene_spread_part3d_popas

end module part3d_popas_class
