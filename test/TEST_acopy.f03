program test_acopy

use field_src_class

use parallel_pipe_class
use grid_class
use sys
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( parallel_pipe ), pointer :: pp => null()
type( grid ), pointer :: gp => null()

type( field_rho ) :: q

integer :: num_modes = 2, part_shape = p_ps_linear
integer :: nr = 32, nz = 1, nrp, nzp, noff1, noff2
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: uq_re => null(), uq_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode
character(len=32) :: filename

allocate( pp, gp )
call pp%new(nst=1)

call init_stdout( pp%getidproc() )
call init_errors( eunit=2, idproc=pp%getlidproc(), monitor=3 )
call write_dbg( 'main', 'test_acopy', 0, 'starts' )

call gp%new( pp, nr, nz )

dr = 1.0 / (nr-0.5)
dxi = 1.0
nrp = gp%get_ndp(1)
nzp = gp%get_ndp(2)
noff1 = gp%get_noff(1)
noff2 = gp%get_noff(2)

call q%new( pp, gp, dr, dxi, num_modes, part_shape )

uq_re => q%get_rf_re()
uq_im => q%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nrp 
    r = (real(i+noff1)-0.5)*dr
    ! uq_re(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    uq_re(mode)%f1(1,i) = 1.0
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nrp
    r = (i+noff1-0.5)*dr
    ! uq_im(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    uq_im(mode)%f1(1,i) = 1.0
  enddo

enddo

call q%acopy_gc_f1()

! output
p => uq_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'q-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )
p => uq_re(1)%get_f1()
write( filename, '(A,I0.3,A)' ) 'q-re-1-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )
p => uq_re(2)%get_f1()
write( filename, '(A,I0.3,A)' ) 'q-re-2-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )
p => uq_im(1)%get_f1()
write( filename, '(A,I0.3,A)' ) 'q-im-1-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )
p => uq_im(2)%get_f1()
write( filename, '(A,I0.3,A)' ) 'q-im-2-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )

call q%del()

call write_dbg( 'main', 'test_acopy', 0, 'ends' )

call end_errors()
call gp%del()
call pp%del()

end program test_acopy