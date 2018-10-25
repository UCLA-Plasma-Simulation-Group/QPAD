program test_field_bperp

use field_b_class
use field_src_class

use parallel_pipe_class
use grid_class
use system
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( parallel_pipe ), pointer :: pp => null()
type( grid ), pointer :: gp => null()

type( field_b ) :: b
type( field_rho ) :: rho

integer :: num_modes = 2, part_shape = p_ps_linear
integer :: nr = 128, nz = 1, nrp, noff
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: uq_re => null(), uq_im => null()
type( ufield ), dimension(:), pointer :: ub_re => null(), ub_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode
character(len=32) :: filename

allocate( pp, gp )
call pp%new(nst=1)

call init_stdout( pp%getidproc() )
call init_errors( eunit=2, idproc=pp%getlidproc(), monitor=3 )
call write_dbg( 'main', 'test_field_bperp', 0, 'starts' )

call gp%new( pp, nr, nz )

dr = 1.0 / (nr-0.5)
dxi = 1.0
nrp = gp%get_ndp(1)
noff = gp%get_noff(1)

call rho%new( pp, gp, dr, dxi, num_modes, part_shape )
call b%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_beam )


uq_re => rho%get_rf_re()
uq_im => rho%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    uq_re(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    uq_im(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
  enddo

enddo

call b%solve(rho)

ub_re => b%get_rf_re()
ub_im => b%get_rf_im()

p => ub_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'br-re-0-', pp%getlidproc(), '.txt'
call write_data( p, filename, 1 )
write( filename, '(A,I0.3,A)' ) 'bphi-re-0-', pp%getlidproc(), '.txt'
call write_data( p, filename, 2 )
p => ub_re(1)%get_f1()
write( filename, '(A,I0.3,A)' ) 'br-re-1-', pp%getlidproc(), '.txt'
call write_data( p, filename, 1 )
write( filename, '(A,I0.3,A)' ) 'bphi-re-1-', pp%getlidproc(), '.txt'
call write_data( p, filename, 2 )
p => ub_re(2)%get_f1()
write( filename, '(A,I0.3,A)' ) 'br-re-2-', pp%getlidproc(), '.txt'
call write_data( p, filename, 1 )
write( filename, '(A,I0.3,A)' ) 'bphi-re-2-', pp%getlidproc(), '.txt'
call write_data( p, filename, 2 )
p => ub_im(1)%get_f1()
write( filename, '(A,I0.3,A)' ) 'br-im-1-', pp%getlidproc(), '.txt'
call write_data( p, filename, 1 )
write( filename, '(A,I0.3,A)' ) 'bphi-im-1-', pp%getlidproc(), '.txt'
call write_data( p, filename, 2 )
p => ub_im(2)%get_f1()
write( filename, '(A,I0.3,A)' ) 'br-im-2-', pp%getlidproc(), '.txt'
call write_data( p, filename, 1 )
write( filename, '(A,I0.3,A)' ) 'bphi-im-2-', pp%getlidproc(), '.txt'
call write_data( p, filename, 2 )

call rho%del()
call b%del()

call write_dbg( 'main', 'test_field_bperp', 0, 'ends' )

call end_errors()
call gp%del()
call pp%del()

end program test_field_bperp