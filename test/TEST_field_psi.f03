program test_field_psi

use field_psi_class
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

type( field_psi ) :: psi
type( field_rho ) :: q

integer :: num_modes = 4, part_shape = p_ps_linear
integer :: nr = 512, nz = 1, nrp, noff
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: uq_re => null(), uq_im => null()
type( ufield ), dimension(:), pointer :: upsi_re => null(), upsi_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode
character(len=32) :: filename

allocate( pp, gp )
call pp%new(nst=1)

call init_stdout( pp%getidproc() )
call init_errors( eunit=2, idproc=pp%getlidproc(), monitor=3 )
call write_dbg( 'main', 'test_field_psi', 0, 'starts' )

call gp%new( pp, nr, nz )

dr = 1.0 / (nr-0.5)
dxi = 1.0
nrp = gp%get_ndp(1)
noff = gp%get_noff(1)

call q%new( pp, gp, dr, dxi, num_modes, part_shape )
call psi%new( pp, gp, dr, dxi, num_modes, part_shape, iter_tol=1.0d-6 )

uq_re => q%get_rf_re()
uq_im => q%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nrp 
    r = (real(i+noff)-0.5)*dr
    uq_re(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nrp
    r = (i+noff-0.5)*dr
    uq_im(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
  enddo

enddo

call psi%solve(q)

upsi_re => psi%get_rf_re()
upsi_im => psi%get_rf_im()


! output
p => upsi_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'psi-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )

do i = 1, num_modes
  p => upsi_re(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'psi-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  p => upsi_im(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'psi-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
enddo

call q%del()
call psi%del()

call write_dbg( 'main', 'test_field_psi', 0, 'ends' )

call end_errors()
call gp%del()
call pp%del()

end program test_field_psi