program test_smooth

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

type( field_rho ) :: q

integer :: num_modes = 1, part_shape = p_ps_linear
integer :: smooth_type = p_smooth_compensated, smooth_order = 5
integer :: nr = 512, nz = 1, nrp, noff
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: uq_re => null(), uq_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode
character(len=32) :: filename

allocate( pp, gp )
call pp%new(nst=1)

call init_stdout( pp%getidproc() )
call init_errors( eunit=2, idproc=pp%getlidproc(), monitor=3 )
call write_dbg( 'main', 'test_smooth', 0, 'starts' )

call gp%new( pp, nr, nz )

dr = 1.0 / (nr-0.5)
dxi = 1.0
nrp = gp%get_ndp(1)
noff = gp%get_noff(1)

call q%new( pp, gp, dr, dxi, num_modes, part_shape, &
  smooth_type, smooth_order )

uq_re => q%get_rf_re()
uq_im => q%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nrp 
    r = (real(i+noff)-0.5)*dr
    uq_re(mode)%f1(1,i) = 0.98 + exp( -(r/0.2)**2 ) + 0.04*rand()
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nrp
    r = (i+noff-0.5)*dr
    uq_im(mode)%f1(1,i) = 0.98 + exp( -(r/0.2)**2 ) + 0.04*rand()
  enddo

enddo

! output the noisy data
p => uq_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'q-noise-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )

do i = 1, num_modes
  p => uq_re(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'q-noise-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  p => uq_im(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'q-noise-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
enddo

call q%copy_gc_f1()
call q%smooth_f1() ! make sure the guard cells have been copied before smooth

! output the smoothed data
p => uq_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'q-smooth-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )

do i = 1, num_modes
  p => uq_re(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'q-smooth-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  p => uq_im(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'q-smooth-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
enddo

call q%del()

call write_dbg( 'main', 'test_smooth', 0, 'ends' )

call end_errors()
call gp%del()
call pp%del()

end program test_smooth