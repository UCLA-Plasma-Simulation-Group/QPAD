program test_field_ez

use parallel_pipe_class
use grid_class
use field_e_class
use field_class
use field_src_class
use sys
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( parallel_pipe ), pointer :: pp => null()
type( grid ), pointer :: gp => null()

type( field_e ) :: e
type( field_jay ) :: jay

integer :: num_modes = 3, part_shape = p_ps_linear
integer :: nr = 128, nz = 1, nrp, noff
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: ujay_re => null(), ujay_im => null()
type( ufield ), dimension(:), pointer :: ue_re => null(), ue_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode
character(len=32) :: filename

allocate( pp, gp )
call pp%new(nst=1)

call init_stdout( pp%getidproc() )
call init_errors( eunit=2, idproc=pp%getlidproc(), monitor=3 )
call write_dbg( 'main', 'test_field_ez', 0, 'starts' )

call gp%new( pp, nr, nz )

dr = 1.0 / (nr-0.5)
dxi = 1.0
nrp = gp%get_ndp(1)
noff = gp%get_noff(1)

call jay%new( pp, gp, dr, dxi, num_modes, part_shape )
call e%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_plasma, iter_tol=1.0d-6 )


ujay_re => jay%get_rf_re()
ujay_im => jay%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    ujay_re(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
    ujay_re(mode)%f1(2,i) = exp( -((r-0.5)/0.05)**2 )
    ujay_re(mode)%f1(3,i) = exp( -((r-0.5)/0.05)**2 )
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    ujay_im(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
    ujay_im(mode)%f1(2,i) = exp( -((r-0.5)/0.05)**2 )
    ujay_im(mode)%f1(3,i) = exp( -((r-0.5)/0.05)**2 )
  enddo

enddo

call e%solve(jay)

ue_re => e%get_rf_re()
ue_im => e%get_rf_im()

p => ue_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'ez-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 3 )
do i = 1, num_modes
  p => ue_re(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'ez-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 3 )
  p => ue_im(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'ez-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 3 )
enddo

call jay%del()
call e%del()

call write_dbg( 'main', 'test_field_ez', 0, 'ends' )

call end_errors()
call gp%del()
call pp%del()

end program test_field_ez