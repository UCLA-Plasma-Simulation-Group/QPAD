program test_field_bt_iter_fast

use parallel_pipe_class
use grid_class
use field_b_class
use field_src_class
use sys
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( parallel_pipe ), pointer :: pp => null()
type( grid ), pointer :: gp => null()

type( field_b ) :: b, b_old
type( field_jay ) :: jay
type( field_djdxi ) :: djdxi

integer :: num_modes = 3, part_shape = p_ps_linear
integer :: nr = 128, nz = 1, nrp, noff
real :: dr, dxi, r
integer, dimension(2,2) :: gc_num

type( ufield ), dimension(:), pointer :: ujay_re => null(), ujay_im => null()
type( ufield ), dimension(:), pointer :: udjdxi_re => null(), udjdxi_im => null()
type( ufield ), dimension(:), pointer :: ub_re => null(), ub_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode, iter
integer :: num_iter = 1
real, dimension(:,:), pointer :: res => null()
character(len=32) :: filename

allocate( pp, gp )
call pp%new(nst=1)

call init_stdout( pp%getidproc() )
call init_errors( eunit=2, idproc=pp%getlidproc(), monitor=3 )
call write_dbg( 'main', 'test_field_bt_iter_fast', 0, 'starts' )

call gp%new( pp, nr, nz )

dr = 1.0 / (nr-0.5)
dxi = 1.0
nrp = gp%get_ndp(1)
noff = gp%get_noff(1)

call jay%new( pp, gp, dr, dxi, num_modes, part_shape )
call djdxi%new( pp, gp, dr, dxi, num_modes, part_shape )
call b%new( pp, gp, dr, dxi, num_modes, part_shape, &
  boundary=p_bnd_open, entity=p_entity_plasma, iter_tol=1.0d-6 )
call b_old%new( pp, gp, dr, dxi, num_modes, part_shape, &
  boundary=p_bnd_open, entity=p_entity_plasma_old, iter_tol=1.0d-6 )

ujay_re => jay%get_rf_re()
ujay_im => jay%get_rf_im()
udjdxi_re => djdxi%get_rf_re()
udjdxi_im => djdxi%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    ujay_re(mode)%f1(1,i) = 0.0
    ujay_re(mode)%f1(2,i) = 0.0
    ujay_re(mode)%f1(3,i) = 0.1*exp( -((r-0.5)/0.05)**2 )
    ! ujay_re(mode)%f1(3,i) = 0.0
    udjdxi_re(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_re(mode)%f1(2,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    ! udjdxi_re(mode)%f1(1,i) = 0.0
    ! udjdxi_re(mode)%f1(2,i) = 0.0
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    ujay_im(mode)%f1(1,i) = 0.0
    ujay_im(mode)%f1(2,i) = 0.0
    ujay_im(mode)%f1(3,i) = 0.1*exp( -((r-0.5)/0.05)**2 )
    ! ujay_im(mode)%f1(3,i) = 0.0
    udjdxi_im(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_im(mode)%f1(2,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    ! udjdxi_im(mode)%f1(1,i) = 0.0
    ! udjdxi_im(mode)%f1(2,i) = 0.0
  enddo

enddo

call jay%copy_gc_f1( bnd_ax = .false. )
do i = 1, num_iter
  call b%solve( djdxi, jay )
  call b_old%solve_old( djdxi, jay )
enddo

ub_re => b%get_rf_re()
ub_im => b%get_rf_im()

p => ub_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'b1-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )
write( filename, '(A,I0.3,A)' ) 'b2-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 2 )
do i = 1, num_modes
  p => ub_re(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'b1-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'b2-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 2 )
  p => ub_im(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'b1-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'b2-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 2 )
enddo

ub_re => b_old%get_rf_re()
ub_im => b_old%get_rf_im()

p => ub_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'b1_old-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )
write( filename, '(A,I0.3,A)' ) 'b2_old-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 2 )
do i = 1, num_modes
  p => ub_re(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'b1_old-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'b2_old-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 2 )
  p => ub_im(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'b1_old-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'b2_old-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 2 )
enddo

call jay%del()
call b%del()
call b_old%del()

! deallocate( ub_re, ub_im )

call write_dbg( 'main', 'test_field_bt_iter_fast', 0, 'ends' )

call end_errors()
call gp%del()
call pp%del()

end program test_field_bt_iter_fast