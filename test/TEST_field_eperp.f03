program test_field_eperp

use parallel_pipe_class
use grid_class
use field_b_class
use field_e_class
use field_psi_class
use field_src_class
use sys
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( parallel_pipe ), pointer :: pp => null()
type( grid ), pointer :: gp => null()

type( field_b ) :: b
type( field_e ) :: e
type( field_psi ) :: psi
type( field_rho ) :: rho
type( field_jay ) :: jay
type( field_djdxi ) :: djdxi

integer :: num_modes = 2, part_shape = p_ps_linear
integer :: nr = 128, nz = 1, nrp, noff, num_iter = 3
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: uq_re => null(), uq_im => null()
type( ufield ), dimension(:), pointer :: ue_re => null(), ue_im => null()
type( ufield ), dimension(:), pointer :: ujay_re => null(), ujay_im => null()
type( ufield ), dimension(:), pointer :: udjdxi_re => null(), udjdxi_im => null()
type( ufield ), dimension(:), pointer :: ub_re => null(), ub_im => null()
type( ufield ), dimension(:), pointer :: upsi_re => null(), upsi_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode, iter
character(len=32) :: filename

allocate( pp, gp )
call pp%new(nst=1)

call init_stdout( pp%getidproc() )
call init_errors( eunit=2, idproc=pp%getlidproc(), monitor=3 )
call write_dbg( 'main', 'test_field_eperp', 0, 'starts' )

call gp%new( pp, nr, nz )

dr = 1.0 / (nr-0.5)
dxi = 1.0
nrp = gp%get_ndp(1)
noff = gp%get_noff(1)

call rho%new( pp, gp, dr, dxi, num_modes, part_shape )
call b%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_plasma, iter_tol=1.0d-6 )
call e%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_plasma, iter_tol=1.0d-6 )
call psi%new( pp, gp, dr, dxi, num_modes, part_shape, iter_tol=1.0d-6 )
call jay%new( pp, gp, dr, dxi, num_modes, part_shape )
call djdxi%new( pp, gp, dr, dxi, num_modes, part_shape )

! solve psi
uq_re => rho%get_rf_re()
uq_im => rho%get_rf_im()

do mode = 0, num_modes
  
  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    uq_re(mode)%f1(1,i) = -0.5*exp( -(r/0.1)**2 )+1.0
    ! uq_re(mode)%f1(1,i) = 0.0
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    uq_im(mode)%f1(1,i) = -0.5*exp( -(r/0.1)**2 )+1.0
    ! uq_im(mode)%f1(1,i) = 0.0
  enddo

enddo

call psi%solve(rho)

! solve b
ujay_re => jay%get_rf_re()
ujay_im => jay%get_rf_im()
udjdxi_re => djdxi%get_rf_re()
udjdxi_im => djdxi%get_rf_im()

do mode = 0, num_modes
  
  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    ujay_re(mode)%f1(1,i) = 0.0
    ujay_re(mode)%f1(2,i) = 0.0
    ujay_re(mode)%f1(3,i) = -0.5*exp( -(r/0.1)**2 )
    udjdxi_re(mode)%f1(1,i) = 0.0
    udjdxi_re(mode)%f1(2,i) = 0.0
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nrp
    r = (real(i+noff)-0.5)*dr
    ujay_im(mode)%f1(1,i) = 0.0
    ujay_im(mode)%f1(2,i) = 0.0
    ujay_im(mode)%f1(3,i) = 0.5*exp( -(r/0.1)**2 )
    udjdxi_im(mode)%f1(1,i) = 0.0
    udjdxi_im(mode)%f1(2,i) = 0.0
  enddo

enddo

do iter = 1, num_iter

  print *, 'iter = ', iter

  call b%solve( djdxi, jay )
  ! ub_re_new => b%get_rf_re()
  ! ub_im_new => b%get_rf_im()

  ! do i = 0, num_modes
  !   ub_res = ub_re_new(i) - ub_re_old(i)
  !   ub_re_old(i) = ub_re_new(i)
  !   p => ub_res%get_f1()
  !   if (i==0) then
  !     res(iter,1) = maxval(abs(p(1,:)))
  !     res(iter,2) = maxval(abs(p(2,:)))
  !     cycle
  !   endif
  !   res(iter,4*i-1) = maxval(abs(p(1,:)))
  !   res(iter,4*i)   = maxval(abs(p(2,:)))

  !   ub_res = ub_im_new(i) - ub_im_old(i)
  !   ub_im_old(i) = ub_im_new(i)
  !   p => ub_res%get_f1()
  !   res(iter,4*i+1) = maxval(abs(p(1,:)))
  !   res(iter,4*i+2) = maxval(abs(p(2,:)))
  ! enddo

enddo

! solve electric field
call e%solve( b, psi )

ue_re => e%get_rf_re()
ue_im => e%get_rf_im()

ub_re => b%get_rf_re()
ub_im => b%get_rf_im()

upsi_re => psi%get_rf_re()
upsi_im => psi%get_rf_im()

! output e-field
p => ue_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'er-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )
write( filename, '(A,I0.3,A)' ) 'ephi-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 2 )
do i = 1, num_modes
  p => ue_re(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'er-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'ephi-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 2 )
  p => ue_im(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'er-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'ephi-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 2 )
enddo

! output b-field
p => ub_re(0)%get_f1()
write( filename, '(A,I0.3,A)' ) 'br-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 1 )
write( filename, '(A,I0.3,A)' ) 'bphi-re-0-', pp%getlidproc(), '.txt'
call write_data( p, trim(filename), 2 )
do i = 1, num_modes
  p => ub_re(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'br-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'bphi-re-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 2 )
  p => ub_im(i)%get_f1()
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'br-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 1 )
  write( filename, '(A,I0.1,A,I0.3,A)' ) 'bphi-im-', i, '-', pp%getlidproc(), '.txt'
  call write_data( p, trim(filename), 2 )
enddo

! output psi-field
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

call rho%del()
call psi%del()
call jay%del()
call djdxi%del()
call b%del()
call e%del()

call write_dbg( 'main', 'test_field_eperp', 0, 'ends' )

call end_errors()
call gp%del()
call pp%del()

end program test_field_eperp