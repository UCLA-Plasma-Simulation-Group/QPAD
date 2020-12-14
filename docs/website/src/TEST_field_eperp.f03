program test_field_eperp

use field_b_class
use field_e_class
use field_psi_class
use field_class
use field_src_class
use system
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( field_b ) :: b
type( field_e ) :: e
type( field_psi ) :: psi
type( field_rho ) :: rho
type( field_jay ) :: jay
type( field_djdxi ) :: djdxi

integer :: num_modes = 3, dim = 3, order = p_fs_2order, part_shape = p_ps_linear
integer, dimension(2) :: nd = (/128, 1/), nvp = (/1, 1/)
integer, dimension(2,2) :: gc_num
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: uq_re => null(), uq_im => null()
type( ufield ), dimension(:), pointer :: ue_re => null(), ue_im => null()
type( ufield ), dimension(:), pointer :: ujay_re => null(), ujay_im => null()
type( ufield ), dimension(:), pointer :: udjdxi_re => null(), udjdxi_im => null()
type( ufield ), dimension(:), pointer :: ub_re => null(), ub_im => null()
type( ufield ), dimension(:), pointer :: upsi_re => null(), upsi_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode

call MPI_INIT( ierr )
call init_errors( 2, 3 )

call write_dbg( 'main', 'test_field_eperp', 0, 'starts' )

dr = 1.0 / (nd(1)-0.5)
dxi = 1.0

gc_num(:,1) = (/0,0/)
gc_num(:,2) = (/0,0/)

call rho%new( num_modes, dr, dxi, nd, nvp, part_shape )
call b%new( num_modes, dr, dxi, nd, nvp, part_shape, p_entity_plasma )
call e%new( num_modes, dr, dxi, nd, nvp, part_shape, p_entity_plasma )
call psi%new( num_modes, dr, dxi, nd, nvp, part_shape )
call jay%new( num_modes, dr, dxi, nd, nvp, part_shape )
call djdxi%new( num_modes, dr, dxi, nd, nvp, part_shape )

! solve psi
uq_re => rho%get_rf_re()
uq_im => rho%get_rf_im()

do mode = 0, num_modes
  
  do i = 1, nd(1)
    r = (i-0.5)*dr
    uq_re(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nd(1)
    r = (i-0.5)*dr
    uq_im(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
  enddo

enddo

call psi%solve(rho)

! solve b
ujay_re => jay%get_rf_re()
ujay_im => jay%get_rf_im()
udjdxi_re => djdxi%get_rf_re()
udjdxi_im => djdxi%get_rf_im()

do mode = 0, num_modes
  
  do i = 1, nd(1)
    r = (i-0.5)*dr
    ujay_re(mode)%f1(1,i) = 0.0
    ujay_re(mode)%f1(2,i) = 0.0
    ujay_re(mode)%f1(3,i) = 0.0
    udjdxi_re(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_re(mode)%f1(2,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_re(mode)%f1(3,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nd(1)
    r = (i-0.5)*dr
    ujay_im(mode)%f1(1,i) = 0.0
    ujay_im(mode)%f1(2,i) = 0.0
    ujay_im(mode)%f1(3,i) = 0.0
    udjdxi_im(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_im(mode)%f1(2,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_im(mode)%f1(3,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
  enddo

enddo

call b%solve( djdxi, jay )

! solve electric field
call e%solve( b, psi )

ue_re => e%get_rf_re()
ue_im => e%get_rf_im()

ub_re => b%get_rf_re()
ub_im => b%get_rf_im()

upsi_re => psi%get_rf_re()
upsi_im => psi%get_rf_im()

p => ue_re(0)%get_f1()
call write_data( p, 'er-re-0.txt', 1 )
call write_data( p, 'ephi-re-0.txt', 2 )
p => ue_re(1)%get_f1()
call write_data( p, 'er-re-1.txt', 1 )
call write_data( p, 'ephi-re-1.txt', 2 )
p => ue_re(2)%get_f1()
call write_data( p, 'er-re-2.txt', 1 )
call write_data( p, 'ephi-re-2.txt', 2 )
p => ue_im(1)%get_f1()
call write_data( p, 'er-im-1.txt', 1 )
call write_data( p, 'ephi-im-1.txt', 2 )
p => ue_im(2)%get_f1()
call write_data( p, 'er-im-2.txt', 1 )
call write_data( p, 'ephi-im-2.txt', 2 )

p => ub_re(0)%get_f1()
call write_data( p, 'br-re-0.txt', 1 )
call write_data( p, 'bphi-re-0.txt', 2 )
p => ub_re(1)%get_f1()
call write_data( p, 'br-re-1.txt', 1 )
call write_data( p, 'bphi-re-1.txt', 2 )
p => ub_re(2)%get_f1()
call write_data( p, 'br-re-2.txt', 1 )
call write_data( p, 'bphi-re-2.txt', 2 )
p => ub_im(1)%get_f1()
call write_data( p, 'br-im-1.txt', 1 )
call write_data( p, 'bphi-im-1.txt', 2 )
p => ub_im(2)%get_f1()
call write_data( p, 'br-im-2.txt', 1 )
call write_data( p, 'bphi-im-2.txt', 2 )

p => upsi_re(0)%get_f1()
call write_data( p, 'psi-re-0.txt', 1 )
p => upsi_re(1)%get_f1()
call write_data( p, 'psi-re-1.txt', 1 )
p => upsi_re(2)%get_f1()
call write_data( p, 'psi-re-2.txt', 1 )
p => upsi_im(1)%get_f1()
call write_data( p, 'psi-im-1.txt', 1 )
p => upsi_im(2)%get_f1()
call write_data( p, 'psi-im-2.txt', 1 )

call rho%del()
call psi%del()
call jay%del()
call djdxi%del()
call b%del()
call e%del()

call write_dbg( 'main', 'test_field_eperp', 0, 'ends' )

call end_errors()
call MPI_FINALIZE( ierr )

end program test_field_eperp