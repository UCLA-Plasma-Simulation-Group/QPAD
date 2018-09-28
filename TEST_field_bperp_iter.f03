program test_field_bperp_iter

use field_b_class
use field_class
use system
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( field_b ) :: b
type( field ) :: jay, djdxi

integer :: num_modes = 2, dim = 3, order = p_fs_2order, part_shape = p_ps_linear
integer, dimension(2) :: nd = (/128, 1/), nvp = (/1, 1/)
integer, dimension(2,2) :: gc_num
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: ujay_re => null(), ujay_im => null()
type( ufield ), dimension(:), pointer :: udjdxi_re => null(), udjdxi_im => null()
type( ufield ), dimension(:), pointer :: ub_re => null(), ub_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode

call MPI_INIT( ierr )
call init_errors( 2, 3 )

call write_dbg( 'main', 'test_field_bperp_iter', 0, 'starts' )

dr = 1.0 / (nd(1)-0.5)
dxi = 1.0

gc_num(:,1) = (/0,0/)
gc_num(:,2) = (/0,0/)

call jay%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )
call djdxi%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )
call b%new( num_modes, dr, dxi, nd, nvp, order, part_shape )


ujay_re => jay%get_rf_re()
ujay_im => jay%get_rf_im()
udjdxi_re => djdxi%get_rf_re()
udjdxi_im => djdxi%get_rf_im()

! set the charge
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

! call HYPRE_StructMatrixGetBoxValues( psi%solver(0)%A, 1, nd(1), 3, &
!   psi%solver(0)%stencil_idx, values, ierr )

! print *, "A = ", values(1:nd(1)*3)

! call HYPRE_StructVectorGetBoxValues( psi%solver(0)%b )

call b%solve( djdxi, jay )

ub_re => b%get_rf_re()
ub_im => b%get_rf_im()

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



call jay%del()
call djdxi%del()
call b%del()

call write_dbg( 'main', 'test_field_bperp_iter', 0, 'ends' )

call end_errors()
call MPI_FINALIZE( ierr )

end program test_field_bperp_iter