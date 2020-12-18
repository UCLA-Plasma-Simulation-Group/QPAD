program test_field_ez

use field_e_class
use field_class
use field_src_class
use system
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( field_e ) :: e
type( field_jay ) :: jay

integer :: num_modes = 2, dim = 3, order = p_fs_2order, part_shape = p_ps_linear
integer, dimension(2) :: nd = (/128, 1/), nvp = (/1, 1/)
integer, dimension(2,2) :: gc_num
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: ujay_re => null(), ujay_im => null()
type( ufield ), dimension(:), pointer :: ue_re => null(), ue_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode

call MPI_INIT( ierr )
call init_errors( 2, 3 )

call write_dbg( 'main', 'test_field_ez', 0, 'starts' )

dr = 1.0 / (nd(1)-0.5)
dxi = 1.0

gc_num(:,1) = (/0,0/)
gc_num(:,2) = (/0,0/)

call jay%new( num_modes, dr, dxi, nd, nvp, part_shape )
call e%new( num_modes, dr, dxi, nd, nvp, part_shape, p_entity_plasma )


ujay_re => jay%get_rf_re()
ujay_im => jay%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nd(1)
    r = (i-0.5)*dr
    ujay_re(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
    ujay_re(mode)%f1(2,i) = exp( -((r-0.5)/0.05)**2 )
    ujay_re(mode)%f1(3,i) = exp( -((r-0.5)/0.05)**2 )
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nd(1)
    r = (i-0.5)*dr
    ujay_im(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
    ujay_im(mode)%f1(2,i) = exp( -((r-0.5)/0.05)**2 )
    ujay_im(mode)%f1(3,i) = exp( -((r-0.5)/0.05)**2 )
  enddo

enddo

! call HYPRE_StructMatrixGetBoxValues( psi%solver(0)%A, 1, nd(1), 3, &
!   psi%solver(0)%stencil_idx, values, ierr )

! print *, "A = ", values(1:nd(1)*3)

! call HYPRE_StructVectorGetBoxValues( psi%solver(0)%b )

call e%solve(jay)

ue_re => e%get_rf_re()
ue_im => e%get_rf_im()

p => ue_re(0)%get_f1()
call write_data( p, 'ez-re-0.txt', 3 )
p => ue_re(1)%get_f1()
call write_data( p, 'ez-re-1.txt', 3 )
p => ue_re(2)%get_f1()
call write_data( p, 'ez-re-2.txt', 3 )
p => ue_im(1)%get_f1()
call write_data( p, 'ez-im-1.txt', 3 )
p => ue_im(2)%get_f1()
call write_data( p, 'ez-im-2.txt', 3 )



call jay%del()
call e%del()

call write_dbg( 'main', 'test_field_ez', 0, 'ends' )

call end_errors()
call MPI_FINALIZE( ierr )

end program test_field_ez