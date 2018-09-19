program test_field_psi

use field_psi_class
use field_class
use system
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( field_psi ) :: psi
type( field ) :: q

integer :: num_modes = 2, dim = 1, order = p_fs_2order, part_shape = p_ps_linear
integer, dimension(2) :: nd = (/128, 1/), nvp = (/1, 1/)
integer, dimension(2,2) :: gc_num
real :: dr = 1.0, dxi = 1.0

type( ufield ), dimension(:), pointer :: uq_re => null(), uq_im => null()
type( ufield ), dimension(:), pointer :: upsi_re => null(), upsi_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode
real, dimension(1024) :: values = 0.0

call MPI_INIT( ierr )
call init_errors( 2, 3 )

call write_dbg( 'main', 'test_field_psi', 0, 'starts' )

gc_num(:,1) = (/0,0/)
gc_num(:,2) = (/0,0/)

call q%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )
call psi%new( num_modes, dr, dxi, nd, nvp, order, part_shape )


uq_re => q%get_rf_re()
uq_im => q%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nd(1)/8
    uq_re(mode)%f1(1,i) = 1.0
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nd(1)/8
    uq_im(mode)%f1(1,i) = 1.0
  enddo

enddo

! call HYPRE_StructMatrixGetBoxValues( psi%solver(0)%A, 1, nd(1), 3, &
!   psi%solver(0)%stencil_idx, values, ierr )

! print *, "A = ", values(1:nd(1)*3)

! call HYPRE_StructVectorGetBoxValues( psi%solver(0)%b )

call psi%solve(q)

upsi_re => psi%get_rf_re()
upsi_im => psi%get_rf_im()

p => upsi_re(0)%get_f1()
call write_data( p, 'psi-re-0.txt', 1 )
p => upsi_re(1)%get_f1()
call write_data( p, 'psi-re-1.txt', 1 )
p => upsi_re(2)%get_f1()
call write_data( p, 'psi-re-2.txt', 1 )



call q%del()
call psi%del()

call write_dbg( 'main', 'test_field_psi', 0, 'ends' )

call end_errors()
call MPI_FINALIZE( ierr )

end program test_field_psi