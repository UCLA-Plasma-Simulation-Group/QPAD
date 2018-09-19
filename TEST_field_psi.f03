program test_field_psi

use field_psi_class
use field_class
use system
use param
use mpi

implicit none

type( field_psi ) :: psi
type( field ) :: q

integer :: num_modes = 2, dim = 1, order = p_fs_2order, part_shape = p_ps_linear
integer, dimension(2) :: nd = (/128, 1/), nvp = (/1, 1/)
integer, dimension(2,2) :: gc_num
real :: dr = 1.0, dxi = 1.0

integer :: ierr, i

call MPI_INIT( ierr )
call init_errors( 2, 3 )

call write_dbg( 'main', 'test_field_psi', 0, 'starts' )

gc_num(:,1) = (/0,0/)
gc_num(:,2) = (/0,0/)

call q%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )
call psi%new( num_modes, dr, dxi, nd, nvp, order, part_shape )

! do something

call q%del()
call psi%del()


call write_dbg( 'main', 'test_field_psi', 0, 'ends' )

call end_errors()
call MPI_FINALIZE( ierr )

end program test_field_psi