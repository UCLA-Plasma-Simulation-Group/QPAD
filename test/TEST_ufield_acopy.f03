program test_ufield_arith

use parallel_pipe_class
use grid_class
use ufield_class
use param
use sysutil

use debug_tool

implicit none

type( parallel_pipe ), pointer :: pp => null()
type( grid ), pointer :: gp => null()
type( ufield ) :: fld

integer :: dim = 1, mode = 0
integer :: nr = 100, nz = 1, nrp, noff
integer, dimension(2,2) :: gc_num = 0

integer :: idproc, i

allocate( pp, gp )
call pp%new(nst=1)
call gp%new( pp, nr, nz, dr=1.0, dxi=1.0 )
idproc = pp%getlidproc()

gc_num(:,1) = (/10, 10/)
gc_num(:,2) = (/0, 0/)
call fld%new( pp, gp, dim, mode, gc_num, has_2d=.true. )

nrp = fld%get_ndp(1)
noff = fld%get_noff(1)
do i = 1, nrp
    fld%f1(1,i) = real(noff+i)
enddo

call fld%copy_gc_f1()
call fld%acopy_gc_f1( dir=p_mpi_bothway, nc=10 )

call write_data( fld%f1, 'fld_'//num2str(idproc)//'.txt', 1 )

call pp%del()
call gp%del()

end program test_ufield_arith