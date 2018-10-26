program test_ufield

use parallel_pipe_class
use grid_class
use ufield_class

implicit none

type( parallel_pipe ), pointer :: pp => null()
type( grid ), pointer :: gp => null()
type( ufield ) :: a1, a2, a3

integer :: dim = 3
integer :: nr = 5, nz = 1
integer, dimension(2,2) :: gc_num = 0
real :: k = 5.0

real :: start, finish
integer :: i

allocate( pp, gp )
call pp%new(nst=1)
call gp%new( pp, nr, nz )


call a1%new( pp, gp, dim, gc_num )
call a2%new( pp, gp, dim, gc_num )
call a3%new( pp, gp, dim, gc_num )

a1%f1 = 1.0
a2%f1 = 2.0

print *, "a1 = ", a1%f1(1,:)
print *, "a2 = ", a2%f1(1,:)

a3 = k*(a1+a2)+(a1+k)*a2

print *, "a3 = ", a3%f1(1,:)

call pp%del()
call gp%del()

end program test_ufield