program test_ufield_arith

use parallel_pipe_class
use grid_class
use ufield_class

implicit none

type( parallel_pipe ), pointer :: pp => null()
type( grid ), pointer :: gp => null()
type( ufield ) :: a1, a2, a3

integer :: dim = 3, mode = 0
integer :: nr = 10, nz = 2
integer, dimension(2,2) :: gc_num = 0
real :: k = 5.0

real :: start, finish
integer :: i

allocate( pp, gp )
call pp%new(nst=1)
call gp%new( pp, nr, nz )


call a1%new( pp, gp, dim, mode, gc_num, has_2d=.true. )
call a2%new( pp, gp, dim, mode, gc_num, has_2d=.true. )
call a3%new( pp, gp, dim, mode, gc_num, has_2d=.true. )

a1%f1 = 0.5
a2%f1 = 0.1

a1%f2 = 1.0
a2%f2 = 2.0

! print *, "a1_f1(dim=1) = ", a1%f1(1,:)
! print *, "a2_f1(dim=1) = ", a2%f1(1,:)

! print *, "a1_f2(dim=1) = ", a1%f2(1,:,:)
! print *, "a2_f2(dim=1) = ", a2%f2(1,:,:)

! ! arithmetic on f1
! a3 = k * ( a1 + a2 ) - a1
! ! arithmetic on f2
! call a3%as( k .dot. (a1 .add. a2) .sub. a1  )

! print *, "a3_f1(dim=1) = ", a3%f1(1,:)
! print *, "a3_f2(dim=1) = ", a3%f2(1,:,:)

! call add_f1( a1, a2, a3, (/1,2/), (/2,1/), (/1,2/) )
! call add_f2( a1, a2, a3, (/1,2/), (/2,1/), (/1,2/) )

! print *, "a3_f1(dim=1) = ", a3%f1(1,:)
! print *, "a3_f2(dim=1) = ", a3%f2(1,:,:)
! print *, "a3_f1(dim=2) = ", a3%f1(2,:)
! print *, "a3_f2(dim=2) = ", a3%f2(2,:,:)

call dot_f1( -1.0, a1 )

print *, "a1_f1(dim=1) = ", a1%f1(1,:)

call pp%del()
call gp%del()

end program test_ufield_arith