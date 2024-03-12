program test_field_arith

use parallel_pipe_class
use grid_class
use field_class
! use ufield_class

implicit none

type( parallel_pipe ), pointer :: pp => null()
type( grid ), pointer :: gp => null()
type( field ) :: a1, a2, a3
! type( ufield ), pointer :: ua1 => null(), ua2 => null(), ua3 => null()

integer :: nr = 2, nz = 2
integer, dimension(2,2) :: gc_num = 0
real :: k = 5.0

real :: start, finish
integer :: i

allocate( pp, gp )
call pp%new(nst=1)
call gp%new( pp, nr, nz )


call a1%new( pp, gp, dim=3, dr=1.0, dxi=1.0, num_modes=1, gc_num=gc_num )
call a2%new( pp, gp, dim=3, dr=1.0, dxi=1.0, num_modes=1, gc_num=gc_num )
call a3%new( pp, gp, dim=3, dr=1.0, dxi=1.0, num_modes=1, gc_num=gc_num )

a1 = 0.5
a2 = 0.1

call a1%as(1.0)
call a2%as(2.0)

print *, "re(a1_f1)(dim=1, mode=0) = ", a1%rf_re(0)%f1(1,:)
print *, "re(a2_f1)(dim=1, mode=0) = ", a2%rf_re(0)%f1(1,:)

print *, "re(a1_f2)(dim=1, mode=0) = ", a1%rf_re(0)%f2(1,:,:)
print *, "re(a2_f2)(dim=1, mode=0) = ", a2%rf_re(0)%f2(1,:,:)

! arithmetic on f1
a3 = k * a1 + a2 - a2 * k
! ! arithmetic on f2
call a3%as( k .dot. a1 .add. a2 .sub. (a2 .dot. k)  )

print *, "re(a3_f1)(dim=1, mode=0) = ", a3%rf_re(0)%f1(1,:)
print *, "re(a3_f2)(dim=1, mode=0) = ", a3%rf_re(0)%f2(1,:,:)

call pp%del()
call gp%del()

end program test_field_arith