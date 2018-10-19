program test_ufield

use ufield_class

implicit none

type( ufield ) :: a1, a2, a3

integer :: dim = 3
integer, dimension(2) :: nd = (/2,2/), nvp = (/1,1/)
integer, dimension(2,2) :: gc_num
real :: k = 5.0

real :: start, finish
integer :: i

gc_num = 0

call a1%new( dim, nd, nvp, gc_num )
call a2%new( dim, nd, nvp, gc_num )
call a3%new( dim, nd, nvp, gc_num )

a1%f1 = 1.0
a2%f1 = 2.0

print *, "a1 = ", a1%f1
print *, "a2 = ", a2%f1


! call cpu_time( start )
! do i = 1, 100000
!     a3 = a1 + a2
! enddo
! call cpu_time( finish )
! print *, '("Time = ",f6.3," seconds.")', finish-start

a3 = k*(a1+a2)+(a1+k)*a2

print *, "a3 = ", a3%f1


end program test_ufield