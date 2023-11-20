program test_random_number

use math_module

implicit none

integer :: n = 10000000, i
real(kind=8) :: x
real :: t1, t2

call cpu_time(t1)
do i = 1, n
    x = rand_norm()
enddo
call cpu_time(t2)

print *, "new random number generator: ", t2 - t1

call cpu_time(t1)
do i = 1, n
    x = ranorm()
enddo
call cpu_time(t2)

print *, "old random number generator: ", t2 - t1



end program test_random_number