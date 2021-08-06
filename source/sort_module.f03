module sort_module

implicit none

private

public :: generate_sort_idx_1d

contains

subroutine generate_sort_idx_1d( ix, ip, len, range )
! count-sort for integer array

  implicit none
  integer, intent(in), dimension(:) :: ix
  integer, intent(out), dimension(:) :: ip
  integer, intent(in) :: len, range

  integer :: i, index
  integer, dimension(:), allocatable :: counter
  
  allocate( counter(range) )
  counter = 0

  do i = 1, len
    index = ix(i)
    counter(index) = counter(index) + 1
  enddo

  do i = 2, range
    counter(i) = counter(i) + counter(i-1)
  enddo

  do i = 1, len
    index = ix(i)
    ip(i) = counter(index)
    counter(index) = counter(index) - 1
  enddo

  deallocate(counter)
    
end subroutine generate_sort_idx_1d
    
end module sort_module