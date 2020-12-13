module debug_tool

use, intrinsic :: ieee_arithmetic

implicit none

public

interface write_data
  module procedure write_data1d
  module procedure write_data2d
  module procedure write_array
end interface write_data

interface check_array
  module procedure check_array_rank1
  module procedure check_array_rank2
end interface check_array

public :: write_data

contains

subroutine write_array( f, fname )

  implicit none

  real, intent(in), dimension(:), pointer :: f
  character(len=*), intent(in) :: fname

  integer :: i1, i2, i
  integer :: unit = 99

  i1 = lbound(f,1)
  i2 = ubound(f,1)

  open( unit, file=trim(fname) )

  do i = i1, i2
    write( unit, '(ES20.6)' ) f(i)
  enddo

  close( unit )

end subroutine write_array

subroutine write_data1d( f, fname, dim )

  implicit none

  real, intent(in), dimension(:,:), pointer :: f
  character(len=*), intent(in) :: fname
  integer, intent(in) :: dim

  integer :: i1, i2, i
  integer :: unit = 99

  i1 = lbound(f,2)
  i2 = ubound(f,2)

  open( unit, file=trim(fname) )

  do i = i1, i2
    write( unit, '(ES20.6)' ) f(dim,i)
  enddo

  close( unit )

end subroutine write_data1d

subroutine write_data2d( f, fname, dim )

  implicit none

  real, intent(in), dimension(:,:,:), pointer :: f
  character(len=*), intent(in) :: fname
  integer, intent(in) :: dim

  integer :: i1, i2, i
  integer :: unit = 2
  character(len=128) :: fmtstr, num2str

  i1 = lbound(f,2)
  i2 = ubound(f,2)

  write( num2str, * ) size(f,3)
  write( fmtstr, * ) '(' // trim(num2str) // 'ES20.6)'
  open( unit, file=trim(fname) )

  do i = i1, i2
    write( unit, trim(fmtstr) ) f(dim,i,:)
  enddo

  close( unit )

end subroutine write_data2d

subroutine check_array_rank1( array, th, stat )

  implicit none

  real, intent(in), dimension(:) :: array
  real, intent(in) :: th
  integer, intent(out) :: stat

  integer :: length, i
  real :: infinity = 1.0d100

  length = size( array )

  stat = 0
  do i = 1, length
    if ( ieee_is_nan( array(i) ) .or. abs( array(i) ) > th ) then
      stat = -1
      print *, 'invalid value at index ', i
      return
    endif
  enddo

end subroutine check_array_rank1

subroutine check_array_rank2( array, th, stat )

  implicit none

  real, intent(in), dimension(:,:) :: array
  real, intent(in) :: th
  integer, intent(out) :: stat

  integer :: n1, n2, i, j
  real :: infinity = 1.0d100

  n1 = size( array, 1 )
  n2 = size( array, 2 )

  stat = 0
  do j = 1, n2
    do i = 1, n1
      if ( ieee_is_nan( array(i,j) ) .or. abs( array(i,j) ) > th ) then
        stat = -1
        print *, 'invalid value at index ', i, j
        return
      endif
    enddo
  enddo

end subroutine check_array_rank2

end module debug_tool