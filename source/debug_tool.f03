module debug_tool

implicit none

public

interface write_data
  module procedure write_data1d
  module procedure write_data2d
  module procedure write_array
end interface write_data

public :: write_data

contains

subroutine write_array( f, fname )

  implicit none

  real, intent(in), dimension(:), pointer :: f
  character(len=*), intent(in) :: fname

  ! integer :: i1, i2, i
  integer :: unit = 99
  character(len=128) :: fmtstr, num2str

  ! i1 = lbound(f,2)
  ! i2 = ubound(f,2)

  write( num2str, * ) size(f)
  write( fmtstr, * ) '(' // trim(num2str) // 'f20.8)'
  open( unit, file=trim(fname) )

  ! do i = i1, i2
  write( unit, trim(fmtstr) ) f(:)
  ! enddo

  close( unit )

end subroutine write_array

subroutine write_data1d( f, fname, dim )

  implicit none

  real, intent(in), dimension(:,:), pointer :: f
  character(len=*), intent(in) :: fname
  integer, intent(in) :: dim

  integer :: i1, i2, i
  integer :: unit = 99
  character(len=128) :: fmtstr, num2str

  i1 = lbound(f,2)
  i2 = ubound(f,2)

  ! write( num2str, * ) size(f,2)
  ! write( fmtstr, * ) '(' // trim(num2str) // 'f20.8)'
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

end module debug_tool