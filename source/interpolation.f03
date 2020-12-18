module interpolation

implicit none

contains

! ----------------------------------------------------------------------------------------
! Linear
! ----------------------------------------------------------------------------------------
subroutine spline_linear( x, s )
  
  implicit none
  
  real, intent(in) :: x
  real, dimension(0:1), intent(inout) :: s
  
  s(0) = 1.0 - x
  s(1) = x

end subroutine spline_linear

! ----------------------------------------------------------------------------------------
! Quadratic
! ----------------------------------------------------------------------------------------
subroutine spline_quadratic( x, s )
  
  implicit none
  
  real, intent(in) :: x
  real, dimension(-1:1), intent(inout) :: s
  
  real :: t0, t1
  
  t0 = 1.0 - x
  t1 = x
  
  s(-1) = 0.5 * t0**2
  s( 0) = 0.5 + t0*t1
  s( 1) = 0.5 * t1**2

end subroutine spline_quadratic

end module interpolation