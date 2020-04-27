module interpolation

implicit none

contains

! !-----------------------------------------------------------------------------------------
! subroutine get_emf_1d( bp, ep, b, e, ix, x, np )
! !-----------------------------------------------------------------------------------------
! !  calculates the values of the electric-magnetic field 
! !  at the positions of the array x on a 1d grid, using ORDER spline 
! !  interpolation
! !---------------------------------------------------

!   implicit none

!   integer, parameter :: rank = 1

!   real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

!   type( t_vdf ), intent(in) :: b, e

!   real(p_k_part), dimension(:,:), intent(in) :: x
!   integer, dimension(:,:), intent(in) :: ix

!   integer, intent(in) :: np


!   real(p_k_fld) :: dx1, f1, f2, f3
!   real(p_k_fld), dimension(LP:UP) :: w1
!   integer :: l, i1, k1


!   do l = 1, np
     
!      i1 = ix(1,l)
!      dx1 = real( x(1,l), p_k_fld )
     
!      ! get spline weitghts for x
!      call SPLINE( dx1, w1 )
     
!      ! Interpolate Fields
!      f1 = 0
!      f2 = 0
!      f3 = 0
     
!      do k1 = LP, UP
!        f1 = f1 + e%f1(1,i1 + k1) * w1(k1)
!        f2 = f2 + e%f1(2,i1 + k1) * w1(k1)
!        f3 = f3 + e%f1(3,i1 + k1) * w1(k1)
!      enddo
     
!      ep( 1, l ) = f1
!      ep( 2, l ) = f2
!      ep( 3, l ) = f3

!      f1 = 0
!      f2 = 0
!      f3 = 0
     
!      do k1 = LP, UP
!        f1 = f1 + b%f1(1,i1 + k1) * w1(k1)
!        f2 = f2 + b%f1(2,i1 + k1) * w1(k1)
!        f3 = f3 + b%f1(3,i1 + k1) * w1(k1)
!      enddo
     
!      bp( 1, l ) = f1
!      bp( 2, l ) = f2
!      bp( 3, l ) = f3
     
!   enddo

! end subroutine get_emf_1d
! !-----------------------------------------------------------------------------------------

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