module math_module

use param

implicit none

contains

! this program calculates a random number y from a gaussian distribution
! with zero mean and unit variance, according to the method of
! mueller and box:
!    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
!    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
! where x is a random number uniformly distributed on (0,1).
! written for the ibm by viktor k. decyk, ucla
function ranorm()
    
integer, save :: r1 = 885098780, r2 = 1824280461
integer, save :: iflg = 0, r4 = 1396483093, r5 = 55318673
real(kind=DB), save :: h1l = 65531.0d0, h1u = 32767.0d0
real(kind=DB), save :: h2l = 65525.0d0,r0 = 0.0d0
real(kind=DB) :: ranorm, r3, asc, bsc, temp
integer :: isc, i1

if (iflg == 0) then
  isc = 65536
  asc = dble(isc)
  bsc = asc*asc
  i1 = r1 - (r1/isc)*isc
  r3 = h1l*dble(r1) + asc*h1u*dble(i1)
  i1 = r3/bsc
  r3 = r3 - dble(i1)*bsc
  bsc = 0.5d0*bsc
  i1 = r2/isc
  isc = r2 - i1*isc
  r0 = h1l*dble(r2) + asc*h1u*dble(isc)
  asc = 1.0d0/bsc
  isc = r0*asc
  r2 = r0 - dble(isc)*bsc
  r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
  isc = r3*asc
  r1 = r3 - dble(isc)*bsc
  temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
  isc = 65536
  asc = dble(isc)
  bsc = asc*asc
  i1 = r4 - (r4/isc)*isc
  r3 = h2l*dble(r4) + asc*h1u*dble(i1)
  i1 = r3/bsc
  r3 = r3 - dble(i1)*bsc
  bsc = 0.5d0*bsc
  i1 = r5/isc
  isc = r5 - i1*isc
  r0 = h2l*dble(r5) + asc*h1u*dble(isc)
  asc = 1.0d0/bsc
  isc = r0*asc
  r5 = r0 - dble(isc)*bsc
  r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
  isc = r3*asc
  r4 = r3 - dble(isc)*bsc
  r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
  ranorm = temp*dsin(r0)
  r0 = temp*dcos(r0)
  iflg = 1
  return
else
  ranorm = r0
  r0 = 0.0d0
  iflg = 0
  return
endif

end function ranorm

function eval_polynomial( x, x0, p ) result(y)
  
  implicit none
  real, intent(in) :: x, x0
  real, dimension(:) :: p
  real :: y

  integer :: i, len_p

  len_p = size(p)
  y = 0.0
  do i = 1, len_p
    y = y + p(i) * (x - x0)**(i-1)
  enddo

end function eval_polynomial

function factorial( n ) result(res)

  implicit none
  integer, intent(in) :: n
  integer :: res, i

  if ( n < 0 ) then
    write(2, *) '[ERROR] Input of factorial function must be positive/zero!'
    write(*, *) 'Input of factorial function must be positive/zero!'
    stop
  endif
  res = 1
  do i = 1, n
    res = res * i
  enddo

end function factorial

! Generalized Laguerre polynomial (recurrence relation)
! https://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials
function laguerre( n, alpha, x ) result(y)

  implicit none
  integer, intent(in) :: n, alpha
  real, intent(in) :: x
  real :: y, l0, l1
  integer :: k

  l0 = 1.0
  l1 = 1.0 + real(alpha) - x

  if ( n == 0 ) then
    y = l0
  elseif ( n == 1 ) then
    y = l1
  elseif ( n > 1 ) then
    do k = 1, n-1
      y = ( ( real( 2*k + 1 + alpha ) - x ) * l1 - real( k + alpha ) * l0 ) / real( k + 1 )
      l0 = l1
      l1 = y
    enddo
  else
    write(2, *) '[ERROR] The first input of generalized Laguerre polynomial must be positive/zero!'
    write(*, *) 'The first input of generalized Laguerre polynomial must be positive/zero!'
    stop
  endif

end function laguerre

end module math_module