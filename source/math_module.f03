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

! A straightforward implementation (Box-Muller method) of random number generator
! for Gaussian distribution with zero mean and unit variance.
function rand_norm()

  implicit none
  logical, save :: iflg = .true.
  real(kind=DB), save :: next_val = 0.0d0
  real(kind=DB) :: u1, u2, theta, amp, rand_norm

  if (iflg) then

    call random_number(u1)
    call random_number(u2)
    theta = 2.0d0 * pi * u1
    amp = dsqrt(-2.0d0 * dlog(u2))
    rand_norm = amp * dcos(theta)
    next_val = amp * dsin(theta)
    iflg = .false.

  else

    rand_norm = next_val
    next_val = 0.0d0
    iflg = .true.

  endif

end function rand_norm

function rand_super_gauss(p)

  implicit none
  real, intent(in) :: p
  real(kind=DB) :: tmp, fac, x, u, sqrt_pi, accept, gam, rand_super_gauss

  if (p < 1.0) then
    write(2, *) '[ERROR] The order p for super-Gaussian distribution must be greater than 1.0'
    stop
  endif

  sqrt_pi = sqrt(pi)
  tmp = 1.0d0 / (1.0d0 - p)
  gam = gamma(1.0d0 + 1.0d0 / (2.0d0 * p))
  fac = 2.0d0 * exp(p ** (p * tmp) - p ** tmp) * gam / sqrt_pi;

  do
    x = rand_norm()
    call random_number(u)

    tmp = 0.5 * x * x
    accept = 0.5 * fac * exp(tmp - tmp ** p) * sqrt_pi / gam
    if (u < accept) then
      rand_super_gauss = x
      return
    endif
  enddo
  
end function rand_super_gauss

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
! function laguerre( n, alpha, x ) result(y)

!   implicit none
!   integer, intent(in) :: n, alpha
!   real, intent(in) :: x
!   real :: y, l0, l1
!   integer :: k

!   l0 = 1.0
!   l1 = 1.0 + real(alpha) - x

!   if ( n == 0 ) then
!     y = l0
!   elseif ( n == 1 ) then
!     y = l1
!   elseif ( n > 1 ) then
!     do k = 1, n-1
!       y = ( ( real( 2*k + 1 + alpha ) - x ) * l1 - real( k + alpha ) * l0 ) / real( k + 1 )
!       l0 = l1
!       l1 = y
!     enddo
!   else
!     write(2, *) '[ERROR] The first input of generalized Laguerre polynomial must be positive/zero!'
!     write(*, *) 'The first input of generalized Laguerre polynomial must be positive/zero!'
!     stop
!   endif

! end function laguerre

function laguerre( p, l, u )

  implicit none

  real :: laguerre
  integer, intent(in) :: p, l
  real, intent(in) :: u

  select case (p)
  case(0)
    laguerre = 1.0
  case(1)
    laguerre = 1.0 + l - u
  case(2)
    laguerre = (2 + 3.0* l + l**2 - 4.0* u - &
        2.0* l*u + u**2)/2.0
  case(3)
    laguerre = (6.0 + 11.0* l + 6.0* l**2 + l**3 - &
          18.0* u - 15.0* l*u - 3.0* l**2*u + &
          9.0* u**2 + 3.0* l*u**2 - u**3)/6.0
  case(4)
    laguerre = (24.0 + 50.0* l + 35.0* l**2 + &
          10.0* l**3 + l**4 - 96.0* u - 104.0* l*u - &
          36.0* l**2*u - 4.0* l**3*u + 72.0* u**2 + &
          42.0* l*u**2 + 6.0* l**2*u**2 - 16.0* u**3 - &
          4.0* l*u**3 + u**4)/24.0
  case(5)
    laguerre = (120.0 + 274.0* l + 225.0* l**2 + &
          85.0* l**3 + 15.0* l**4 + l**5 - 600.0* u - &
          770.0* l*u - 355.0* l**2*u - 70.0* l**3*u - &
          5.0* l**4*u + 600.0* u**2 + 470.0* l*u**2 + &
          120.0* l**2*u**2 + 10.0* l**3*u**2 - &
          200.0* u**3 - 90.0* l*u**3 - 10.0* l**2*u**3 + &
          25.0* u**4 + 5.0* l*u**4 - u**5)/120.0
  case(6)
    laguerre = (720.0 + 1764.0* l + 1624.0* l**2 + &
          735.0* l**3 + 175.0* l**4 + 21.0* l**5 + &
          l**6 - 4320.0* u - 6264.0* l*u - &
          3480.0* l**2*u - 930.0* l**3*u - &
          120.0* l**4*u - 6.0* l**5*u + 5400.0* u**2 + &
          5130.0* l*u**2 + 1785.0* l**2*u**2 + &
          270.0* l**3*u**2 + 15.0* l**4*u**2 - &
          2400.0* u**3 - 1480.0* l*u**3 - &
          300.0* l**2*u**3 - 20.0* l**3*u**3 + &
          450.0* u**4 + 165.0* l*u**4 + 15.0* l**2*u**4 - &
          36.0* u**5 - 6.0* l*u**5 + u**6)/720.0
  case(7)
    laguerre = (5040.0 + 13068.0* l + 13132.0* l**2 + &
          6769.0* l**3 + 1960.0* l**4 + 322.0* l**5 + &
          28.0* l**6 + l**7 - 35280.0* u - &
          56196.0* l*u - 35728.0* l**2*u - &
          11655.0* l**3*u - 2065.0* l**4*u - &
          189.0* l**5*u - 7.0* l**6*u + 52920.0* u**2 + &
          57834.0* l*u**2 + 24675.0* l**2*u**2 + &
          5145.0* l**3*u**2 + 525.0* l**4*u**2 + &
          21.0* l**5*u**2 - 29400.0* u**3 - &
          22330.0* l*u**3 - 6265.0* l**2*u**3 - &
          770.0* l**3*u**3 - 35.0* l**4*u**3 + &
          7350.0* u**4 + 3745.0* l*u**4 + &
          630.0* l**2*u**4 + 35.0* l**3*u**4 - &
          882.0* u**5 - 273.0* l*u**5 - 21.0* l**2*u**5 + &
          49.0* u**6 + 7.0* l*u**6 - u**7)/5040.0
  case default
    laguerre = 0
  end select

end function laguerre

! solve tri-diagonal matrix equation using Gaussian elimination
! https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
subroutine tdma(n, a, b, c, d, x)
  implicit none
  integer, intent(in) :: n
  real, intent(in), dimension(n) :: a, c
  real, intent(inout), dimension(n) :: b, d
  real, intent(out), dimension(n) :: x
  !  --- Local variables ---
  integer :: i
  real :: q
  !  --- Elimination ---
  do i = 2, n
    q = a(i) / b(i - 1)
    b(i) = b(i) - c(i - 1) * q
    d(i) = d(i) - d(i - 1) * q
  enddo
  ! --- Backsubstitution ---
  q = d(n) / b(n)
  x(n) = q
  do i = n - 1, 1, -1
    q = (d(i) - c(i) * q) / b(i)
    x(i) = q
  enddo
end subroutine tdma

end module math_module