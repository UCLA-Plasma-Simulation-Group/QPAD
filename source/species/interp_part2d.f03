module interp_part2d

use ufield_class
use interpolation
use param

implicit none

private

! interface gen_interp_info
!     module procedure gen_interp_info_part2d_s1
!     module procedure gen_interp_info_part2d_s2
! end interface

interface interp_field
    module procedure interp_field_part2d_scalar
    module procedure interp_field_part2d_vector

end interface

interface transform_to_cartesian
    module procedure transform_to_cartesian
end interface

public :: gen_interp_info, interp_field, transform_to_cartesian

contains

subroutine gen_interp_info(x, dr, noff, np, ptrcur, weight, idx, pcos, psin)
    implicit none
    real, dimension(:, :), intent(in) :: x
    real, intent(in) :: dr
    integer, intent(in) :: noff, np
    integer(kind=LG), intent(in) :: ptrcur
    real, dimension(:, :), pointer, intent(inout) :: weight
    integer, dimension(:), intent(out) :: idx
    real, dimension(:), intent(out) :: pcos, psin
    
    integer :: i
    integer(kind=LG) :: pp
    real :: pos, idr
    real, dimension(:), allocatable :: wt
    
    idr = 1.0 / dr
    allocate(wt(size(weight)))
    pp = ptrcur
    if (lbound(weight, 1) == 0) then
        do i = 1, np
    
            pos = sqrt(x(1, pp)**2 + x(2, pp)**2)
    
            ! cosine and sine
            pcos(i) = x(1, pp) / pos
            psin(i) = x(2, pp) / pos
    
            pos = pos * idr
            idx(i) = floor(pos) - noff + 1
            pos = pos - floor(pos)
            call spline_linear(pos, wt)
            weight(0:1, i) = wt
    
            pp = pp + 1
    
        enddo

    else if (lbound(weight, 1) == -1) then
        do i = 1, np
    
            pos = sqrt(x(1, pp)**2 + x(2, pp)**2)
    
            ! cosine and sine
            pcos(i) = x(1, pp) / pos
            psin(i) = x(2, pp) / pos
    
            pos = pos * idr
            idx(i) = floor(pos) - noff + 1
            pos = pos - floor(pos)
            idx(i) = idx(i) + floor(pos - 0.5) + 1
            pos = pos - floor(pos-0.5) - 1
            
            if (noff == 0 .and. (idx(i) == 1)) then
                weight(-1, i) = 0
                weight( 0, i) = 1 - 0.5*(0.5+pos)**2
                weight( 1, i) = 0.5*(0.5+pos)**2
            else
                call spline_quadratic(pos, wt)
                weight(-1:1, i) = wt
            endif
    
            pp = pp + 1
    
        enddo
    else
        print *, 'Invalid lower bound of weight'
    end if
    
end subroutine gen_interp_info


subroutine interp_field_part2d_vector(f_re, f_im, max_mode, weight, idx, pcos, psin, np, fp)

    implicit none
    type(ufield), dimension(:), pointer, intent(in) :: f_re, f_im
    integer, intent(in) :: max_mode, np
    real, dimension(:, :), pointer, intent(inout) :: weight
    integer, dimension(:), intent(in) :: idx
    real, dimension(:), intent(in) :: pcos, psin
    real, dimension(:, :), intent(inout) :: fp

    real, dimension(:, :), pointer :: ptr_f_re => null(), ptr_f_im => null()
    integer :: i, j, mode
    real :: phase_re, phase_im
    complex(kind=DB), dimension(np) :: phase

    fp = 0.0
    phase = cmplx(1.0, 0.0)
    
    if (lbound(weight, 1) == 0) then
        ! interpolate m = 0 mode
        ptr_f_re => f_re(0)%get_f1()
        do i = 1, np
            do j = 0, 1
                fp(:, i) = fp(:, i) + ptr_f_re(:, idx(i) + j) * weight(j, i)
            enddo
        enddo

        ! interpolate m > 0 modes
        do mode = 1, max_mode
            ptr_f_re => f_re(mode)%get_f1()
            ptr_f_im => f_im(mode)%get_f1()

            do i = 1, np
                phase(i) = phase(i) * cmplx(pcos(i), psin(i))
                phase_re = 2.0 * real(phase(i))
                phase_im = 2.0 * aimag(phase(i))
                do j = 0, 1
                    fp(:, i) = fp(:, i) + (ptr_f_re(:, idx(i) + j) * phase_re - &
                                        ptr_f_im(:, idx(i) + j) * phase_im) * weight(j, i)
                enddo
            enddo
        enddo
    else if (lbound(weight, 1) == -1) then
        ! interpolate m = 0 mode
        ptr_f_re => f_re(0)%get_f1()
        do i = 1, np
            do j = -1, 1
                fp(:, i) = fp(:, i) + ptr_f_re(:, idx(i) + j) * weight(j, i)
            enddo
        enddo

        ! interpolate m > 0 modes
        do mode = 1, max_mode
            ptr_f_re => f_re(mode)%get_f1()
            ptr_f_im => f_im(mode)%get_f1()

            do i = 1, np
                phase(i) = phase(i) * cmplx(pcos(i), psin(i))
                phase_re = 2.0 * real(phase(i))
                phase_im = 2.0 * aimag(phase(i))
                do j = -1, 1
                    fp(:, i) = fp(:, i) + (ptr_f_re(:, idx(i) + j) * phase_re - &
                                        ptr_f_im(:, idx(i) + j) * phase_im) * weight(j, i)
                enddo
            enddo
        enddo
    else
        print *, 'Invalid lower bound of weight'
    end if
    
end subroutine interp_field_part2d_vector

subroutine interp_field_part2d_scalar(f_re, f_im, max_mode, weight, idx, pcos, psin, np, fp)

    implicit none
    type(ufield), dimension(:), pointer, intent(in) :: f_re, f_im
    integer, intent(in) :: max_mode, np
    real, dimension(:, :), pointer, intent(inout) :: weight
    integer, dimension(:), intent(in) :: idx
    real, dimension(:), intent(in) :: pcos, psin
    real, dimension(:), intent(inout) :: fp

    real, dimension(:, :), pointer :: ptr_f_re => null(), ptr_f_im => null()
    integer :: i, j, mode
    real :: phase_re, phase_im
    complex(kind=DB), dimension(np) :: phase

    fp = 0.0
    phase = cmplx(1.0, 0.0)

    if (lbound(weight, 1) == 0) then
        ! interpolate m = 0 mode
        ptr_f_re => f_re(0)%get_f1()
        do i = 1, np
            do j = 0, 1
                fp(i) = fp(i) + ptr_f_re(1, idx(i) + j) * weight(j, i)
            enddo
        enddo

        ! interpolate m > 0 modes
        do mode = 1, max_mode
            ptr_f_re => f_re(mode)%get_f1()
            ptr_f_im => f_im(mode)%get_f1()

            do i = 1, np
                phase(i) = phase(i) * cmplx(pcos(i), psin(i))
                phase_re = 2.0 * real(phase(i))
                phase_im = 2.0 * aimag(phase(i))
                do j = 0, 1                                                             
                    fp(i) = fp(i) + (ptr_f_re(1, idx(i) + j) * phase_re - &
                                    ptr_f_im(1, idx(i) + j) * phase_im) * weight(j, i)
                enddo
            enddo
        enddo
    else if (lbound(weight, 1) == -1) then
    ! interpolate m = 0 mode
        ptr_f_re => f_re(0)%get_f1()
        do i = 1, np
            do j = -1, 1
                fp(i) = fp(i) + ptr_f_re(1, idx(i) + j) * weight(j, i)
            enddo
        enddo
    
        ! interpolate m > 0 modes
        do mode = 1, max_mode
            ptr_f_re => f_re(mode)%get_f1()
            ptr_f_im => f_im(mode)%get_f1()
    
            do i = 1, np
                phase(i) = phase(i) * cmplx(pcos(i), psin(i))
                phase_re = 2.0 * real(phase(i))
                phase_im = 2.0 * aimag(phase(i))
                do j = -1, 1                                                             
                    fp(i) = fp(i) + (ptr_f_re(1, idx(i) + j) * phase_re - &
                                     ptr_f_im(1, idx(i) + j) * phase_im) * weight(j, i)
                enddo
            enddo
        enddo
    else
        print *, 'Invalid lower bound of weight'
    end if
    
end subroutine interp_field_part2d_scalar


subroutine transform_to_cartesian(fp, np, pcos, psin)

    implicit none
    real, dimension(:, :), intent(inout) :: fp
    integer, intent(in) :: np
    real, dimension(:), intent(in) :: pcos, psin

    integer :: i
    real :: tmp

    do i = 1, np
        tmp      = fp(1, i) * pcos(i) - fp(2, i) * psin(i)
        fp(2, i) = fp(1, i) * psin(i) + fp(2, i) * pcos(i)
        fp(1, i) = tmp
    enddo

end subroutine transform_to_cartesian

end module interp_part2d